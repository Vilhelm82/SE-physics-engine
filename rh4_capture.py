"""RH-4 -- the register at the absorbing interface: a coherent CAPTURE in the hold, built on FD's own instrument.

Everything from FD is kept: loop_endpoints_with_effects (finite ramps, inverses, stretches, detuning, background
absorption gamma integrated inside the no-jump propagators, exact lost-effect operators), the absorbing hold
diag(e^{-kappa t_d/2}, e^{-kappa t_d/2}, e^{-i delta_d t_d}, e^{-i delta_r t_d}) with kappa t_d = 40, monitored (heralded).
ONE thing is added.  Before each absorbing hold, a capture pulse of duration t_c exchanges the bright pair Q = (p, q)
with a fresh register pair R_k = (r_kp, r_kq):  H_cap = g sum_m (|m><r_km| + h.c.),  angle theta = g t_c = (pi/2)(1 + e_c),
e_c a calibration error unknown to the recovery.  Background gamma acts on Q and detuning phases on P during t_c.
Register modes decohere from their fill time to readout: loss gamma_r (no-jump; the lost population is HERALDED if the
register is monitored, else replaced by I/2) and dephasing gamma_phi (coherences between distinct register modes, and
register-P coherences, decay at gamma_phi per filled mode).  Whatever the capture misses is absorbed and heralded
exactly as in FD: the fallback is FD.
Readout after the last hold: projective 'register occupied'.  Register branch -> U_cal^dag (10 -> 2) then G, with U_cal
the polar isometry of the stacked register map of the IDEAL-register model (gamma_r = gamma_phi = e_c = 0) at the
operating point; register amplitude outside range(U_cal) is heralded leftover.  Residual Q -> I/2 (FD's convention).
Output: herald h and the exact Pauli twirl of the accepted map (from its Pauli transfer matrix), in qec_distance_stack's
channel format.  Regression: capture off must reproduce native_channel(row) exactly.
Space: single-excitation P + Q + R_1..R_5 = 14 dims; density matrices propagated stage by stage.
Written 2026-09-06 by Claire.  Not a hardware claim: t_c, gamma_r, gamma_phi, e_c are declared model parameters.
"""
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import expm, polar
import reflection_loop_finite_dump as FD
import qec_distance_stack as Q

N = 5; DIM = 4 + 2*N
CASE = (0.01, 1e-4, 0.0); GAMMA = 1e-6; KAPPA, ACTION = 10., 40.
PAULI2 = Q.PAULIS
CHAR = np.array([[1, 1, 1, 1], [1, 1, -1, -1], [1, -1, 1, -1], [1, -1, -1, 1]])   # chi(P, a): +1 commute, -1 anticommute
IP, IQ = [2, 3], [0, 1]
IR = lambda k: [4 + 2*k, 5 + 2*k]
ALLR = [m for k in range(N) for m in IR(k)]


def embed4(m):
    out = np.eye(DIM, dtype=complex); out[:4, :4] = m; return out


def embed4z(m):
    """Zero-padded embedding for EFFECT operators (an effect on P+Q must not act on the register)."""
    out = np.zeros((DIM, DIM), dtype=complex); out[:4, :4] = m; return out


class Model:
    """Steps are no-jump matrices M (rho -> M rho M^dag); dephasing is a separate elementwise factor."""
    def __init__(self, e_c=0., t_c=1., gamma_r=0., gamma_phi=0., register_heralded=True, capture=True):
        self.th = (np.pi/2)*(1 + e_c); self.t_c = t_c if capture else 0.; self.capture = capture
        self.gr, self.gphi, self.rh = gamma_r, gamma_phi, register_heralded
        loops, effects = FD.loop_endpoints_with_effects([CASE], [GAMMA])
        self.loops, self.effects = loops[0], effects[0]
        stages = FD.primitive_stages(1); T0 = sum(s.duration for s in stages)
        _, _, stretches = FD.controls('five'); self.loop_T = [T0*float(s) for s in stretches]
        self.td = ACTION/KAPPA; self.chi = np.exp(-ACTION/2); self.delta = np.array([CASE[1], CASE[2]])
        self.T_total = sum(self.loop_T) + N*(self.td + self.t_c)

    # ---- matrices
    def M_phases_P(self, t):
        d = np.ones(DIM, dtype=complex); d[IP] = np.exp(-1j*self.delta*t); return np.diag(d)
    def M_regdecay(self, t, filled):
        d = np.ones(DIM); d[[m for k in filled for m in IR(k)]] = np.exp(-self.gr*t/2); return np.diag(d)
    def M_loop(self, k): return embed4(self.loops[k])
    def M_capture(self, k):
        M = np.eye(DIM, dtype=complex)
        if not self.capture or self.t_c == 0: return M
        gen = np.array([[-GAMMA/2, -1j*self.th/self.t_c], [-1j*self.th/self.t_c, -self.gr/2]])
        blk = expm(gen*self.t_c)
        for m, r in zip(IQ, IR(k)): M[np.ix_([m, r], [m, r])] = blk
        return M@self.M_phases_P(self.t_c)
    def M_absorb(self):
        d = np.ones(DIM, dtype=complex); d[IQ] = self.chi; return np.diag(d)@self.M_phases_P(self.td)
    def dephase(self, rho, t, filled):
        if not self.gphi or not filled: return rho
        f = np.zeros(DIM); f[[m for k in filled for m in IR(k)]] = 1.
        decay = np.exp(-self.gphi*t*(f[:, None] + f[None, :])); np.fill_diagonal(decay, 1.)
        return rho*decay

    # ---- pure linear map (valid for gamma_phi = 0): total no-jump matrix, chronological
    def total_matrix(self):
        M = np.eye(DIM, dtype=complex); filled = []
        for k in range(N):
            M = self.M_regdecay(self.loop_T[k], filled)@self.M_loop(k)@M
            M = self.M_regdecay(self.t_c, filled)@self.M_capture(k)@M; filled.append(k)
            M = self.M_regdecay(self.td, filled)@self.M_absorb()@M
        return M

    # ---- density-matrix propagation with loss attribution
    def run(self, rho_P):
        rho = np.zeros((DIM, DIM), dtype=complex); rho[np.ix_(IP, IP)] = rho_P
        her = unh = 0.; filled = []
        def apply(M, rho, heralded):
            new = M@rho@M.conj().T; lost = float(np.trace(rho).real - np.trace(new).real)
            return new, (lost if heralded else 0.), (0. if heralded else lost)
        for k in range(N):
            unh += float(np.trace(embed4z(self.effects[k])@rho).real)                     # background during transport (exact effect)
            rho = self.M_loop(k)@rho@self.M_loop(k).conj().T
            rho, a, b = apply(self.M_regdecay(self.loop_T[k], filled), rho, self.rh); her += a; unh += b
            rho = self.dephase(rho, self.loop_T[k], filled)
            rho, a, b = apply(self.M_capture(k), rho, False); her += a; unh += b            # Q background + register loss during t_c: unheralded
            rho, a, b = apply(self.M_regdecay(self.t_c, filled), rho, self.rh); her += a; unh += b
            rho = self.dephase(rho, self.t_c, filled); filled = filled + [k]
            rho, a, b = apply(self.M_absorb(), rho, True); her += a; unh += b               # monitored absorption: heralded (FD)
            rho, a, b = apply(self.M_regdecay(self.td, filled), rho, self.rh); her += a; unh += b
            rho = self.dephase(rho, self.td, filled)
        return rho, her, unh


def channel(model, U_cal):
    """Herald and exact Pauli twirl of the accepted map, qec_distance_stack format."""
    Gm = Q.G; ptm = np.zeros(4); heralds = []
    for a, sig in enumerate(PAULI2):
        rho_in = sig/2 if a == 0 else sig/2 + np.eye(2)/2          # trace-1 inputs: (I + sig)/2, and I/2 for a = 0
        rho, her, unh = model.run(rho_in)
        rP = rho[np.ix_(IP, IP)]; rR = rho[np.ix_(ALLR, ALLR)]; rQ = float(np.trace(rho[np.ix_(IQ, IQ)]).real)
        rec = U_cal.conj().T@rR@U_cal; leftover = float(np.trace(rR).real - np.trace(rec).real)
        out = rP + Gm@rec@Gm.conj().T + (unh + rQ)*np.eye(2)/2      # accepted map: survivor + recovered + unheralded replaced
        her_tot = her + leftover; heralds.append(her_tot)
        # PTM diagonal relative to the target G: R_aa = tr(sig_a G^dag out G)/tr(sig_a rho_in) ... use linearity on (I+sig)/2 - I/2
        if a == 0: base = Gm.conj().T@out@Gm; ptm[0] = float(np.trace(base).real)
        else: ptm[a] = float(np.trace(sig@(Gm.conj().T@out@Gm - base)).real)
    h = float(heralds[0]); spread = float(np.ptp(heralds))                  # twirl folds state-dependence of the flag into its average
    w = CHAR@ptm/4; assert abs(w.sum() + h - 1) < 1e-9, (w.sum(), h)
    cond = w/(1 - h); cond[0] = 1 - cond[1:].sum()
    return dict(herald=h, herald_state_spread=spread, pauli=cond.tolist(), unheralded_pauli_rate=float(w[1:].sum()), total_time=model.T_total)


def calibrate(t_c):
    ideal = Model(e_c=0., t_c=t_c, gamma_r=0., gamma_phi=0.)
    F = ideal.total_matrix()[np.ix_(ALLR, IP)]; U, H = polar(F)
    s = np.linalg.svd(F, compute_uv=False)
    return U, dict(captured_fraction=float((s**2).sum()/2), DB_stack=float((s[0]**2 - s[1]**2)/(s**2).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--t-c', type=float, default=1.0)
    ap.add_argument('--qec', action='store_true', help='run d=3,5 surface codes on selected points')
    ap.add_argument('--shots', type=int, default=200000)
    ap.add_argument('--json', type=Path, default=Q.ROOT/'docs/rh4_capture-checks.json')
    args = ap.parse_args(); t0 = time.time(); checks = []
    def check(name, value, tol):
        ok = bool(np.isfinite(value) and abs(value) <= tol); checks.append(dict(name=name, value=float(value), tolerance=tol, passed=ok))
        print(f"  [{'ok' if ok else 'FAIL'}] {name}: {value:.3e} (tol {tol:.0e})", flush=True)
    row = Q.operating_point(background=GAMMA); fd = Q.native_channel(row)
    print("== regression: capture off must be FD ==", flush=True)
    off = Model(capture=False); U0 = np.zeros((2*N, 2), dtype=complex)
    ch_off = channel(off, U0)
    check("capture off: herald == FD herald", ch_off['herald'] - fd['herald'], 1e-9)
    check("capture off: unheralded Pauli == FD (4.45092e-6)", ch_off['unheralded_pauli_rate'] - sum(fd['pauli'][1:])*(1 - fd['herald']), 1e-11)
    check("capture off: conditional Pauli vector == FD", float(np.abs(np.array(ch_off['pauli']) - np.array(fd['pauli'])).max()), 1e-9)
    print("\n== calibration (ideal register, e_c = 0) ==", flush=True)
    U, cal = calibrate(args.t_c)
    print(f"  captured fraction of the input = {cal['captured_fraction']:.6e} (FD herald {fd['herald']:.6e});  DB_stack of the register map = {cal['DB_stack']:.3e}", flush=True)
    check("calibration: captured fraction == FD herald (capture takes what FD absorbed)", cal['captured_fraction'] - fd['herald'], 2e-6)
    print(f"\n== sweep (all at n=1, gain 1%, delta_d/a=1e-4, gamma/a=1e-6, t_c = {args.t_c:.2f}; T = {Model(t_c=args.t_c).T_total:.2f}) ==", flush=True)
    print(f"  {'e_c':>6} {'gamma_r':>8} {'gamma_phi':>9} {'reg.her':>7} | {'herald':>11} {'unher. Pauli':>13} {'r_eff':>8}   note", flush=True)
    grid = []
    base_unh = sum(fd['pauli'][1:])*(1 - fd['herald'])
    def r_eff(unh): return (unh - base_unh)/(0.75*fd['herald'])     # fraction of the captured shots replaced by I/2, RH-3's r
    for e_c in (0., 0.01, 0.05):
        for gr in (0., 1e-4, 1e-3, 1e-2):
            for gphi in (0., 1e-4, 1e-3, 1e-2):
                for rh in (True, False):
                    if (gr == 0. and not rh): continue
                    if (e_c and (gr or gphi)) and not (gr in (0., 1e-3) and gphi in (0., 1e-3)): continue   # keep the mixed corner small
                    m = Model(e_c=e_c, t_c=args.t_c, gamma_r=gr, gamma_phi=gphi, register_heralded=rh)
                    ch = channel(m, U); ch.update(e_c=e_c, gamma_r=gr, gamma_phi=gphi, register_heralded=rh, r_eff=r_eff(ch['unheralded_pauli_rate']))
                    grid.append(ch)
                    note = 'ideal' if not (e_c or gr or gphi) else ''
                    print(f"  {e_c:6.3f} {gr:8.0e} {gphi:9.0e} {str(rh):>7} | {ch['herald']:11.4e} {ch['unheralded_pauli_rate']:13.4e} {ch['r_eff']:8.4f}   {note}", flush=True)
    ideal = next(g for g in grid if not (g['e_c'] or g['gamma_r'] or g['gamma_phi']))
    check("ideal capture: herald < 1e-2 x FD herald", ideal['herald']/fd['herald'], 1e-2)
    check("ideal capture: unheralded Pauli within 10% of FD's background floor", ideal['unheralded_pauli_rate']/base_unh - 1, 0.1)
    res = dict(scope='RH-4 capture model at the absorbing interface; FD instrument with a coherent capture before each absorbing hold',
               operating_point=dict(case=CASE, gamma=GAMMA, kappa=KAPPA, action=ACTION, t_c=args.t_c), fd_channel=fd, capture_off=ch_off,
               calibration=cal, grid=grid, checks=checks, qec=[])
    if args.qec:
        print("\n== decoder: d=3,5 surface, both bases, selected points ==", flush=True)
        sel = [g for g in grid if (g['e_c'], g['gamma_r'], g['gamma_phi'], g['register_heralded']) in
               [(0., 0., 0., True), (0., 1e-3, 0., False), (0., 0., 1e-2, True), (0., 1e-2, 0., False), (0., 1e-2, 1e-2, False), (0., 1e-2, 1e-2, True), (0.05, 1e-3, 1e-3, False)]]
        codes = [Q.surface_circuit(d, b, 1e-3) for d in (3, 5) for b in ('X', 'Z')]
        runs = 0
        for name, chn in [('fd', fd)] + [(f"cap e_c={g['e_c']:g} gr={g['gamma_r']:g} gphi={g['gamma_phi']:g} her={g['register_heralded']}", g) for g in sel]:
            fails = []
            for code in codes:
                r = Q.run_case(code, chn, args.shots, 20260906 + runs*1009); runs += 1
                fails.append((code.name[-2:] + code.basis, r['flags_visible']['failures'])); res['qec'].append(dict(scenario=name, code=code.name, basis=code.basis, **{k: r[k] for k in ('shots', 'seed', 'flags_visible', 'flags_hidden')}))
            print(f"  {name:48s} " + '  '.join(f"{c}={f}" for c, f in fails), flush=True)
            args.json.write_text(json.dumps(res, indent=2) + '\n')
    res['check_count'] = len(checks); res['checks_passed'] = int(sum(c['passed'] for c in checks))
    res['provenance'] = dict(date='2026-09-06', runtime_s=round(time.time() - t0, 1), head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                             sha256={p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in ('rh4_capture.py', 'reflection_loop_finite_dump.py', 'qec_distance_stack.py')})
    args.json.write_text(json.dumps(res, indent=2) + '\n')
    print(f"\nRESULT: {res['checks_passed']}/{res['check_count']} checks passed.  runtime {res['provenance']['runtime_s']}s", flush=True)
    if res['checks_passed'] != res['check_count']: sys.exit(1)


if __name__ == '__main__':
    main()
