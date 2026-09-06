"""RH-2P -- precision-scaling separation test on the eps^3 law of DB_stack (tight-frame trine, gain-only).

Applies Cella U-0784 (thm:precision-separation, Theorem 2 of transfer_function_exponent_family_v6.tex):
    C_{p,k} = a + u_p b_k
with a the path-invariant GEOMETRIC amplitude, b_k the path-dependent CONDITIONING amplitude, u_p the unit
roundoff at precision p.  Also U-0786 (log-log transfer slope) and U-0799 (two-term finite-window bias).
Object:  C(eps) := DB_stack(eps)/eps^3, DB_stack the tight-frame deviation of the stacked leak (RH-2).
Algebra: the SAME closed-form loop as rh1_common.trine_loop at gain-only (segment/frame/rotation), re-implemented
over two numeric backends so precision is a parameter:  numpy complex64 (u=2^-24), numpy complex128 (2^-53),
mpmath 30 dps, mpmath 50 dps.  float64 is checked against rh1_common.trine_loop to 1e-13 before anything else.
Paths (same real quantity, different arithmetic routes):
  k1  singular values of the stacked leak F (6x2)                                   [svd]
  k2  eigenvalues of F^dag F by the stable 2x2 formula  sqrt(((a-d)/2)^2 + |b|^2) / ((a+d)/2)
  k3  eigenvalues of F^dag F by  sqrt(tr^2 - 4 det) / tr                            [cancellation in the radicand]
  k4  eigenvalues of I - K^dag K by the stable 2x2 formula                          [cancellation forming the matrix]
PREDICTIONS (Claire, written before the run, commit precedes log):
  (P1) a(eps) is finite, nonzero, path-invariant across k1..k4 and eps-independent: the cubic law is geometric.
  (P2) b_k3, b_k4 >> b_k2 ~ b_k1; k4 fails outright in float64 (relative error > 0.1) at eps <= 1.25e-3.
  (P3) local log-log slope at 50 dps -> 3.000 from above or below with a two-term correction of exponent rho ~ 1.
  (P4) float32 is predicted by a + u_32 b_k to within a factor 2 only where the linear regime holds (eps >= 1e-2).
Written 2026-09-06 by Claire.  Runs in seconds.
"""
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path
import numpy as np
from mpmath import mp, mpf, mpc, matrix as mpm
from rh1_common import trine_loop, TRINE, compose_discard

EPS_GRID = [2e-2, 1e-2, 5e-3, 2.5e-3, 1.25e-3, 6.25e-4, 3.125e-4]
U32, U64 = 2.0**-24, 2.0**-53


# ---------------- numpy backend (dtype = complex64 / complex128), identical formula sequence ----------------
def loop_np(alpha, eps, dt):
    rt = np.float32 if dt == np.complex64 else np.float64
    ANGLE = rt(np.pi/2); TAU = ANGLE*np.sqrt(rt(15))
    J = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dt)
    G1 = np.array([[0, 0, 1j], [0, 0, 0], [-1j, 0, 0]], dt); G2 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], dt)
    GS = (G1, G2, G1); gap = rt(1) + rt(eps); speed = ANGLE/TAU; omega = np.hypot(gap, speed)
    u = np.eye(3, dtype=dt)
    for arc in range(3):
        khat = ((gap*J + speed*GS[arc])/omega).astype(dt); phase = omega*TAU
        seg = (np.eye(3, dtype=dt) - 1j*np.sin(phase)*khat - 2*np.sin(phase/2)**2*(khat@khat)).astype(dt)
        u = (seg@u).astype(dt)
    s, c = np.sin(ANGLE), np.cos(ANGLE)                                   # frame(2, pi/2), same as rh1_common
    fr = np.array([[0, 1, 0], [s, 0, c], [c, 0, -s]], dt)
    out = np.eye(4, dtype=dt); out[:3, :3] = fr@u
    ca, sa = rt(np.cos(alpha)), rt(np.sin(alpha)); rot = np.eye(4, dtype=dt); rot[2:, 2:] = [[ca, -sa], [sa, ca]]
    return (rot@out@rot.T).astype(dt)


# ---------------- mpmath backend, identical formula sequence ----------------
def loop_mp(alpha, eps):
    I = mpc(0, 1); ANGLE = mp.pi/2; TAU = ANGLE*mp.sqrt(15)
    J = mpm([[0, 1, 0], [1, 0, 0], [0, 0, 0]]); G1 = mpm([[0, 0, I], [0, 0, 0], [-I, 0, 0]]); G2 = mpm([[0, 0, 0], [0, 0, -I], [0, I, 0]])
    GS = (G1, G2, G1); gap = mpf(1) + mpf(eps); speed = ANGLE/TAU; omega = mp.sqrt(gap*gap + speed*speed)
    u = mp.eye(3)
    for arc in range(3):
        khat = (gap*J + speed*GS[arc])/omega; phase = omega*TAU
        seg = mp.eye(3) - I*mp.sin(phase)*khat - 2*mp.sin(phase/2)**2*(khat*khat)
        u = seg*u
    s, c = mp.sin(ANGLE), mp.cos(ANGLE)
    fr = mpm([[0, 1, 0], [s, 0, c], [c, 0, -s]])
    m = fr*u; out = mp.eye(4)
    for i in range(3):
        for j in range(3): out[i, j] = m[i, j]
    ca, sa = mp.cos(alpha), mp.sin(alpha); rot = mp.eye(4); rot[2, 2] = ca; rot[2, 3] = -sa; rot[3, 2] = sa; rot[3, 3] = ca
    return rot*out*rot.T


def sub(M, rows, cols):
    out = mpm(len(rows), len(cols))
    for a, i in enumerate(rows):
        for b, j in enumerate(cols): out[a, b] = M[i, j]
    return out


def stack_mp(alpha_list, eps):
    K = mp.eye(2); leaks = []
    for al in alpha_list:
        L = loop_mp(al, eps); S = sub(L, (2, 3), (2, 3)); E = sub(L, (0, 1), (2, 3))
        leaks.append(E*K); K = S*K
    F = mpm(6, 2)
    for k, Fk in enumerate(leaks):
        for i in range(2):
            for j in range(2): F[2*k+i, j] = Fk[i, j]
    return F, K


def stack_np(alpha_list, eps, dt):
    loops = [loop_np(al, eps, dt) for al in alpha_list]
    K = np.eye(2, dtype=dt); leaks = []
    for L in loops:
        leaks.append(L[:2, 2:]@K); K = (L[2:, 2:]@K).astype(dt)
    return np.vstack(leaks).astype(dt), K, loops


# ---------------- the four arithmetic paths to DB_stack ----------------
def stable2_np(H):
    a, d, b = H[0, 0].real, H[1, 1].real, H[0, 1]
    return float(np.sqrt(((a-d)/2)**2 + abs(b)**2)/((a+d)/2))


def paths_np(F, K):
    dt = F.dtype; out = {}
    s = np.linalg.svd(F, compute_uv=False); p = s*s
    out['k1'] = float((p[0]-p[1])/(p[0]+p[1]))
    H = (F.conj().T@F).astype(dt); out['k2'] = stable2_np(H)
    a, d, b = H[0, 0].real, H[1, 1].real, H[0, 1]; tr = a+d; det = a*d - abs(b)**2; rad = tr*tr - 4*det
    out['k3'] = float(np.sqrt(rad)/tr) if rad >= 0 else float('nan')
    Hp = (np.eye(2, dtype=dt) - (K.conj().T@K)).astype(dt); out['k4'] = stable2_np(Hp)
    out['c'] = float((a+d).real/2)
    return out


def stable2_mp(H):
    a, d, b = H[0, 0].real, H[1, 1].real, H[0, 1]
    return mp.sqrt(((a-d)/2)**2 + abs(b)**2)/((a+d)/2)


def paths_mp(F, K):
    out = {}
    try:
        s = mp.svd_c(F, compute_uv=False); p = [x*x for x in s]
        out['k1'] = (p[0]-p[1])/(p[0]+p[1])
    except Exception as ex:                                   # svd path unavailable at this precision: record, do not fake
        out['k1'] = None; out['k1_error'] = repr(ex)
    H = F.H*F; out['k2'] = stable2_mp(H)
    a, d, b = H[0, 0].real, H[1, 1].real, H[0, 1]; tr = a+d; det = a*d - abs(b)**2; rad = tr*tr - 4*det
    out['k3'] = mp.sqrt(rad)/tr if rad >= 0 else None
    out['k4'] = stable2_mp(mp.eye(2) - K.H*K)
    out['c'] = (a+d)/2
    return out                                                # mpf values; stringified only at receipt time


def check(ch, name, value, tol):
    ok = bool(np.isfinite(value) and abs(value) <= tol)
    ch.append(dict(name=name, value=float(value), tolerance=tol, passed=ok))
    print(f"  [{'ok' if ok else 'FAIL'}] {name}: {value:.3e} (tol {tol:.0e})", flush=True)


def run():
    t0 = time.time(); ch = []
    res = dict(eps=EPS_GRID, cella=dict(applied=['U-0784 thm:precision-separation', 'U-0786 thm:robust-transfer',
               'U-0799 prop:bias', 'U-0776 lem:nonident'], source='cella:Papers_Library/01_completed_papers/'
               'precision_flow_and_transfer_functions/transfer_function_exponent_family_v6.tex'), checks=ch)
    print("== consistency: float64 closed form vs rh1_common.trine_loop ==", flush=True)
    for eps in (1e-2, 1e-3):
        for al in TRINE:
            check(ch, f"loop_np(complex128) == trine_loop  alpha={float(al):.4f} eps={eps}",
                  float(np.abs(loop_np(float(al), eps, np.complex128) - trine_loop(float(al), eps)).max()), 1e-13)
    PATHS = ('k1', 'k2', 'k3', 'k4'); table = {}
    print("\n== DB_stack by precision and path ==", flush=True)
    print(f"{'eps':>9} {'prec':>5} {'u':>9} | {'k1 svd':>13} {'k2 stable':>13} {'k3 tr/det':>13} {'k4 I-KK':>13}", flush=True)
    for eps in EPS_GRID:
        row = {}
        for name, dt, u in (('f32', np.complex64, U32), ('f64', np.complex128, U64)):
            F, K, _ = stack_np(TRINE, eps, dt); r = paths_np(F, K); r['u'] = u; row[name] = r
        for dps in (30, 50):
            mp.dps = dps
            F, K = stack_mp([mp.pi/3, 2*mp.pi/3, mp.pi/3], eps); r = paths_mp(F, K); r['u'] = float(mp.eps); row[f'mp{dps}'] = r
        table[eps] = row
        for name in ('f32', 'f64', 'mp30', 'mp50'):
            r = row[name]; vals = [r[k] for k in PATHS]
            fmt = lambda v: '   (none)    ' if v is None else f"{float(v):13.6e}"
            print(f"{eps:9.4e} {name:>5} {r['u']:9.1e} | " + ' '.join(fmt(v) for v in vals), flush=True)
    # ---- Theorem 2 (U-0784): C_{p,k} = a + u_p b_k.  a from the two high precisions; b_k from float64; predict float32.
    print("\n== Theorem 2 fit: C = DB/eps^3 = a + u b_k ==", flush=True)
    print(f"{'eps':>9} {'a (mp50)':>12} {'|a30-a50|/a':>11} | " + ' '.join(f"{'b_'+k:>11} {'f32 pred/meas':>14}" for k in PATHS), flush=True)
    fit = {}
    for eps in EPS_GRID:
        mp.dps = 50; e3 = mpf(eps)**3; row = table[eps]; f = {}
        a50 = row['mp50']['k2']/e3; a30 = row['mp30']['k2']/e3; ident = float(abs(a30-a50)/abs(a50))
        check(ch, f"eps={eps}: mp30 and mp50 agree on a (identifiability)", ident, 1e-20)
        f['a'] = float(a50); f['a_str'] = mp.nstr(a50, 30); f['ident_rel'] = ident; f['paths'] = {}
        if row['mp50']['k1'] is not None:
            check(ch, f"eps={eps}: k1 == k2 at 50 dps (paths agree when conditioning is gone)", float(abs(row['mp50']['k1']-row['mp50']['k2'])/row['mp50']['k2']), 1e-20)
        for k in PATHS:
            c64 = row['f64'][k]/eps**3; c32 = row['f32'][k]/eps**3
            b = (c64 - f['a'])/U64 if np.isfinite(c64) else float('nan')
            pred32 = f['a'] + U32*b if np.isfinite(b) else float('nan')
            f['paths'][k] = dict(C64=c64, C32=c32, b=b, rel_err_64=float(abs(c64-f['a'])/f['a']) if np.isfinite(c64) else None,
                                 f32_pred=pred32, f32_pred_over_meas=(pred32/c32 if np.isfinite(pred32) and np.isfinite(c32) and c32 != 0 else None))
        fit[eps] = f
        print(f"{eps:9.4e} {f['a']:12.6e} {ident:11.1e} | " + ' '.join(
            f"{f['paths'][k]['b']:11.2e} {str(round(f['paths'][k]['f32_pred_over_meas'], 3)) if f['paths'][k]['f32_pred_over_meas'] is not None else 'n/a':>14}" for k in PATHS), flush=True)
    # ---- U-0786 / U-0799: local log-log slope of DB(eps) at 50 dps and its drift toward the integer exponent
    print("\n== log-log slope of DB_stack (k2, 50 dps) between adjacent eps; and of the floor c DB^2/6 ==", flush=True)
    mp.dps = 50; slopes = []
    for e1, e2 in zip(EPS_GRID[:-1], EPS_GRID[1:]):
        d1, d2 = table[e1]['mp50']['k2'], table[e2]['mp50']['k2']; c1, c2 = table[e1]['mp50']['c'], table[e2]['mp50']['c']
        al = float(mp.log(d2/d1)/mp.log(mpf(e2)/mpf(e1))); fl = float(mp.log((c2*d2*d2)/(c1*d1*d1))/mp.log(mpf(e2)/mpf(e1)))
        mid = float(mp.sqrt(mpf(e1)*mpf(e2))); slopes.append(dict(eps_mid=mid, slope_DB=al, slope_floor=fl, slope_c=float(mp.log(c2/c1)/mp.log(mpf(e2)/mpf(e1)))))
        print(f"  eps_mid {mid:9.3e}: slope DB = {al:.9f}   slope c = {slopes[-1]['slope_c']:.9f}   slope floor = {fl:.9f}", flush=True)
    dev = np.array([s['slope_DB']-3 for s in slopes]); mids = np.array([s['eps_mid'] for s in slopes])
    rho = None
    if np.all(dev > 0) or np.all(dev < 0):
        p = np.polyfit(np.log(mids), np.log(np.abs(dev)), 1); rho = float(p[0])
        print(f"  two-term bias (U-0799): slope - 3 = kappa eps^rho with rho = {rho:.4f}, sign {'+' if dev[0] > 0 else '-'}", flush=True)
    else:
        print("  slope - 3 changes sign across the grid: no single-power correction fitted (recorded)", flush=True)
    res.update(table={str(e): {p: {k: (mp.nstr(v, 40) if isinstance(v, (mpf, mpc)) else v) for k, v in r.items()} for p, r in row.items()}
                      for e, row in table.items()},
               theorem2={str(e): f for e, f in fit.items()}, slopes=slopes, rho=rho)
    # ---- outcomes for P1..P4 (reported, not checked)
    a_vals = np.array([fit[e]['a'] for e in EPS_GRID])
    res['outcomes'] = dict(
        P1_a_range=[float(a_vals.min()), float(a_vals.max())], P1_a_rel_spread=float(np.ptp(a_vals)/a_vals.mean()),
        P2_b_by_path_at_1e_2={k: fit[1e-2]['paths'][k]['b'] for k in PATHS},
        P2_k4_rel_err_64={str(e): fit[e]['paths']['k4']['rel_err_64'] for e in EPS_GRID},
        P3_slope_limit=slopes[-1]['slope_DB'], P3_rho=rho,
        P4_f32_pred_over_meas={str(e): {k: fit[e]['paths'][k]['f32_pred_over_meas'] for k in PATHS} for e in EPS_GRID})
    files = ['rh2_precision.py', 'rh1_common.py', 'reflection_loop_dynamics.py', 'reflection_loop_composite.py', 'reflection_loop_reference_echo.py']
    res['provenance'] = dict(test='rh2_precision', date='2026-09-06', runtime_s=round(time.time()-t0, 1),
        head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        blinding='designer-run; predictions P1-P4 in the module docstring, committed before the log',
        sha256={p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in files})
    res['check_count'] = len(ch); res['checks_passed'] = int(sum(c['passed'] for c in ch))
    ap = argparse.ArgumentParser(); ap.add_argument('--json', nargs='?', const='docs/rh2_precision-checks.json', default=None); a = ap.parse_args()
    if a.json: Path(a.json).write_text(json.dumps(res, indent=2, allow_nan=True)+'\n')
    print(f"\nRESULT: {res['checks_passed']}/{res['check_count']} checks passed.  runtime {res['provenance']['runtime_s']}s", flush=True)
    if res['checks_passed'] != res['check_count']: sys.exit(1)


if __name__ == '__main__':
    run()
