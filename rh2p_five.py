"""RH-2P/5 -- precision-scaling separation (Cella U-0784) on the five-loop word, plus the exact closed form.
Five-loop gain-only word (rh1_common.word_loops('five') at pure gain): loops = rot(a) f^(+-) rot(a)^T with f = the
n = 1 primitive (= trine_loop at alpha = 0, checked to 6e-17), angles pi*(1/6, 5/6, 7/6, 3/4, 1/4), signs (+,+,+,-,-),
inverse = conjugate transpose (stretches rescale time and leave the propagator unchanged).
Backends and paths reused from rh2_precision.py.  Precisions complex64, complex128, 30 dps, 50 dps.
CLOSED FORM (RH-2A machinery, same sigma series, five angles, sigma -> conj(sigma) on inverse loops):
    c0 = 5625 pi^2/2048 = (25/2) x^2,   a0 = 2025 sqrt(5) pi^3/32768 = (3 sqrt5/5) x^3 = (3/5) a0_trine,
    kappa1 = -55/32 (same as trine),   floor = (15/4) x^8 = 82.9393 eps^8,   x = (15 pi/32) eps.
Written 2026-09-06 by Claire after the closed form was derived (see docs/2026-09-06-RH-2A-a0-closed-form.md addendum).
"""
import json, hashlib, subprocess, sys, time
from pathlib import Path
import numpy as np
from mpmath import mp, mpf, matrix as mpm
import rh2_precision as RP
from rh1_common import word_loops

ANG = [1/6, 5/6, 7/6, 3/4, 1/4]; ANG_FRAC = [(1, 6), (5, 6), (7, 6), (3, 4), (1, 4)]; SGN = [1, 1, 1, -1, -1]
EPS_GRID = RP.EPS_GRID; U32, U64 = RP.U32, RP.U64
A0_SYM = lambda: 2025*mp.sqrt(5)*mp.pi**3/32768
C0_SYM = lambda: 5625*mp.pi**2/2048


def five_np(eps, dt):
    rt = np.float32 if dt == np.complex64 else np.float64
    f = RP.loop_np(0.0, eps, dt); loops = []
    for a, s in zip(ANG, SGN):
        ca, sa = rt(np.cos(np.pi*a)), rt(np.sin(np.pi*a))
        rot = np.eye(4, dtype=dt); rot[2:, 2:] = [[ca, -sa], [sa, ca]]
        loops.append((rot@(f if s > 0 else f.conj().T)@rot.T).astype(dt))
    K = np.eye(2, dtype=dt); leaks = []
    for L in loops:
        leaks.append(L[:2, 2:]@K); K = (L[2:, 2:]@K).astype(dt)
    return np.vstack(leaks).astype(dt), K, loops


def five_mp(eps):
    f = RP.loop_mp(0, eps); K = mp.eye(2); leaks = []
    for (nu, de), s in zip(ANG_FRAC, SGN):
        th = mp.pi*mpf(nu)/de; ca, sa = mp.cos(th), mp.sin(th)
        rot = mp.eye(4); rot[2, 2] = ca; rot[2, 3] = -sa; rot[3, 2] = sa; rot[3, 3] = ca
        L = rot*(f if s > 0 else f.H)*rot.T
        S = RP.sub(L, (2, 3), (2, 3)); E = RP.sub(L, (0, 1), (2, 3))
        leaks.append(E*K); K = S*K
    F = mpm(10, 2)
    for k, Fk in enumerate(leaks):
        for i in range(2):
            for j in range(2): F[2*k+i, j] = Fk[i, j]
    return F, K


def run():
    t0 = time.time(); ch = []
    def check(name, value, tol):
        ok = bool(np.isfinite(value) and abs(value) <= tol); ch.append(dict(name=name, value=float(value), tolerance=tol, passed=ok))
        print(f"  [{'ok' if ok else 'FAIL'}] {name}: {value:.3e} (tol {tol:.0e})", flush=True)
    print("== consistency: float64 five-loop stack vs rh1_common.word_loops('five') ==", flush=True)
    for eps in (1e-2, 1e-3):
        ref = word_loops('five', [(eps, 0., 0.)])[0]; _, _, loops = five_np(eps, np.complex128)
        check(f"five_np(complex128) == word_loops five  eps={eps}", float(max(np.abs(a - b).max() for a, b in zip(loops, ref))), 1e-13)
    PATHS = ('k1', 'k2', 'k3', 'k4'); table = {}
    print("\n== DB_stack (five-loop) by precision and path ==", flush=True)
    print(f"{'eps':>9} {'prec':>5} {'u':>9} | {'k1 svd':>13} {'k2 stable':>13} {'k3 tr/det':>13} {'k4 I-KK':>13}", flush=True)
    for eps in EPS_GRID:
        row = {}
        for name, dt, u in (('f32', np.complex64, U32), ('f64', np.complex128, U64)):
            F, K, _ = five_np(eps, dt); r = RP.paths_np(F, K); r['u'] = u; row[name] = r
        for dps in (30, 50):
            mp.dps = dps; F, K = five_mp(eps); r = RP.paths_mp(F, K); r['u'] = float(mp.eps); row[f'mp{dps}'] = r
        table[eps] = row
        for name in ('f32', 'f64', 'mp30', 'mp50'):
            r = row[name]; fmt = lambda v: '   (none)    ' if v is None else f"{float(v):13.6e}"
            print(f"{eps:9.4e} {name:>5} {r['u']:9.1e} | " + ' '.join(fmt(r[k]) for k in PATHS), flush=True)
    print("\n== Theorem 2 (U-0784): a identified from mp30/mp50, b_k from float64 ==", flush=True)
    fit = {}
    for eps in EPS_GRID:
        mp.dps = 50; e3 = mpf(eps)**3; row = table[eps]
        a50 = row['mp50']['k2']/e3; a30 = row['mp30']['k2']/e3; ident = float(abs(a30 - a50)/abs(a50))
        check(f"eps={eps}: mp30 and mp50 agree on a", ident, 1e-20)
        if row['mp50']['k1'] is not None:
            check(f"eps={eps}: k1 == k2 at 50 dps", float(abs(row['mp50']['k1'] - row['mp50']['k2'])/row['mp50']['k2']), 1e-20)
        f = dict(a=float(a50), ident_rel=ident, paths={})
        for k in PATHS:
            c64 = row['f64'][k]/eps**3; b = (c64 - f['a'])/U64 if np.isfinite(c64) else float('nan')
            f['paths'][k] = dict(C64=c64, b=b, rel_err_64=float(abs(c64 - f['a'])/f['a']) if np.isfinite(c64) else None)
        fit[eps] = f
        print(f"  eps={eps:9.4e}  a={f['a']:.9f}  " + '  '.join(f"b_{k}={f['paths'][k]['b']:9.2e}" for k in PATHS), flush=True)
    # slopes at 50 dps and Richardson limit against the closed form
    mp.dps = 50; slopes = []
    for e1, e2 in zip(EPS_GRID[:-1], EPS_GRID[1:]):
        d1, d2 = table[e1]['mp50']['k2'], table[e2]['mp50']['k2']; c1, c2 = table[e1]['mp50']['c'], table[e2]['mp50']['c']
        slopes.append(dict(eps_mid=float(mp.sqrt(mpf(e1)*mpf(e2))), slope_DB=float(mp.log(d2/d1)/mp.log(mpf(e2)/mpf(e1))),
                           slope_floor=float(mp.log((c2*d2*d2)/(c1*d1*d1))/mp.log(mpf(e2)/mpf(e1)))))
        print(f"  slope DB = {slopes[-1]['slope_DB']:.9f}   slope floor = {slopes[-1]['slope_floor']:.9f}   (eps_mid {slopes[-1]['eps_mid']:.3e})", flush=True)
    mp.dps = 70; eg = [mpf(2)**(-k) for k in range(8, 17)]; a = []; cc = []
    for e in eg:
        F, K = five_mp(e); r = RP.paths_mp(F, K); a.append(r['k2']/e**3); cc.append(r['c']/e**2)
    def rich(v):
        R1 = [2*v[i+1] - v[i] for i in range(len(v) - 1)]
        R2 = [(4*R1[i+1] - R1[i])/3 for i in range(len(R1) - 1)]
        return [(8*R2[i+1] - R2[i])/7 for i in range(len(R2) - 1)]
    a0_num, c0_num = rich(a)[-1], rich(cc)[-1]
    a0s, c0s = A0_SYM(), C0_SYM()
    k1_num = (a[-1]/a0s - 1)/eg[-1]
    print(f"\n  a0 Richardson = {mp.nstr(a0_num, 18)}   closed form 2025 sqrt5 pi^3/32768 = {mp.nstr(a0s, 18)}   rel {mp.nstr(abs(a0_num-a0s)/a0s, 4)}")
    print(f"  c0 Richardson = {mp.nstr(c0_num, 18)}   closed form 5625 pi^2/2048      = {mp.nstr(c0s, 18)}   rel {mp.nstr(abs(c0_num-c0s)/c0s, 4)}")
    print(f"  kappa1 numerical = {mp.nstr(k1_num, 8)}   closed form -55/32 = {-55/32}")
    check("a0 Richardson == closed form (rel 1e-12)", float(abs(a0_num - a0s)/a0s), 1e-12)
    check("c0 Richardson == closed form (rel 1e-12)", float(abs(c0_num - c0s)/c0s), 1e-12)
    check("kappa1 numerical within 0.1% of -55/32", float(abs(k1_num/mpf(-55)*32 - 1)), 1e-3)
    check("a0_five / a0_trine == 3/5 (RH-2P value 7.1409906174871697658)", float(a0s/mpf('7.1409906174871697658') - mpf(3)/5), 1e-15)
    # RH-2 regression: measured DB_stack at 1e-2 was 4.204e-6
    rec = json.load(open('docs/rh2-checks.json'))
    db_rh2 = rec['words']['five']['rows'][1]['DB_stack']
    check("RH-2 five DB_stack(1e-2) == float64 k2 here", float(abs(db_rh2 - table[1e-2]['f64']['k2'])/db_rh2), 1e-9)
    files = ['rh2p_five.py', 'rh2_precision.py', 'rh1_common.py', 'reflection_loop_dynamics.py', 'reflection_loop_composite.py']
    res = dict(word='five', angles_over_pi=ANG, signs=SGN, eps=EPS_GRID, closed_form=dict(a0='2025*sqrt(5)*pi^3/32768', c0='5625*pi^2/2048', kappa1='-55/32', floor='(15/4)(15 pi/32)^8 = 38443359375 pi^8/4398046511104'),
               table={str(e): {p: {k: (mp.nstr(v, 40) if not isinstance(v, (float, type(None), str)) else v) for k, v in r.items()} for p, r in row.items()} for e, row in table.items()},
               theorem2={str(e): f for e, f in fit.items()}, slopes=slopes,
               richardson=dict(a0=mp.nstr(a0_num, 25), c0=mp.nstr(c0_num, 25), kappa1=mp.nstr(k1_num, 12)),
               checks=ch, check_count=len(ch), checks_passed=int(sum(c['passed'] for c in ch)),
               provenance=dict(test='rh2p_five', date='2026-09-06', runtime_s=round(time.time()-t0, 1),
                               head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                               sha256={p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in files}))
    Path('docs/rh2p_five-checks.json').write_text(json.dumps(res, indent=2, allow_nan=True) + '\n')
    print(f"\nRESULT: {res['checks_passed']}/{res['check_count']} checks passed.  runtime {res['provenance']['runtime_s']}s", flush=True)
    if res['checks_passed'] != res['check_count']: sys.exit(1)


if __name__ == '__main__':
    run()
