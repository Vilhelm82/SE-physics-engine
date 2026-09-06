"""RH-2A -- closed form for a0, the cubic coefficient of the tight-frame deviation DB_stack = a0 eps^3 + ...
Derived 2026-09-06 by Claire.  Runs in about a minute (SymPy series of the primitives dominates).

REDUCTION.  At exact return every arc is a 2pi spin-1 rotation, so at eps = 0 each loop is the frame flip: a
reflection I - 2 u u^T on the code plane P.  At eps != 0 the loop is unitary on the 3-frame (p, q, D) and fixes the
spectator, so on P it is EXACTLY  S = I - (1 - sigma) u u^T  and into Q it is  E = ell u^T, with one complex number
sigma(eps) = <D|L|D>, one bright vector ell(eps), |sigma|^2 + |ell|^2 = 1.  Hence
    F^dag F = |ell|^2 sum_k conj(w_k) w_k^T,   w_1 = u_1, w_2 = S_1 u_2, w_3 = S_1 S_2 u_3,
and DB_stack is a function of sigma and the three angles only.
ROUTE A (exact): sigma(eps) as a series from the closed-form loop (rh1_common.trine_loop's algebra: segment/frame),
then the 2x2 frame algebra at (pi/3, 2pi/3, pi/3), truncated-series arithmetic, leading order of t^2 + |w|^2.
ROUTE B (numerical): third-order Richardson on DB/eps^3 from rh2_precision's 70-digit closed form.
RESULT.  With  x := (15 pi / 32) eps  (= q0 * dphi, the first-order arc over-rotation projected by the mixing ratio
q0 = 1/4):   c = (15/2) x^2 + O(eps^3),   DB = sqrt(5) x^3 (1 - (55/32) eps + ...),   floor c DB^2/6 = (25/4) x^8.
So  a0 = 3375 sqrt(5) pi^3 / 32768 = 7.14099061748717...,  c0 = 3375 pi^2 / 2048,  and  (15 pi / 32) is FD-5's b_1.
Why cubic: sigma = -1 + (1125 pi^2/2048) eps^2 + (1125 pi^2/32768)(-29 + 6 sqrt(15) i pi) eps^3 + ...  The eps^2 term
is REAL -- a real shift of sigma is a partial reflection and leaves the tight frame tight; the first IMAGINARY part
(the phase of sigma) is at eps^3, and only the phase breaks the frame.  sqrt(5) is |(2,-1)|: the bright direction
from arcs (G1, G2, G1) couples D to p twice and to q once.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))  # repo root on sys.path (reorg 2026-09-06)
import sympy as sp, time, json, hashlib, subprocess, sys
from pathlib import Path
from mpmath import mp, mpf

t0 = time.time()
eps = sp.symbols('epsilon', positive=True)
S1, C1, p, q = sp.symbols('S1 C1 p q')
J = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
G1 = sp.Matrix([[0, 0, sp.I], [0, 0, 0], [-sp.I, 0, 0]]); G2 = sp.Matrix([[0, 0, 0], [0, 0, -sp.I], [0, sp.I, 0]])
def seg(G):
    k = p*J + q*G
    return sp.eye(3) - sp.I*S1*k + (C1 - 1)*(k*k)
U = seg(G1)*seg(G2)*seg(G1)                                   # arcs 0,1,2 with GS = (G1, G2, G1)
R = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, -1]])             # frame(2, pi/2)
M = (R*U).applyfunc(sp.expand)
sigma_sym = sp.expand(M[2, 2]); ell_sym = [sp.expand(M[0, 2]), sp.expand(M[1, 2])]

NP = 7                                                         # primitive series order
g = 1 + eps; s = 1/sp.sqrt(15); tau = sp.pi/2*sp.sqrt(15)
omega = sp.sqrt(g**2 + s**2); phi = omega*tau
delta = sp.series(sp.expand(phi/(2*sp.pi)) - 1, eps, 0, NP).removeO()
tr = lambda e: sp.expand(sp.series(sp.expand(e), eps, 0, NP).removeO())
z = 2*sp.pi*delta
S1s = tr(sum((-1)**k*z**(2*k + 1)/sp.factorial(2*k + 1) for k in range(4)))
C1s = tr(sum((-1)**k*z**(2*k)/sp.factorial(2*k) for k in range(4)))
ps = tr(sp.series(g/omega, eps, 0, NP).removeO()); qs = tr(sp.series(s/omega, eps, 0, NP).removeO())
print(f"[{time.time()-t0:5.1f}s] primitives: delta = {delta}", flush=True)

# ---- truncated series arithmetic: dict {power: coeff}
N = 14                                                         # D-stage truncation; sigma exact to eps^6 -> D2 trustworthy to eps^13
def todict(expr):
    P = sp.Poly(sp.expand(expr), eps); return {m[0]: sp.expand(c) for m, c in P.terms() if m[0] < N}
def add(a, b):
    out = dict(a)
    for k, v in b.items(): out[k] = sp.expand(out.get(k, 0) + v)
    return out
def mul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            if i + j < N: out[i + j] = sp.expand(out.get(i + j, 0) + x*y)
    return out
def scal(c, a): return {k: sp.expand(c*v) for k, v in a.items()}
def conj(a): return {k: sp.expand(sp.conjugate(v)) for k, v in a.items()}
def const(c): return {0: sp.sympify(c)}
def powr(a, n):
    out = const(1)
    for _ in range(n): out = mul(out, a)
    return out
prim = {S1: todict(S1s), C1: todict(C1s), p: todict(ps), q: todict(qs)}
def evalpoly(expr):
    P = sp.Poly(expr, S1, C1, p, q); acc = {}
    for (a1, a2, a3, a4), c in P.terms():
        term = const(c)
        for sym, e in ((S1, a1), (C1, a2), (p, a3), (q, a4)):
            if e: term = mul(term, powr(prim[sym], e))
        acc = add(acc, term)
    return acc
sigma = evalpoly(sigma_sym)
sigma = {k: v for k, v in sigma.items() if k <= NP - 1}         # only the first NP orders are complete
print(f"[{time.time()-t0:5.1f}s] sigma:", flush=True)
for k in range(NP): print(f"    sigma[{k}] = {sp.simplify(sigma.get(k, 0))}", flush=True)
ell2 = add(const(1), scal(-1, mul(sigma, conj(sigma))))        # |ell|^2 = 1 - |sigma|^2 (unitarity of the 3-frame column)
ell2 = {k: v for k, v in ell2.items() if k <= NP - 1}

# ---- 2x2 frame algebra at the trine angles
angles = [sp.pi/3, 2*sp.pi/3, sp.pi/3]
u = [(sp.cos(a), sp.sin(a)) for a in angles]
oms = add(const(1), scal(-1, sigma))                           # 1 - sigma
def S_apply(uk, vec):
    dot = add(scal(uk[0], vec[0]), scal(uk[1], vec[1])); f = mul(oms, dot)
    return [add(vec[0], scal(-uk[0], f)), add(vec[1], scal(-uk[1], f))]
w1 = [const(u[0][0]), const(u[0][1])]
w2 = S_apply(u[0], [const(u[1][0]), const(u[1][1])])
w3 = S_apply(u[0], S_apply(u[1], [const(u[2][0]), const(u[2][1])]))
H = [[{}, {}], [{}, {}]]
for w in (w1, w2, w3):
    for i in range(2):
        for j in range(2): H[i][j] = add(H[i][j], mul(conj(w[i]), w[j]))
H = [[mul(ell2, H[i][j]) for j in range(2)] for i in range(2)]
c = scal(sp.Rational(1, 2), add(H[0][0], H[1][1]))
t = scal(sp.Rational(1, 2), add(H[0][0], scal(-1, H[1][1])))
wv = H[0][1]
D2 = add(mul(t, t), mul(wv, conj(wv)))
print(f"[{time.time()-t0:5.1f}s] frame algebra done", flush=True)

# ---- results and checks
c0 = sp.simplify(c[2]); c3 = sp.simplify(c[3])
lead = min(k for k in D2 if sp.simplify(D2[k]) != 0)
coef = sp.simplify(D2[lead]); nxt = sp.simplify(D2.get(lead + 1, 0))
a0 = sp.simplify(sp.radsimp(sp.sqrt(coef)/c0))
kappa1 = sp.simplify(nxt/(2*coef) - c3/c0)
x_over_eps = sp.Rational(15, 32)*sp.pi                         # x = (15 pi/32) eps = q0 * dphi
claims = {
    'sigma2_real': sp.simplify(sp.im(sigma[2])) == 0,
    'sigma3_imag_nonzero': sp.simplify(sp.im(sigma[3])) != 0,
    'c0 = 3375 pi^2/2048': sp.simplify(c0 - sp.Rational(3375, 2048)*sp.pi**2) == 0,
    'c0 = (15/2) x^2': sp.simplify(c0 - sp.Rational(15, 2)*x_over_eps**2) == 0,
    'D2 leading order = 10': lead == 10,
    'a0 = 3375 sqrt5 pi^3/32768': sp.simplify(a0 - sp.Rational(3375, 32768)*sp.sqrt(5)*sp.pi**3) == 0,
    'a0 = sqrt5 x^3': sp.simplify(a0 - sp.sqrt(5)*x_over_eps**3) == 0,
    'floor c0 a0^2/6 = (25/4) x^8': sp.simplify(c0*a0**2/6 - sp.Rational(25, 4)*x_over_eps**8) == 0,
    'kappa1 = -55/32': sp.simplify(kappa1 + sp.Rational(55, 32)) == 0,
}
print(f"\n[{time.time()-t0:5.1f}s] RESULTS", flush=True)
print(f"  c0     = {c0}  = {sp.N(c0, 18)}")
print(f"  a0     = {a0}  = {sp.N(a0, 18)}")
print(f"  kappa1 = {kappa1}   (DB = a0 eps^3 (1 + kappa1 eps + ...))")
print(f"  floor  = c0 a0^2/6 = {sp.simplify(c0*a0**2/6)} = {sp.N(c0*a0**2/6, 12)} eps^8")
for k, v in claims.items(): print(f"  [{'ok' if v else 'FAIL'}] {k}")

# ---- route B: 70-digit Richardson on the closed-form propagator (rh2_precision), and the RH-2P receipt
import rlq.rh2_precision as RP
mp.dps = 70
eg = [mpf(2)**(-k) for k in range(8, 17)]; a = []
for e in eg:
    F, K = RP.stack_mp([mp.pi/3, 2*mp.pi/3, mp.pi/3], e); a.append(RP.paths_mp(F, K)['k2']/e**3)
R1 = [2*a[i+1]-a[i] for i in range(len(a)-1)]; R2 = [(4*R1[i+1]-R1[i])/3 for i in range(len(R1)-1)]; R3 = [(8*R2[i+1]-R2[i])/7 for i in range(len(R2)-1)]
a0_num = R3[-1]; a0_sym = mpf(sp.N(a0, 40).__str__())
rel = abs(a0_num - a0_sym)/a0_sym
k1_num = (a[-1]/a0_sym - 1)/eg[-1]
rec = json.load(open('docs/receipts/rh2_precision-checks.json'))
a0_rh2p = float(rec['theorem2']['0.0003125']['a'])           # smallest-eps a(eps) from RH-2P, not the limit
print(f"\n  route B (Richardson, 70 dps): {mp.nstr(a0_num, 18)}   rel diff to route A: {mp.nstr(rel, 4)}")
print(f"  kappa1 numerical at eps=2^-16: {mp.nstr(k1_num, 8)}   (series: {float(kappa1):.6f})")
print(f"  RH-2P a(3.1e-4) = {a0_rh2p:.7f}; series prediction a0 (1 + kappa1 eps) = {float(a0*(1 + kappa1*sp.Rational(1, 3200))):.7f}")
claims['route A == route B to 1e-12'] = bool(rel < mpf('1e-12'))
claims['kappa1 numerical within 0.1% of -55/32'] = bool(abs(k1_num/mpf(-55)*32 - 1) < mpf('1e-3'))
claims['RH-2P a(3.1e-4) reproduced to 1e-5 by a0(1+kappa1 eps)'] = abs(a0_rh2p - float(a0*(1 + kappa1*sp.Rational(1, 3200)))) < 1e-5
files = ['experiments/2026-09-06/rh2_a0_closed_form.py', 'rlq/rh2_precision.py', 'rlq/rh1_common.py', 'rlq/reflection_loop_dynamics.py']
out = dict(a0=str(a0), a0_decimal=str(sp.N(a0, 30)), c0=str(c0), kappa1=str(kappa1), floor=str(sp.simplify(c0*a0**2/6)),
           x_over_eps='15*pi/32', sigma={k: str(sp.simplify(v)) for k, v in sigma.items()},
           D2={k: str(sp.simplify(v)) for k, v in D2.items() if k >= 10}, c={k: str(sp.simplify(v)) for k, v in c.items()},
           route_B_richardson=mp.nstr(a0_num, 25), route_rel_diff=mp.nstr(rel, 5),
           cella_pin=dict(name='a0_trine_tight_frame_deviation_cubic_coefficient', digest='566bbaea565eabb4', certified_digits=14),
           checks={k: bool(v) for k, v in claims.items()}, check_count=len(claims), checks_passed=int(sum(claims.values())),
           provenance=dict(date='2026-09-06', head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                           runtime_s=round(time.time()-t0, 1), sha256={f: hashlib.sha256(Path(f).read_bytes()).hexdigest() for f in files}))
Path('docs/receipts/rh2_a0-checks.json').write_text(json.dumps(out, indent=2) + '\n')
print(f"\nRESULT: {out['checks_passed']}/{out['check_count']} checks passed.  runtime {out['provenance']['runtime_s']}s")
if out['checks_passed'] != out['check_count']: sys.exit(1)
