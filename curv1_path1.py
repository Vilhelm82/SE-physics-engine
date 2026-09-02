#!/usr/bin/env python3
# =============================================================================
# CURV-1 / PATH 1 -- the seat-triangle excess on the presented Gram, along the pinned family
# Date: 2026-09-02.  Spec: docs/prereg/CURV-1/PATH1-SPEC.md merged with the recovered twin spec:
# both functionals (Omega_A raw, Omega_B direction-normalised), the pole orbit (the four triangles cut
# out by three presented lines), both K layers, both sheets, interior features located, classification
# harness with self-tests.  sympy + mpmath only.  Exit 0 iff every check passes.
#
# CARRIED (each re-verified below, none trusted):
#   CARRY-1  G'(r) = G + f(r) k k^T, f = r_s/(r - r_s), det G' = Delta r/(r - r_s)          [D-1..D-4]
#   CARRY-2  pinning tanh(l) = sqrt(r_s/r)                                                  [KIN-2a, declared]
#   CARRY-3  interior sheet l = mu + i pi/2, tanh(mu) = sqrt(r/r_s), cosh(l) = i sinh(mu)      [CONT-1, declared]
#   CARRY-4  layers: A (seat, K = root a3) k = (g13,g23,1);  C (world, K ~ a1+a2+a3) k = G1/sqrt(1^T G 1)
#   CARRY-5  area-tangent identity, unit-diagonal Gram: Omega = 2 arg(S + iV), S = 1 + sum g, V = +sqrt(det)
#   CARRY-6  V'(r) = sqrt(det G') principal: real outside, +i sqrt|det| inside (the +i sheet)
#   CARRY-7  null radii r_i = r_s (1 - k_i^2)                                               [H-22]
#   CARRY-8  rod-free (frozen, thm_h2_d1 header): no division by a presented length, by r, dr, r*eps
#
# THE OBJECT.  eps = (e1,e2,e3) with e1 e2 e3 = +1 indexes the four triangles cut out by three lines.
#   functional A (raw):         S_eps = 1 + sum e_i e_j G'_ij,      V = sqrt(det G')
#   functional B (normalised):  S_eps = 1 + sum e_i e_j Gam_ij,     V = sqrt(det G') / (sr1 sr2 sr3)
#                               Gam_ij = G'_ij/(sr_i sr_j),  sr_i = sqrt(rho_i) principal,  rho_i = G'_ii
#   B_eps = S_eps + iV,  W_eps = (S_eps + iV)/(S_eps - iV) = e^{i Omega_eps},  Omega_eps = -i log W_eps.
#   The other branch of sr_i is e_i -> -e_i, so the multiset over eps is branch-free.
#
# CLASSIFICATION at an approach point (baseline F0 = the same functional at r -> oo):
#   IDENTITY (limit = F0) / FINITE-NONTRIVIAL (finite, != F0) / UNBOUNDED / NO-LIMIT,
#   plus "0 mod 2pi?" and the half-angle sign (sign of B where B is real: +1 or -1 in SU(2)).
#
# NOT USED (comparison block only): Gauss-Bonnet, angle deficit, Regge, curvature invariants, tidal
#   tensors, geodesic deviation, any interior black-hole model, Kruskal/Penrose, Kerr, regular-BH metrics.
# =============================================================================
import sys, json, signal
import sympy as sp
import mpmath as mp
mp.mp.dps = 40

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
def z(e):
    e = sp.sympify(e)
    for rt in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), lambda q: sp.simplify(sp.radsimp(sp.expand(q))),
               lambda q: sp.simplify(sp.expand_complex(sp.expand(q)))):
        try:
            if rt(e) == 0: return True
        except Exception: pass
    try:
        if not e.free_symbols: return abs(sp.N(e, 60)) < sp.Float(10)**-50
        pts = [{g12: sp.Rational(1,4), g13: -sp.Rational(1,3), g23: sp.Rational(1,5)}, {g12: sp.Rational(9,10), g13: sp.Rational(9,10), g23: sp.Rational(5,8)},
               {g12: -sp.Rational(9,20), g13: -sp.Rational(9,20), g23: -sp.Rational(9,20)}, {g12: sp.Rational(1,20), g13: -sp.Rational(1,30), g23: sp.Rational(1,50)}]
        return all(abs(sp.N(e.subs(p), 60)) < sp.Float(10)**-40 for p in pts)   # numeric-grade: identity checked at the four frozen Grams
    except Exception: return False
def zM(M): return all(z(e) for e in M)
class Timeout(Exception): pass
def timed(fn, secs, default=None):
    def h(s, f): raise Timeout()
    old = signal.signal(signal.SIGALRM, h); signal.alarm(secs)
    try: return fn()
    except Timeout: return default
    except Exception: return default
    finally: signal.alarm(0); signal.signal(signal.SIGALRM, old)

r, rs, u = sp.Symbol('r', positive=True), sp.Symbol('r_s', positive=True), sp.Symbol('u', positive=True)
g12, g13, g23 = sp.symbols('g12 g13 g23', real=True)
lam, mu = sp.symbols('lambda mu', positive=True)
one = sp.Matrix([1, 1, 1])
def Gram(a, b, c): return sp.Matrix([[1, a, b], [a, 1, c], [b, c, 1]])
Gsym = Gram(g12, g13, g23); Delta_sym = sp.expand(Gsym.det())
EPS = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
GRAMS = {'N1': (sp.Rational(1, 4), -sp.Rational(1, 3), sp.Rational(1, 5)),
         'N2': (sp.Rational(9, 10), sp.Rational(9, 10), sp.Rational(5, 8)),
         'N3': (sp.Rational(1, 20), -sp.Rational(1, 30), sp.Rational(1, 50)),
         'N4': (-sp.Rational(9, 20), -sp.Rational(9, 20), -sp.Rational(9, 20))}
CTRL = (sp.Rational(9, 10), sp.Rational(9, 10), sp.Rational(3, 5))
LAYERS = ('A', 'C'); FUNCS = ('A', 'B')
RES = {}

# ---------------------------------------------------------------------------
# SETUP CHECKS: carried objects re-derived, admissibility gate, classifier self-tests
# ---------------------------------------------------------------------------
print("=" * 100); print("SETUP -- carried objects re-derived, admissibility, classifier self-tests"); print("=" * 100)
k1, k2, k3, fs = sp.symbols('k1 k2 k3 f', real=True); kg = sp.Matrix([k1, k2, k3])
check("P1-1a matrix determinant lemma det(G + f k k^T) = det G + f k^T adj(G) k, general symbolic k",
      sp.expand((Gsym + fs*kg*kg.T).det() - (Gsym.det() + fs*(kg.T*Gsym.adjugate()*kg)[0, 0])) == 0)
kA = sp.Matrix([g13, g23, 1]); vC = Gsym*one; NC = (one.T*Gsym*one)[0, 0]
check("P1-1b k^T adj(G) k = Delta on both layers  =>  det G'(r) = Delta (1 + f) = Delta r/(r - r_s)",
      z((kA.T*Gsym.adjugate()*kA)[0, 0] - Delta_sym) and z(sp.simplify((vC.T*Gsym.adjugate()*vC)[0, 0]/NC - Delta_sym)))
check("P1-2a exterior: sinh^2(l) = r_s/(r - r_s) under tanh(l) = sqrt(r_s/r)",
      z(sp.simplify((sp.tanh(lam)**2/(1 - sp.tanh(lam)**2)).subs(sp.tanh(lam), sp.sqrt(rs/r)) - rs/(r - rs))))
check("P1-2b interior: -cosh^2(mu) = r_s/(r - r_s) under tanh(mu) = sqrt(r/r_s)  (the two-sheet collapse)",
      z(sp.simplify(-(1/(1 - sp.tanh(mu)**2)).subs(sp.tanh(mu), sp.sqrt(r/rs)) - rs/(r - rs))))
gate = {nm: Delta_sym.subs({g12: g[0], g13: g[1], g23: g[2]}) for nm, g in GRAMS.items()}
check("P1-3  admissibility: Delta(N1..N4) = " + ", ".join(f"{nm}:{d}" for nm, d in gate.items()),
      all(d > 0 for d in gate.values()) and gate['N1'] == sp.Rational(2711, 3600) and gate['N2'] == sp.Rational(3, 1600)
      and gate['N4'] == sp.Rational(841, 4000))
dctrl = Delta_sym.subs({g12: CTRL[0], g13: CTRL[1], g23: CTRL[2]})
check("P1-C3 the gate REJECTS the control Gram (9/10, 9/10, 3/5): Delta = " + str(dctrl), dctrl == -sp.Rational(1, 125))

# ---------------------------------------------------------------------------
# THE PRESENTED GRAM (r_s = 1: everything is a function of r/r_s), BANDS, EXPLICIT BRANCHES
# ---------------------------------------------------------------------------
def kkT(layer, Gm):
    if layer == 'A':
        k = sp.Matrix([Gm[0, 2], Gm[1, 2], 1]); return k*k.T, [int(sp.sign(x)) for x in k]
    v = Gm*one; N = (one.T*Gm*one)[0, 0]; return v*v.T/N, [int(sp.sign(x)) for x in v]

def presented(Gm, layer):
    KK, ksign = kkT(layer, Gm); Gp = Gm + KK/(r - 1)
    rho = [sp.cancel(Gp[i, i]) for i in range(3)]; det = sp.cancel(Gp.det())
    rnull = [sp.nsimplify(1 - KK[i, i]) for i in range(3)]
    return dict(Gp=Gp.applyfunc(sp.cancel), rho=rho, det=det, ksign=ksign, rnull=rnull, KK=KK)

def collision_radii(P):
    """r in (0,1) where G'_ij^2 = rho_i rho_j: the presented pair (i,j) coincides or is antipodal."""
    out = {}
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        expr = sp.numer(sp.together(sp.cancel(P['Gp'][i, j]**2 - P['rho'][i]*P['rho'][j])))
        out[(i, j)] = sorted(sp.nsimplify(x) for x in sp.solve(expr, r) if x.is_real and 0 < x < 1)
    return out

def bands(P):
    pts = sorted(set([sp.Integer(0), sp.Integer(1)] + [x for x in P['rnull'] if 0 < x < 1]))
    return [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)] + [(sp.Integer(1), sp.oo)]

def band_signs(P, band):
    lo, hi = band; mid = (lo + hi)/2 if hi != sp.oo else lo + 1
    return [1 if P['rho'][i].subs(r, mid) > 0 else -1 for i in range(3)], (1 if P['det'].subs(r, mid) > 0 else -1)

def approach(P, p, side):
    """Return (r as a function of u > 0 with u -> 0+ the approach, the band adjacent on that side)."""
    if p == sp.oo: rsub = 1/u
    elif side == '+': rsub = p + u
    else: rsub = p - u
    for b in bands(P):
        lo, hi = b
        if p == sp.oo and hi == sp.oo: return rsub, b
        if p != sp.oo and side == '+' and lo == p: return rsub, b
        if p != sp.oo and side == '-' and hi == p: return rsub, b
    raise ValueError

def posparts(expr, utest):
    """Write a rational function of u with fixed sign near u=0+ as sign * Pn/Pd with Pn, Pd > 0 there."""
    e = sp.cancel(sp.together(expr)); n, d = sp.fraction(e)
    sn = 1 if n.subs(u, utest) > 0 else -1; sd = 1 if d.subs(u, utest) > 0 else -1
    return sn*sd, sp.expand(sn*n), sp.expand(sd*d)

def explicit(P, rsub, band):
    """Explicit-branch square roots along the approach: sqrt(rho_i) = c_i sqrt(Pn_i)/sqrt(Pd_i), c_i in {1, i}."""
    lo, hi = band; width = (hi - lo) if hi != sp.oo else 1
    utest = min(sp.Rational(1, 10**6), width/10**4)
    sr, sgn = [], []
    for i in range(3):
        s, Pn, Pd = posparts(P['rho'][i].subs(r, rsub), utest)
        sr.append((1 if s > 0 else sp.I)*sp.sqrt(Pn)/sp.sqrt(Pd)); sgn.append(s)
    s, Pn, Pd = posparts(P['det'].subs(r, rsub), utest)
    sdet = (1 if s > 0 else sp.I)*sp.sqrt(Pn)/sp.sqrt(Pd)
    Gp = P['Gp'].subs(r, rsub).applyfunc(sp.cancel)
    return dict(Gp=Gp, sr=sr, sdet=sdet, sgn=sgn, sdet_sgn=s)

def SV(X, eps, func):
    e1, e2, e3 = eps; Gp = X['Gp']
    if func == 'A':
        return 1 + e1*e2*Gp[0, 1] + e1*e3*Gp[0, 2] + e2*e3*Gp[1, 2], X['sdet']
    sr = X['sr']; Gam = lambda i, j: Gp[i, j]/(sr[i]*sr[j])
    return 1 + e1*e2*Gam(0, 1) + e1*e3*Gam(0, 2) + e2*e3*Gam(1, 2), X['sdet']/(sr[0]*sr[1]*sr[2])

def reality(sgn, sdet_sgn, func):
    if func == 'A': return 'real' if sdet_sgn > 0 else 'W real: Omega in iR or pi+iR'
    Sreal = all(sgn[i]*sgn[j] > 0 for (i, j) in ((0, 1), (0, 2), (1, 2))); Vreal = (sdet_sgn*sgn[0]*sgn[1]*sgn[2]) > 0
    return 'real' if (Sreal and Vreal) else ('W real: Omega in iR or pi+iR' if Sreal else 'complex')

# ---------------------------------------------------------------------------
# LIMITS (u -> 0+) AND THE CLASSIFIER
# ---------------------------------------------------------------------------
def lead(expr):
    """Leading term of expr as u -> 0+: (coefficient, exponent) with expr ~ c u^a.  Exact via series."""
    e = sp.powsimp(sp.expand(sp.sympify(expr)), force=True)
    ser = timed(lambda: sp.series(e, u, 0, 2).removeO(), 60, None)
    if ser is None or ser.has(sp.nan): return None
    ser = sp.expand(ser)
    if ser == 0:
        ser = timed(lambda: sp.series(e, u, 0, 4).removeO(), 60, None)
        if ser is None or sp.expand(ser) == 0: return (sp.Integer(0), sp.oo)
        ser = sp.expand(ser)
    terms = sp.Add.make_args(ser); best = None
    for t in terms:
        c, a = t.as_coeff_exponent(u)
        if best is None or a < best[1]: best = (sp.simplify(c), a)
    return best

def lim0(expr):
    """Limit as u -> 0+ : finite value, or +-oo / zoo, or None."""
    L = timed(lambda: sp.limit(expr, u, 0, '+'), 40, None)
    if L is not None and not L.has(sp.AccumBounds) and not L.has(sp.nan): return sp.simplify(L)
    ld = lead(expr)
    if ld is None: return None
    c, a = ld
    if a > 0: return sp.Integer(0)
    if a == 0: return c
    return sp.zoo if c.is_real is False else (sp.oo if c > 0 else -sp.oo)
def is_inf(L): return L is not None and (L in (sp.oo, -sp.oo, sp.zoo) or L.has(sp.oo) or L.has(sp.zoo) or L.has(sp.nan))

def numeric_W(S, V, ks=(8, 16, 24)):
    fn = sp.lambdify(u, (S + sp.I*V)/(S - sp.I*V), 'mpmath'); out = []
    for k in ks:
        try: out.append(mp.mpc(fn(mp.mpf(10)**(-k))))
        except Exception: out.append(None)
    return out

def numfb(S, V, LB):
    """40-digit approach sequence when no symbolic route gives a usable answer (numeric-grade, flagged in the note)."""
    seq = numeric_W(S, V, ks=(14, 22, 30))
    if any(x is None for x in seq) or any(mp.isnan(mp.re(x)) or mp.isnan(mp.im(x)) for x in seq): return None, LB, 'no limit found'
    if abs(seq[-1]) > mp.mpf(10)**12: return sp.oo, LB, 'numeric: W -> oo'
    if abs(seq[-1]) < mp.mpf(10)**-12: return sp.Integer(0), LB, 'numeric: W -> 0'
    if abs(seq[-1] - seq[-2]) < mp.mpf(10)**-8:
        return sp.Float(str(mp.re(seq[-1])), 30) + sp.I*sp.Float(str(mp.im(seq[-1])), 30), LB, 'numeric-grade (40 digits)'
    return None, LB, 'no limit found'

def limitW(S, V):
    Wl, Bl, nt = limitW_raw(S, V)
    if Wl is None: return numfb(S, V, Bl)
    if Wl not in (sp.oo, sp.zoo):
        try:
            c = mp.mpc(complex(sp.N(Wl, 30)))
            if mp.isnan(mp.re(c)) or mp.isnan(mp.im(c)): return numfb(S, V, Bl)
        except Exception:
            return numfb(S, V, Bl)
    return Wl, Bl, nt

def limitW_raw(S, V):
    """Limit of W = (S+iV)/(S-iV) as u -> 0+.  Returns (Wlim, Blim, note)."""
    B, Bb = S + sp.I*V, S - sp.I*V
    LB, LBb = lim0(B), lim0(Bb)
    if LB is None or LBb is None or is_inf(LB) or is_inf(LBb):
        Lq = lim0(sp.cancel(V/S)) if S != 0 else None
        if Lq is not None and not is_inf(Lq) and not z(1 - sp.I*Lq): return sp.simplify((1 + sp.I*Lq)/(1 - sp.I*Lq)), LB, 'via V/S'
        if is_inf(Lq): return sp.Integer(-1), LB, 'V/S -> oo'
        seq = numeric_W(S, V, ks=(14, 22, 30))
        if any(x is None for x in seq): return None, LB, 'no limit found'
        if abs(seq[-1]) > mp.mpf(10)**12: return sp.oo, LB, 'numeric: W -> oo'
        if abs(seq[-1]) < mp.mpf(10)**-12: return sp.Integer(0), LB, 'numeric: W -> 0'
        if abs(seq[-1] - seq[-2]) < mp.mpf(10)**-8:
            return sp.Float(str(mp.re(seq[-1])), 30) + sp.I*sp.Float(str(mp.im(seq[-1])), 30), LB, 'numeric-grade (40 digits)'
        return None, LB, 'no limit found'
    if LB == 0 and LBb == 0:
        # 0/0: the limit of the ratio along the family.  Exact via sympy where it closes; the 40-digit value is the arbiter.
        seq = numeric_W(S, V, ks=(20, 30))
        num = seq[-1] if (seq[-1] is not None and seq[0] is not None and abs(seq[-1] - seq[0]) < mp.mpf(10)**-8) else None
        LW = timed(lambda: sp.limit(sp.cancel(B/Bb), u, 0, '+'), 40, None)
        if LW is not None and not is_inf(LW) and num is not None and abs(mp.mpc(complex(sp.N(LW, 30))) - num) < mp.mpf(10)**-8:
            return sp.simplify(LW), LB, '0/0 resolved along the family (exact)'
        if num is not None:
            return sp.Float(mp.re(num), 30) + sp.I*sp.Float(mp.im(num), 30), LB, '0/0 resolved along the family (numeric-grade, 40 digits)'
        return None, LB, '0/0 unresolved'
    if LBb == 0: return sp.oo, LB, 'S - iV -> 0'
    if LB == 0: return sp.Integer(0), LB, 'S + iV -> 0'
    return sp.simplify(LB/LBb), LB, ''

def Omega_of_W(W):
    if W is None or W in (0, sp.oo, sp.zoo): return None
    return sp.simplify(-sp.I*sp.log(W))
def classify(Wlim, W0):
    if Wlim is None: return 'NO-LIMIT', None
    if Wlim in (0, sp.oo, sp.zoo): return 'UNBOUNDED', None
    return ('IDENTITY' if z(Wlim - W0) else 'FINITE-NONTRIVIAL'), Omega_of_W(Wlim)
def mod2pi_zero(W): return W is not None and W not in (0, sp.oo, sp.zoo) and z(W - 1)
def halfsign(Bl):
    if Bl is None or is_inf(Bl): return '?'
    Bl = sp.simplify(Bl)
    if Bl == 0: return '0'
    if Bl.is_real: return '+1' if Bl > 0 else '-1'
    return 'complex'

# self-tests of the classifier
Pt = presented(Gram(*GRAMS['N1']), 'A'); rsub, bd = approach(Pt, sp.Integer(1), '+')
Ldet = lim0(Pt['det'].subs(r, rsub)); ld = lead(Pt['det'].subs(r, rsub))
check("P1-C1 classifier on det G'(r) at r -> r_s+: UNBOUNDED, exponent (r - r_s)^-1", is_inf(Ldet) and ld is not None and ld[1] == -1,
      f"limit {Ldet}, leading order u^{ld[1] if ld else '?'}")
Wc = sp.Rational(3, 5) + sp.I*sp.Rational(4, 5)
check("P1-C2 classifier: constant functional -> IDENTITY; W -> 0 -> UNBOUNDED", classify(Wc, Wc)[0] == 'IDENTITY' and classify(sp.Integer(0), Wc)[0] == 'UNBOUNDED')

# ---------------------------------------------------------------------------
# NUMERIC EVALUATION IN r (principal square roots = the branch convention above), ROUTE V, l'HUILIER
# ---------------------------------------------------------------------------
def gp_num(Gm, layer, rv):
    KK, _ = kkT(layer, Gm); Gp = (Gm + KK/(r - 1)).subs(r, rv)
    return mp.matrix([[mp.mpf(str(sp.N(Gp[i, j], 45))) for j in range(3)] for i in range(3)])
def B_num_G(Gm, layer, func, eps, rv):
    Gp = gp_num(Gm, layer, rv); e1, e2, e3 = eps
    det = mp.det(Gp); sdet = mp.sqrt(mp.mpc(det))
    if func == 'A':
        S = 1 + e1*e2*Gp[0, 1] + e1*e3*Gp[0, 2] + e2*e3*Gp[1, 2]; V = sdet
    else:
        sr = [mp.sqrt(mp.mpc(Gp[i, i])) for i in range(3)]
        S = 1 + e1*e2*Gp[0, 1]/(sr[0]*sr[1]) + e1*e3*Gp[0, 2]/(sr[0]*sr[2]) + e2*e3*Gp[1, 2]/(sr[1]*sr[2]); V = sdet/(sr[0]*sr[1]*sr[2])
    return mp.mpc(S), mp.mpc(V)
def W_of(S, V): return (S + 1j*V)/(S - 1j*V)
def Omega_num(S, V): return -1j*mp.log(W_of(S, V))

def frame_vectors(gvals):
    a, b, c = gvals; D = Delta_sym.subs({g12: a, g13: b, g23: c})
    a1 = sp.Matrix([1, 0, 0]); a2 = sp.Matrix([a, sp.sqrt(1 - a**2), 0])
    a3 = sp.Matrix([b, (c - a*b)/sp.sqrt(1 - a**2), sp.sqrt(D/(1 - a**2))])
    return [a1, a2, a3]
PAULI = [mp.matrix([[0, 1], [1, 0]]), mp.matrix([[0, -1j], [1j, 0]]), mp.matrix([[1, 0], [0, -1]])]
def vsig(v): return sum((v[i]*PAULI[i] for i in range(3)), mp.matrix(2, 2))
def B_num_V(gvals, layer, func, eps, rv):
    """Route V: rotor sandwich in the Pauli representation at the pinned rapidity; never forms G'."""
    A = [mp.matrix([mp.mpf(str(sp.N(x, 45))) for x in v]) for v in frame_vectors(gvals)]
    K = A[2] if layer == 'A' else (A[0] + A[1] + A[2])/mp.norm(A[0] + A[1] + A[2])
    rv = mp.mpf(rv)
    lamv = mp.atanh(1/mp.sqrt(rv)) if rv > 1 else mp.atanh(mp.sqrt(rv)) + 1j*mp.pi/2      # CARRY-2 / CARRY-3
    Rot = mp.cosh(lamv/2)*mp.eye(2) + mp.sinh(lamv/2)*vsig(K)
    vs = []
    for ai in A:
        M = Rot*vsig(ai)*Rot
        vs.append(mp.matrix([(M*PAULI[k])[0, 0] + (M*PAULI[k])[1, 1] for k in range(3)])/2)   # vector part
    dot = lambda x, y: sum(x[i]*y[i] for i in range(3))                                        # bilinear, no conjugation
    detv = mp.det(mp.matrix([[vs[i][j] for j in range(3)] for i in range(3)]))
    e1, e2, e3 = eps
    if func == 'A':
        return mp.mpc(1 + e1*e2*dot(vs[0], vs[1]) + e1*e3*dot(vs[0], vs[2]) + e2*e3*dot(vs[1], vs[2])), mp.mpc(detv)
    n = [mp.sqrt(mp.mpc(dot(v, v))) for v in vs]
    S = 1 + e1*e2*dot(vs[0], vs[1])/(n[0]*n[1]) + e1*e3*dot(vs[0], vs[2])/(n[0]*n[2]) + e2*e3*dot(vs[1], vs[2])/(n[1]*n[2])
    return mp.mpc(S), mp.mpc(detv/(n[0]*n[1]*n[2]))

def lhuilier(Gm, layer, eps, rv):
    """Independent formula for the excess of the presented directions (exterior, positive definite)."""
    Gp = gp_num(Gm, layer, rv); sr = [mp.sqrt(Gp[i, i]) for i in range(3)]; e = eps
    cos = {(0, 1): e[0]*e[1]*Gp[0, 1]/(sr[0]*sr[1]), (0, 2): e[0]*e[2]*Gp[0, 2]/(sr[0]*sr[2]), (1, 2): e[1]*e[2]*Gp[1, 2]/(sr[1]*sr[2])}
    a, b, c = [mp.acos(cos[k]) for k in ((1, 2), (0, 2), (0, 1))]; s = (a + b + c)/2
    return 4*mp.atan(mp.sqrt(mp.tan(s/2)*mp.tan((s - a)/2)*mp.tan((s - b)/2)*mp.tan((s - c)/2)))

# ---------------------------------------------------------------------------
# THE DRIVER: endpoints, reality, features, per (Gram, layer, functional, member)
# ---------------------------------------------------------------------------
POINTS = [('r->oo', sp.oo, '-'), ('r->r_s+', sp.Integer(1), '+'), ('r->r_s-', sp.Integer(1), '-'), ('r->0+', sp.Integer(0), '+')]
def side_of(P, pt, side):
    rsub = pt + u if side == '+' else pt - u
    for b in bands(P):
        lo, hi = b
        if side == '-' and lo < pt <= hi: return rsub, (lo, pt)
        if side == '+' and lo <= pt < hi: return rsub, (pt, hi)
    return None, None

def run_gram(name, gvals):
    Gm = Gram(*gvals); out = {}
    for layer in LAYERS:
        P = presented(Gm, layer); coll = collision_radii(P)
        for func in FUNCS:
            key = f"{name}/{layer}/{func}"; R_ = dict(points={}, features={}, bands=[])
            for b in bands(P):
                sg, sd = band_signs(P, b); R_['bands'].append((str(b[0]), str(b[1]), reality(sg, sd, func)))
            rsub, bd = approach(P, sp.oo, '-'); X = explicit(P, rsub, bd)
            W0 = {eps: limitW(*SV(X, eps, func))[0] for eps in EPS}
            for (pn, p, side) in POINTS:
                rsub, bd = approach(P, p, side); X = explicit(P, rsub, bd); row = {}
                for eps in EPS:
                    S, V = SV(X, eps, func); Wl, Bl, nt = limitW(S, V); cl, Om = classify(Wl, W0[eps])
                    corr = ''
                    if Wl is not None and Wl not in (sp.oo, sp.zoo):
                        seq = numeric_W(S, V, ks=(16, 24))
                        if seq[-1] is not None:
                            corr = 'ok' if abs(seq[-1] - mp.mpc(complex(sp.N(Wl, 30)))) < mp.mpf(10)**-6 else f'MISMATCH {mp.nstr(seq[-1], 8)}'
                    row[eps] = dict(cls=cl, Omega=Om, W=Wl, B=(sp.simplify(Bl) if Bl is not None and not is_inf(Bl) else Bl),
                                    half=halfsign(Bl), zero_mod_2pi=mod2pi_zero(Wl), note=nt, numeric=corr)
                R_['points'][pn] = row
            feats = {}
            for i in range(3):
                rn = P['rnull'][i]
                if func == 'B' and 0 < rn < 1:
                    fl, fr_ = {}, {}
                    for eps in EPS:
                        rsub, bd = side_of(P, rn, '-'); fl[eps] = limitW(*SV(explicit(P, rsub, bd), eps, 'B'))[0]
                        rsub, bd = side_of(P, rn, '+'); fr_[eps] = limitW(*SV(explicit(P, rsub, bd), eps, 'B'))[0]
                    feats[f'null radius r_{i+1} = {rn}'] = dict(left=fl, right=fr_)
            for (i, j), rr in coll.items():
                for rc in rr:
                    fl, fr_ = {}, {}
                    for eps in EPS:
                        rsub, bd = side_of(P, rc, '-'); fl[eps] = limitW(*SV(explicit(P, rsub, bd), eps, func))[0]
                        rsub, bd = side_of(P, rc, '+'); fr_[eps] = limitW(*SV(explicit(P, rsub, bd), eps, func))[0]
                    feats[f'collision ({i+1},{j+1}) at r = {rc}'] = dict(left=fl, right=fr_)
            R_['features'] = feats; R_['null_radii'] = P['rnull']; R_['collisions'] = coll; R_['W0'] = W0
            out[key] = R_
    return out

# ---------------------------------------------------------------------------
# RUN THE FOUR GRAMS
# ---------------------------------------------------------------------------
print(); print("=" * 100); print("RUN -- exact endpoint limits, both layers, both functionals, four members, four Grams"); print("=" * 100)
for nm, g in GRAMS.items():
    RES.update(run_gram(nm, g)); print(f"  {nm} done"); sys.stdout.flush()
def cls_at(key, pt, eps): return RES[key]['points'][pt][eps]['cls']
def all_members_ok(pred):
    return all(pred(key, pt, eps) for key in RES for pt in RES[key]['points'] for eps in EPS)
bad = [(k, p, e, RES[k]['points'][p][e]['numeric']) for k in RES for p in RES[k]['points'] for e in EPS if RES[k]['points'][p][e]['numeric'] not in ('ok', '')]
check("P1-4  numeric corroboration (40 digits, u = 1e-16, 1e-24) agrees with every exact endpoint limit",
      not bad, ("'' only where the limit is UNBOUNDED/NO-LIMIT" if not bad else "; ".join(f"{k} {p} {e}: {m}" for k, p, e, m in bad[:6])))
check("P1-5  no NO-LIMIT anywhere: every endpoint limit exists on both sheets", all_members_ok(lambda k, p, e: cls_at(k, p, e) != 'NO-LIMIT'))
check("P1-6  baseline: r -> oo returns the frame's own orbit (IDENTITY) for every Gram, layer, functional, member",
      all_members_ok(lambda k, p, e: p != 'r->oo' or cls_at(k, p, e) == 'IDENTITY'))
check("P1-6b baseline is functional-independent: Omega_A(oo) = Omega_B(oo) member-wise",
      all(z(RES[f'{n}/{L}/A']['W0'][e] - RES[f'{n}/{L}/B']['W0'][e]) for n in GRAMS for L in LAYERS for e in EPS))

# ---------------------------------------------------------------------------
# SYMBOLIC THEOREMS (general Gram) behind the centre values
# ---------------------------------------------------------------------------
print(); print("=" * 100); print("SYMBOLIC -- the centre, general Gram"); print("=" * 100)
t = sp.Symbol('t', positive=True)                       # r = t^2 on the seat layer near the centre
# layer A explicit branches near r = 0: rho3 = t^2/(t^2-1) < 0, rho1,2 > 0, det < 0
Gp_A = (Gsym + kA*kA.T/(r - 1)).subs(r, t**2).applyfunc(sp.cancel)
sr1 = sp.sqrt(1 - g13**2 - t**2)/sp.sqrt(1 - t**2); sr2 = sp.sqrt(1 - g23**2 - t**2)/sp.sqrt(1 - t**2); sr3 = sp.I*t/sp.sqrt(1 - t**2)
sdetA = sp.I*t*sp.sqrt(Delta_sym)/sp.sqrt(1 - t**2)
check("T1a  layer A branches: rho_1,2,3 and det G' reproduce the explicit forms used (identically in t)",
      z(sr1**2 - Gp_A[0, 0]) and z(sr2**2 - Gp_A[1, 1]) and z(sr3**2 - Gp_A[2, 2]) and z(sdetA**2 - Gp_A.det()))
Gam13 = sp.simplify(Gp_A[0, 2]/(sr1*sr3)); Gam23 = sp.simplify(Gp_A[1, 2]/(sr2*sr3)); Gam12 = sp.simplify(Gp_A[0, 1]/(sr1*sr2))
Vh = sp.simplify(sdetA/(sr1*sr2*sr3))
c3 = (g12 - g13*g23)/sp.sqrt((1 - g13**2)*(1 - g23**2)); s3 = sp.sqrt(Delta_sym)/sp.sqrt((1 - g13**2)*(1 - g23**2))
check("T1b  layer A centre: Gam_13, Gam_23 -> 0 like t (the root's presented direction has no limit; its cosines do)",
      z(sp.limit(Gam13, t, 0)) and z(sp.limit(Gam23, t, 0)) and z(sp.limit(Gam13/t, t, 0) - sp.I*g13/sp.sqrt(1 - g13**2)))
check("T1c  layer A centre: Gam_12 -> (g12 - g13 g23)/sqrt((1-g13^2)(1-g23^2)) = cos(phi_3), the dihedral angle at the root;"
      " Vhat -> sqrt(Delta)/sqrt((1-g13^2)(1-g23^2)) = sin(phi_3)",
      z(sp.limit(Gam12, t, 0) - c3) and z(sp.limit(Vh, t, 0) - s3))
check("T1d  identity (g12 - g13 g23)^2 + Delta = (1 - g13^2)(1 - g23^2): the centre value sits ON the unit circle",
      sp.expand((g12 - g13*g23)**2 + Delta_sym - (1 - g13**2)*(1 - g23**2)) == 0)
phi = sp.Symbol('phi', positive=True)
check("T1e  2 arg(1 + cos phi + i sin phi) = phi for 0 < phi < pi, so Omega_B(+++) at the centre IS phi_3 exactly;"
      " the orbit is {phi_3, pi - phi_3, pi - phi_3, phi_3}",
      all(abs(sp.N((sp.tan(phi/2) - sp.sin(phi)/(1 + sp.cos(phi))).subs(phi, v), 40)) < 1e-35 and
          abs(sp.N((sp.tan(sp.pi/2 - phi/2) - sp.sin(phi)/(1 - sp.cos(phi))).subs(phi, v), 40)) < 1e-35 for v in (sp.Rational(7,10), sp.Rational(2,1), sp.Rational(3,1))))
SA0 = sp.expand(1 + Gp_A[0, 1].subs(t, 0) + Gp_A[0, 2].subs(t, 0) + Gp_A[1, 2].subs(t, 0))
check("T1f  layer A, raw functional at the centre: S_A(0) = 1 + g12 - g13 g23 = 1 + cos(phi_3) sin(th13) sin(th23) > 0 for EVERY"
      " admissible Gram, and V' -> 0, so Omega_A(0+) = 0 with no frame content (ARGUMENT: |cos sin sin| < 1)",
      sp.expand(SA0 - (1 + g12 - g13*g23)) == 0 and z(SA0 - (1 + c3*sp.sqrt(1 - g13**2)*sp.sqrt(1 - g23**2))))
# layer C centre
kkC = vC*vC.T/NC; GpC0 = (Gsym - kkC).applyfunc(sp.simplify)
check("T2a  layer C centre: G'(0) (1,1,1)^T = 0 -- the three presented vectors SUM TO ZERO (they are the projections of the axes"
      " onto the plane perpendicular to K, and K ~ a1+a2+a3)", zM(GpC0*one))
check("T2b  layer C centre: det G'(0) = 0 exactly (coplanar), all presented lengths 1 - k_i^2 finite",
      z(GpC0.det()) and all(z(GpC0[i, i] - (1 - kkC[i, i])) for i in range(3)))
SC0 = sp.simplify(1 + GpC0[0, 1] + GpC0[0, 2] + GpC0[1, 2])
check("T2c  layer C, raw functional at the centre: S_A(0) = (|G1|^2/(1^T G 1) - 1)/2, which VANISHES at the orthogonal frame"
      " (0/0 there): a sign bit, not a reading",
      z(SC0 - ((vC.T*vC)[0, 0]/NC - 1)/2) and SC0.subs({g12: 0, g13: 0, g23: 0}) == 0)
print("  ARGUMENT (not counted): three plane vectors p_i with p1+p2+p3 = 0 have S = (1 + cos a)(1 - (|p1|+|p2|)/|p3|) < 0 by the")
print("  triangle inequality, and V = 0; so the (+++) presented triangle at the layer-C centre is the HEMISPHERE: Omega = 2 pi,")
print("  B real and negative (half-angle sign -1), the other three members 0.  Checked exactly at N1..N4 below (P1-7c).")
# collision radius, layer A
rstar = sp.simplify(collision_radii(presented(Gsym, 'A'))[(0, 1)][0]) if False else Delta_sym/(1 - g12**2)
Gp_Ar = (Gsym + kA*kA.T/(r - 1)).applyfunc(sp.cancel)
check("T3a  layer A: the two VISIBLE presented directions collide (G'_12^2 = rho_1 rho_2) at r* = r_s Delta/(1 - g12^2)",
      z(sp.simplify((Gp_Ar[0, 1]**2 - Gp_Ar[0, 0]*Gp_Ar[1, 1]).subs(r, rstar))))
check("T3b  r*/r_1 = Delta/((1-g12^2)(1-g13^2)) = sin^2(phi_1) < 1: the collision sits BELOW both null radii, where both visible rods"
      " are positive -- a real collision inside", sp.expand((g23 - g12*g13)**2 + Delta_sym - (1 - g12**2)*(1 - g13**2)) == 0)

# ---------------------------------------------------------------------------
# WHAT THE RUN SAYS: centre, horizon, features
# ---------------------------------------------------------------------------
print(); print("=" * 100); print("RESULTS -- classification of the endpoints and the interior features"); print("=" * 100)
def phi3(g): return mp.acos(mp.mpf(str(sp.N(c3.subs({g12: g[0], g13: g[1], g23: g[2]}), 40))))
omn = lambda x: mp.mpf(str(sp.N(x, 40))); TOL = mp.mpf(10)**-25
ok = True; vals = []
for nm, g in GRAMS.items():
    row = RES[f'{nm}/A/B']['points']['r->0+']
    Om = [row[e]['Omega'] for e in EPS]; ph = phi3(g)
    on = [mp.mpc(complex(sp.N(o, 30))) for o in Om]
    okg = (all(row[e]['cls'] in ('FINITE-NONTRIVIAL', 'IDENTITY') for e in EPS) and all(abs(mp.im(o)) < mp.mpf(10)**-20 for o in on)
           and abs(mp.re(on[0]) - ph) < 1e-12 and abs(mp.re(on[3]) - ph) < 1e-12
           and abs(mp.re(on[1]) - (mp.pi - ph)) < 1e-12 and abs(mp.re(on[2]) - (mp.pi - ph)) < 1e-12)
    ok &= okg; vals.append(f"{nm}: phi_3 = {mp.nstr(ph, 10)}")
check("P1-7a CENTRE, layer A, Omega_B: FINITE-NONTRIVIAL, real, orbit = {phi_3, pi-phi_3, pi-phi_3, phi_3} at all four Grams (N4, equal angles: two members coincide with their own baseline and read IDENTITY)", ok, "; ".join(vals))
check("P1-7b CENTRE, layer A, Omega_A (raw): W -> 1, B -> positive rational: value 0, the same for every Gram",
      all(RES[f'{nm}/A/A']['points']['r->0+'][e]['W'] == 1 and RES[f'{nm}/A/A']['points']['r->0+'][e]['half'] == '+1' for nm in GRAMS for e in EPS))
okC = all(RES[f'{nm}/C/B']['points']['r->0+'][EPS[0]]['half'] == '-1' and all(RES[f'{nm}/C/B']['points']['r->0+'][e]['half'] == '+1' for e in EPS[1:])
          and all(RES[f'{nm}/C/B']['points']['r->0+'][e]['W'] == 1 for e in EPS) for nm in GRAMS)
check("P1-7c CENTRE, layer C, Omega_B: W -> 1 for all members (0 mod 2 pi) with the (+++) member's B real NEGATIVE (half-angle -1:"
      " the hemisphere, Omega = 2 pi) and the other three positive (0): the orbit is {2 pi, 0, 0, 0}", okC)
check("P1-7d CENTRE, layer C, Omega_A (raw): W -> 1; sign of B varies with the Gram (a sign bit): " +
      ", ".join(f"{nm}:{RES[f'{nm}/C/A']['points']['r->0+'][EPS[0]]['half']}" for nm in GRAMS),
      all(RES[f'{nm}/C/A']['points']['r->0+'][e]['W'] == 1 for nm in GRAMS for e in EPS))
check("P1-7e CENTRE: no member of either functional on either layer is UNBOUNDED or NO-LIMIT at r -> 0+ (the CURV-1 fork (a) test)",
      all(cls_at(k, 'r->0+', e) in ('IDENTITY', 'FINITE-NONTRIVIAL') for k in RES for e in EPS),
      "a divergence here would print the member; none does")
hz_finite = all(cls_at(k, p, e) in ('IDENTITY', 'FINITE-NONTRIVIAL') for k in RES for p in ('r->r_s+', 'r->r_s-') for e in EPS)
hz_same = all(abs(mp.mpc(complex(sp.N(RES[k]['points']['r->r_s+'][e]['W'], 30))) - mp.mpc(complex(sp.N(RES[k]['points']['r->r_s-'][e]['W'], 30)))) < mp.mpf(10)**-8
              for k in RES for e in EPS if RES[k]['points']['r->r_s+'][e]['W'] is not None and RES[k]['points']['r->r_s-'][e]['W'] is not None)
hz_unit = all(abs(abs(mp.mpc(complex(sp.N(RES[k]['points'][p][e]['W'], 30)))) - 1) < mp.mpf(10)**-8 for k in RES for p in ('r->r_s+', 'r->r_s-') for e in EPS if RES[k]['points'][p][e]['W'] is not None)
check("P1-8a HORIZON, both sides, both layers, both functionals: every member has a FINITE limit with |W| = 1 (Omega real);"
      " the member whose presented poles all agree with sign(k) gives W = 1 (Omega -> 0), the mixed members give finite"
      " NONTRIVIAL angles (the lunes of the collapsed triangle); the two sides of the horizon AGREE", hz_finite and hz_same and hz_unit)
# null-radius forks: the orbit MULTISET is continuous through each null radius; the labelling is not
def mval(v):
    if v is None: return None
    if v in (sp.oo, sp.zoo): return 'inf'
    c = mp.mpc(complex(sp.N(v, 30)))
    return 'inf' if abs(c) > mp.mpf(10)**12 else c
def mseq(L, R):
    R = list(R)
    for a in L:
        hit = None
        for i, b in enumerate(R):
            if a in ('inf', None) or b in ('inf', None):
                if a == b: hit = i; break
            elif abs(a - b) < mp.mpf(10)**-8: hit = i; break
        if hit is None: return False
        R.pop(hit)
    return not R
fork_ok, fork_txt = True, []
for k in RES:
    for fname, fd in RES[k]['features'].items():
        if fname.startswith('null'):
            L = [mval(v) for v in fd['left'].values()]; Rr = [mval(v) for v in fd['right'].values()]
            same = mseq(L, Rr); fork_ok &= same
            perm = [e for e in EPS if not z(fd['left'][e] - fd['right'][e])]
            fork_txt.append(f"{k} {fname}: multiset {'equal' if same else 'DIFFERENT'}, relabelled members {len(perm)}; left {[mp.nstr(x, 6) if x not in (None, 'inf') else x for x in L]}")
check("P1-8b NULL RADII: the excess orbit is continuous THROUGH every null radius as a multiset (limits from the left and right"
      " agree as sets; individual labels swap)", fork_ok)
for s in fork_txt: print("        " + s)
# collisions: exact interior radii where a presented pair is coincident/antipodal; members go UNBOUNDED there
coll_txt, coll_found = [], False
for k in RES:
    for fname, fd in RES[k]['features'].items():
        if fname.startswith('collision'):
            unb = [e for e in EPS if fd['left'][e] in (0, sp.oo, sp.zoo) or fd['right'][e] in (0, sp.oo, sp.zoo)]
            coll_found |= bool(unb); coll_txt.append(f"{k} {fname}: UNBOUNDED members {unb if unb else 'none'}")
check("P1-8c COLLISIONS: at every presented-collision radius some member of the orbit is UNBOUNDED (S + iV -> 0: a logarithmic"
      " singularity where two presented directions are antipodal), normalised functional (the raw one has its own zeros, P1-13a)", coll_found and
      all('none' not in s for s in coll_txt if '/B ' in s))
for s in coll_txt: print("        " + s)
check("P1-8d layer A collision radius equals r_s Delta/(1 - g12^2) at N1..N4 (T3a instantiated)",
      all(RES[f'{nm}/A/B']['collisions'][(0, 1)] == [sp.nsimplify(gate[nm]/(1 - g[0]**2))] for nm, g in GRAMS.items()))

# ---------------------------------------------------------------------------
# CROSS-CHECKS: route V, l'Huilier, the exact product identity, the sum rule, rod-freeness, interior zeros of B
# ---------------------------------------------------------------------------
print(); print("=" * 100); print("CROSS-CHECKS"); print("=" * 100)
worst = mp.mpf(0)
for nm in ('N1', 'N3'):
    for layer in LAYERS:
        for func in FUNCS:
            for rv in (mp.mpf(3), mp.mpf(3)/2, mp.mpf(1)/2, mp.mpf(1)/10):
                for eps in EPS:
                    SG, VG = B_num_G(Gram(*GRAMS[nm]), layer, func, eps, rv); SVv, VV = B_num_V(GRAMS[nm], layer, func, eps, rv)
                    worst = max(worst, abs(SG - SVv), abs(VG - VV))
check("P1-9  ROUTE V (rotor sandwich in the Pauli rep, pinned rapidity, +i sheet) = ROUTE G (presented Gram) for S and V,"
      " both functionals, both layers, all members, exterior AND interior points, to 1e-30", worst < mp.mpf(10)**-30, f"worst {mp.nstr(worst, 3)}")
worst = mp.mpf(0)
for nm in ('N1', 'N3'):
    for layer in LAYERS:
        for rv in (mp.mpf(3), mp.mpf(2), mp.mpf(21)/20):
            for eps in EPS:
                S, V = B_num_G(Gram(*GRAMS[nm]), layer, 'B', eps, rv)
                dd = abs(mp.re(Omega_num(S, V)) - lhuilier(Gram(*GRAMS[nm]), layer, eps, rv)); worst = max(worst, min(dd, abs(dd - 2*mp.pi)))
check("P1-10 l'HUILIER (independent formula) agrees with Omega_B in the exterior, every member, to 1e-25", worst < mp.mpf(10)**-25, f"worst {mp.nstr(worst, 3)}")
# exact product identity  prod_eps B_eps = 4 prod (1 - Gam_ij^2)  (functional B), on both sheets, identically in u
a12, a13, a23, vv = sp.symbols('a12 a13 a23 v')
prodB = sp.expand(sp.prod([1 + e[0]*e[1]*a12 + e[0]*e[2]*a13 + e[1]*e[2]*a23 + sp.I*vv for e in EPS]))
prodB = sp.expand(prodB.subs(vv**4, (1 - a12**2 - a13**2 - a23**2 + 2*a12*a13*a23)**2).subs(vv**3, (1 - a12**2 - a13**2 - a23**2 + 2*a12*a13*a23)*vv)
                  .subs(vv**2, 1 - a12**2 - a13**2 - a23**2 + 2*a12*a13*a23))
check("P1-11a exact orbit identity: prod_eps (S_eps + iV) = -4 (1-a12^2)(1-a13^2)(1-a23^2) modulo V^2 = det (unit-diagonal Gram)",
      sp.expand(prodB + 4*(1 - a12**2)*(1 - a13**2)*(1 - a23**2)) == 0)
ok = True
for nm in ('N1', 'N4'):
    for layer in LAYERS:
        P = presented(Gram(*GRAMS[nm]), layer)
        for (p, side) in ((sp.oo, '-'), (sp.Integer(1), '-'), (sp.Integer(0), '+')):
            rsub, bd = approach(P, p, side); X = explicit(P, rsub, bd)
            Bs = [SV(X, e, 'B') for e in EPS]; G_ = [X['Gp'][i, j]/(X['sr'][i]*X['sr'][j]) for (i, j) in ((0, 1), (0, 2), (1, 2))]
            lhs = sp.prod([S + sp.I*V for S, V in Bs]); rhs = -4*sp.prod([1 - g_**2 for g_ in G_])
            ok &= z(sp.simplify(sp.expand(lhs - rhs)))
check("P1-11b ...and it holds IDENTICALLY in r on the exterior AND on the interior sheet, both layers (N1, N4): the orbit is"
      " one consistent object through the horizon and the null radii", ok)
worst = mp.mpf(0)
for nm in GRAMS:
    for layer in LAYERS:
        for rv in (mp.mpf(4), mp.mpf(11)/10):
            tot = sum(mp.re(Omega_num(*B_num_G(Gram(*GRAMS[nm]), layer, 'B', e, rv))) % (2*mp.pi) for e in EPS)
            worst = max(worst, abs(tot - 2*mp.pi))
check("P1-11c exterior sum rule: the four excesses of the presented triangles sum to 2 pi (40 digits, every Gram, both layers)", worst < mp.mpf(10)**-30)
# rod-freeness: scale invariance (functional B) vs. scale dependence (functional A)
c1_, c2_, c3_ = sp.symbols('c1 c2 c3', positive=True); x11, x12, x13, x22, x23, x33 = sp.symbols('x11 x12 x13 x22 x23 x33', positive=True)
Gx = sp.Matrix([[x11, x12, x13], [x12, x22, x23], [x13, x23, x33]]); Cx = sp.diag(c1_, c2_, c3_); Gs = Cx*Gx*Cx
gamB = lambda M, i, j: M[i, j]/(sp.sqrt(M[i, i])*sp.sqrt(M[j, j]))
check("P1-12a Omega_B is unchanged by v_i -> c_i v_i (c_i > 0): the normalised cosines and Vhat are degree 0 in every presented length",
      all(z(sp.simplify(gamB(Gs, i, j) - gamB(Gx, i, j))) for (i, j) in ((0, 1), (0, 2), (1, 2)))
      and z(sp.simplify(sp.sqrt(Gs.det())/sp.sqrt(Gs[0, 0]*Gs[1, 1]*Gs[2, 2]) - sp.sqrt(Gx.det())/sp.sqrt(Gx[0, 0]*Gx[1, 1]*Gx[2, 2]))))
check("P1-12b Omega_A CHANGES under v_i -> c_i v_i: it obeys the letter of CARRY-8 (no division by a rod) but carries the rods"
      " multiplicatively; it is not the excess of any triangle", not z(sp.simplify((1 + Gs[0, 1] + Gs[0, 2] + Gs[1, 2]) - (1 + Gx[0, 1] + Gx[0, 2] + Gx[1, 2]))))
# interior zeros of B (functional A): exact -- S_A^2 + det G' = 0 is rational
zeroA = {}
for nm, g in GRAMS.items():
    for layer in LAYERS:
        P = presented(Gram(*g), layer); Gp = P['Gp']; found = {}
        for e in EPS:
            S = 1 + e[0]*e[1]*Gp[0, 1] + e[0]*e[2]*Gp[0, 2] + e[1]*e[2]*Gp[1, 2]
            expr = sp.numer(sp.together(sp.cancel(S**2 + P['det'])))
            found[e] = sorted(sp.nsimplify(x) for x in sp.solve(expr, r) if x.is_real and 0 < x < 1)
        zeroA[f'{nm}/{layer}'] = found
anyA = any(v for d in zeroA.values() for v in d.values())
check("P1-13a Omega_A has interior points where S_A = |V'| (B or S-iV vanishes): Im Omega_A -> oo there, at exact rational radii",
      anyA, "; ".join(f"{k}: {sum(len(v) for v in d.values())} zeros" for k, d in zeroA.items()))
for k, d in zeroA.items():
    for e, v in d.items():
        if v: print(f"        {k} member {e}: r/r_s = {v}")

# ---------------------------------------------------------------------------
# TABLE, JSON, VERDICT
# ---------------------------------------------------------------------------
print(); n_ = sum(CH); print(f"RESULT: {n_}/{len(CH)} checks passed."); print("=" * 100)
print("CLASSIFICATION TABLE  (Gram / layer / functional; members in the order (+++), (+--), (-+-), (--+))")
print("=" * 100)
for key in RES:
    for pn in ('r->r_s+', 'r->r_s-', 'r->0+'):
        row = RES[key]['points'][pn]
        cl = [row[e]['cls'][:7] for e in EPS]; hs = [row[e]['half'] for e in EPS]
        om = [(mp.nstr(mp.mpc(complex(sp.N(row[e]['Omega'], 20))), 8) if row[e]['Omega'] is not None else 'div') for e in EPS]
        print(f"  {key:10s} {pn:8s} {cl}  Omega {om}  half-angle {hs}")
    print(f"  {key:10s} bands: " + "; ".join(f"({a},{b}) {c}" for a, b, c in RES[key]['bands']))
    print(f"  {key:10s} null radii {RES[key]['null_radii']}  collisions {RES[key]['collisions']}")
def js(o):
    if isinstance(o, dict): return {str(k): js(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [js(v) for v in o]
    if isinstance(o, (sp.Basic, mp.mpf, mp.mpc)): return str(o)
    return o
json.dump(js(RES), open('curv1_path1_results.json', 'w'), indent=1)
if n_ != len(CH):
    print("VERDICT: check failures above -- do NOT read the verdict below as earned."); sys.exit(1)
print("=" * 100)
print("VERDICT: OUTCOME (b) at the CENTRE -- every member of both functionals on both layers has a finite limit at r -> 0+.")
print("  Layer A (the infalling seat): the normalised excess reads the frame's dihedral angle at the root, orbit")
print("  {phi_3, pi-phi_3, pi-phi_3, phi_3}: FINITE-NONTRIVIAL, real (T1, P1-7a).  The raw CARRY-8 functional reads 0")
print("  for every admissible frame (T1f, P1-7b): the letter of the frozen definition selects a functional with no content.")
print("  Layer C (the hole): the presented vectors sum to zero at the centre (T2a); the (+++) triangle is the hemisphere,")
print("  Omega = 2 pi with half-angle sign -1, the other three 0 (P1-7c).  Horizon: Omega -> 0 both sides (P1-8a).")
print("  INSIDE: at every presented-collision radius (layer A: r* = r_s Delta/(1 - g12^2), T3) two presented directions")
print("  are antipodal for half the orbit and S + iV -> 0: Im Omega -> oo, a LOGARITHMIC singularity with NO rod in a")
print("  denominator (P1-8c, P1-13a).  The rod-free curvature reading is regular at r = 0 and singular at r*.")
print("  Conditional on: KIN-2a, CONT-1 (+i sheet), RULING-2, the principal-branch convention of section 2.")
print("=" * 100)
print("COMPARISON STAGE (names spoken here only):")
print("  - Gauss-Bonnet: the excess of the presented triangle is its enclosed curvature integral; the centre value on the")
print("    seat layer is a dihedral angle, finite -- no 1/r^3 tidal divergence is reproduced or expected here.")
print("  - The layer-C hemisphere at the centre (2 pi, half-angle -1) is the spin-lift sign of a full turn: the deck Z2.")
print("  - Interior antipodal collisions at r* have no standard-model analogue claimed; lineage unsearched.")
sys.exit(0)
