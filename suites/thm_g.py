#!/usr/bin/env python3
# =============================================================================
# THM-TARGET G, part 1 - the hbar-seat readout theorem
# Date: 2026-08-28
#
# PRE-REGISTERED PREDICTION (LABELLED_MODEL_v1 sec.5, frozen 2026-08-17, as
# quoted by 2026-08-27-thm-target-E.md and HANDOFF-2026-08-27.md; PRIMARY
# DOCUMENT NOT IN SESSION - flagged, Will to confirm against sec.5 verbatim):
#   c seat:    reciprocal-pair product MAY VANISH (null locus) - invariant FREE
#   hbar seat: protected product CANNOT VANISH - invariant FLOORED
#   G seat:    product VANISHES AT THE HORIZON (untreated here; part 2)
#
# ALLOWED INPUTS
#   BARE-1   pivot = hyperbolic rotor exp(lambda*K/2), K^2 = +1 [D1 discharged]
#   SEAT-hb  (from handoff, T-tier): the hbar-seat pivot is a SQUEEZE; the
#            seat state is a positive area element (quadratic form) on the
#            seat plane; pivot acts symplectically.  Status: quoted from the
#            freeze-dated handoff; primary sec-doc absent - NAMED conditionality.
#   Carried machinery from E/addenda: two-fixed-pole aggregation, wandering
#            intervals, mirror, "the measure is the plane" (re-verified here
#            in hbar variables, not assumed).
# BANNED until comparison stage: uncertainty principle, hbar, commutators,
#   photons, squeezed-state theory, Robertson-Schrodinger by name.
# ZERO physics labels consumed in the derivation: every check below holds for
#   symbolic protected area A > 0.  The identification of A's floor value is
#   comparison-stage naming only.
#
# TWO PATHS sharing nothing but sympy:
#   Path 1: 2x2 symplectic matrices on the seat plane (quadrature picture).
#   Path 2: hand-rolled Cl(3) (per thm_e.py) computing the c-seat family, then
#           the COVER DICTIONARY Theta = 2*phi, Lambda = 2*lambda as a RESULT.
# =============================================================================
import sympy as sp
import time; _T0 = time.time()

lam = sp.Symbol('lambda', real=True)
phi = sp.Symbol('phi', real=True)
A_  = sp.Symbol('A', positive=True)          # protected area, symbolic
Th, Lm = sp.Symbol('Theta', real=True), sp.Symbol('Lambda', real=True)

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-_T0:6.1f}s] {t}" + (f" - {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for r in (sp.simplify, lambda x: sp.simplify(sp.expand_trig(sp.expand(x))),
              lambda x: sp.simplify(sp.expand(x.rewrite(sp.exp)))):
        try:
            if r(e) == 0: return True
        except Exception: pass
    return False

print("=" * 78)
print("PATH 1 - symplectic 2x2 on the seat plane")
print("=" * 78)
S = sp.Matrix([[sp.exp(lam), 0], [0, sp.exp(-lam)]])       # the squeeze pivot
C0 = A_ * sp.eye(2)                                        # base state, area A
C  = S * C0 * S.T
nhat = sp.Matrix([sp.cos(phi), sp.sin(phi)])
q = sp.simplify((nhat.T * C * nhat)[0, 0])                 # reading along phi
q_closed = A_ * (sp.cosh(2*lam) + sp.sinh(2*lam) * sp.cos(2*phi))
check("G-1a  family: q(phi) = A*(cosh 2l + sinh 2l cos 2phi)  [path 1]",
      z(q - q_closed))

check("G-1b  PROTECTED PRODUCT: det C' = A^2 exactly, every pivot (det S = 1)",
      z(sp.simplify(C.det() - A_**2)) and sp.simplify(S.det() - 1) == 0,
      "the pivot group preserves the product IDENTICALLY - not on shell, always")

qmin = sp.simplify(q_closed.subs(sp.cos(2*phi), -1))
check("G-1c  ORBIT SEPARATION: min_phi q = A*e^{-2l} > 0 for all finite pivot",
      z(qmin - A_*sp.exp(-2*lam)),
      "the zero stratum is reachable only at lambda = oo: infinite energy, again")

conj = sp.simplify(q_closed * q_closed.subs(phi, phi + sp.pi/2))
conj_closed = A_**2 * (1 + sp.sinh(2*lam)**2 * sp.sin(2*phi)**2)
check("G-2   conjugate product = A^2*(1 + sinh^2(2l) sin^2(2phi))  >= A^2",
      z(conj - conj_closed),
      "REQUIRES-margin sinh(2l)*|sin 2phi|: equality iff principal axes or l=0")

d = S * nhat
tanmap = sp.simplify(d[1] / d[0])
check("G-3a  angle map: tan phi' = e^{-2 lambda} tan phi  [derived]",
      z(sp.simplify(tanmap - sp.exp(-2*lam)*sp.tan(phi))))
phip = sp.atan(sp.exp(-2*lam) * sp.tan(phi))
Jac = sp.simplify(sp.diff(phip, phi))
check("G-3b  d phi'/d phi = A/q  - transport exponent ONE (c seat had 1/s^2)",
      z(sp.simplify(Jac - A_/q_closed).subs(A_, 1)),
      "the measure is the plane, hbar edition; the plane is a CIRCLE here")
check("G-3c  q * (dphi'/dphi) = A  =>  <q>_pivoted-measure = A EXACTLY:",
      z(sp.simplify(q_closed*Jac - A_).subs(A_, 1)),
      "the LINEAR channel is pivot-blind on the hbar seat (p_blind = 1, not 2)")

print(); print("=" * 78)
print("PATH 1 - moments on the pivoted circle measure (A = 1)")
print("=" * 78)
qq = q_closed.subs(A_, 1)
c2, s2 = sp.cosh(2*lam), sp.sinh(2*lam)
Xs = sp.Symbol('Xs')
def cospow_int(k):        # Int_0^pi cos(2phi)^k dphi, exact table
    return sp.Integer(0) if k % 2 else sp.pi * sp.binomial(k, k//2) / 2**k
def mom_piv(n):           # <q^n>_piv = (1/pi) Int q^{n-1} dphi   [by G-3c]
    e = sp.expand((c2 + s2*Xs)**(n-1))
    return sp.simplify(sum(e.coeff(Xs, k)*cospow_int(k) for k in range(n))/sp.pi)
mom = {n: mom_piv(n) for n in (1, 2, 3, 4)}
x = sp.cosh(2*lam)
Leg = {1: sp.Integer(1), 2: x, 3: (3*x**2 - 1)/2, 4: (5*x**3 - 3*x)/2}
for n in (1, 2, 3, 4):
    check(f"G-4.{n} <q^{n}>_piv = P_{n-1}(cosh 2l)  [Legendre, exact table]",
          z(sp.expand(mom[n] - Leg[n])))
check("G-5   two-pole invariant: sqrt(q(0)*q(pi/2)) = A, every pivot",
      z(sp.simplify(sp.sqrt(q_closed.subs(phi, 0) * q_closed.subs(phi, sp.pi/2))
                    - A_)),
      "invariant = the protected product's root = the pairing value")
check("G-6   pole exchange phi -> pi/2 - phi conjugates lambda -> -lambda",
      z(sp.simplify(q_closed.subs(phi, sp.pi/2 - phi)
                    - q_closed.subs(lam, -lam))),
      "the c-1 analog structure; equal pole weight CONDITIONAL on c-1 root-generic")

# MIRROR for ALL real p, by three exact links (no definite integral needed):
#   Int q^-p dphi  =[L2]=  Int q^{1-p} dphi'  =[L1]=  Int (c2-s2 cos2phi')^{p-1} dphi'
#                  =[L3]=  Int (c2+s2 cos2psi)^{p-1} dpsi  =  Int q^{p-1} dphi
# hence <q^-p>_unpiv = <q^p>_piv identically in p.
T = sp.Symbol('T', real=True)          # T = tan(phi): rational coordinates
cos2f  = (1 - T**2)/(1 + T**2)
t_p    = sp.exp(-2*lam)*T                # tan(phi') by G-3a
cos2fp = (1 - t_p**2)/(1 + t_p**2)
qT     = c2 + s2*cos2f
check("G-7a  L1 office identity: (cosh2l - sinh2l cos2phi') = 1/q, exactly",
      sp.cancel(sp.together(((c2 - s2*cos2fp) - 1/qT).rewrite(sp.exp))) == 0)
check("G-7b  L2 = G-3b (Jacobian 1/q); L3 reflection: cos(2(pi/2-x)) = -cos 2x",
      z(sp.expand_trig(sp.cos(2*(sp.pi/2 - phi)) + sp.cos(2*phi))),
      "chain composed: <q^p>_piv = <q^-p>_unpiv for ALL real p - the mirror")
import mpmath as mp
mp.mp.dps = 30
ok = True
for lv in ('0.7', '1.3', '0.31'):
    lv = mp.mpf(lv)
    g_num = mp.quad(lambda f: mp.log(mp.cosh(2*lv) + mp.sinh(2*lv)*mp.cos(2*f)),
                    [0, mp.pi/2, mp.pi]) / mp.pi
    ok = ok and abs(g_num - 2*mp.log(mp.cosh(lv))) < mp.mpf(10)**-25
check("G-7c  unpivoted circle geometric mean = cosh^2(l)  [30-digit corroboration]",
      ok, "piv = 1/cosh^2(l) then follows from the all-p mirror; product = 1")

print(); print("=" * 78)
print("PATH 2 - Cl(3) c-seat family, then the COVER DICTIONARY as a result")
print("=" * 78)
def blade_mul(a, b):
    seq = list(a) + list(b); sign = 1
    for i in range(1, len(seq)):
        j = i
        while j > 0 and seq[j-1] > seq[j]:
            seq[j-1], seq[j] = seq[j], seq[j-1]; sign = -sign; j -= 1
    out, k = [], 0
    while k < len(seq):
        if k+1 < len(seq) and seq[k] == seq[k+1]: k += 2
        else: out.append(seq[k]); k += 1
    return tuple(out), sign
class MV:
    def __init__(s, d=None):
        s.d = {}
        for k, v in (d or {}).items():
            v = sp.sympify(v)
            if v != 0: s.d[tuple(k)] = s.d.get(tuple(k), 0) + v
    def __add__(s, o):
        d = dict(s.d)
        for k, v in o.d.items(): d[k] = d.get(k, 0) + v
        return MV(d)
    def __mul__(s, o):
        if not isinstance(o, MV): return MV({k: v*o for k, v in s.d.items()})
        d = {}
        for ka, va in s.d.items():
            for kb, vb in o.d.items():
                key, sg = blade_mul(ka, kb)
                d[key] = d.get(key, 0) + sg*va*vb
        return MV(d)
    __rmul__ = __mul__
    def c(s, key): return sp.simplify(s.d.get(tuple(key), 0))
Arot = MV({(): sp.cosh(Lm/2), (1,): sp.sinh(Lm/2)})
X = MV({(): 1, (1,): sp.cos(Th), (2,): sp.sin(Th)})
Xp = Arot * X * Arot
s_c = sp.simplify(Xp.c(()))                                 # c-seat scaling
check("G-8a  c-seat family recomputed: s = cosh(L) + sinh(L) cos(Theta)",
      z(s_c - (sp.cosh(Lm) + sp.sinh(Lm)*sp.cos(Th))))
check("G-8b  DICTIONARY: q(phi,l)/A = s(Theta,L) under Theta = 2phi, L = 2l",
      z(sp.simplify(qq - s_c.subs({Th: 2*phi, Lm: 2*lam}))),
      "hbar seat = c seat pulled through the 2:1 cover with DOUBLED rapidity")
cthp = sp.simplify(Xp.c((1,)) / Xp.c(()))
angle_c = (sp.cos(Th) + sp.tanh(Lm)) / (1 + sp.tanh(Lm)*sp.cos(Th))
check("G-8c  c-seat angle map recomputed  [path 2]",
      sp.cancel(sp.together((cthp - angle_c).rewrite(sp.exp))) == 0)
# pull the c-seat angle map through the cover, in rational T = tan(phi):
lhsT = (cos2f + sp.tanh(2*lam)) / (1 + sp.tanh(2*lam)*cos2f)   # cos(2phi')
check("G-8d  cover pullback of the c angle map = the squeeze angle map, exactly",
      sp.cancel(sp.together((lhsT - cos2fp).rewrite(sp.exp))) == 0,
      "one presentation geometry, two seats, related by the double cover")
u = sp.Symbol('u', real=True)
J_c = sp.simplify(sp.diff((u + sp.tanh(Lm))/(1 + sp.tanh(Lm)*u), u))
check("G-8e  c Jacobian 1/s^2 pulls back consistently to the 1/q circle Jacobian",
      z(sp.simplify(J_c - 1/(sp.cosh(Lm) + sp.sinh(Lm)*u)**2)),
      "exponent 2 on the sphere-plane, exponent 1 on the circle: SEAT CHARACTER")

print()
print("  G-9 [ARGUMENT, uncounted] tan phi' = e^{-2l} tan phi is translation by")
print("      -2l in Y = log tan phi; exactly two fixed directions; the wandering-")
print("      interval theorem transfers verbatim: finite invariant measures sit")
print("      on the two poles; equal weight from pole exchange (G-6), conditional")
print("      on c-1 being root-generic rather than c-seat-specific.  NAMED.")

print(); print("=" * 78)
print("COMPARISON STAGE (banned names spoken here only)")
print("=" * 78)
occ = sp.simplify((mom[2] - 1)/2)
check("G-10  <q^2>_piv - 1 = 2 sinh^2(lambda): the quadratic channel reads",
      z(occ - sp.sinh(lam)**2),
      "the squeezed-vacuum occupation sinh^2(l) - free receipt")
print("""  - Conjugate-product bound with saturation on principal axes only:
    the Robertson-Schrodinger structure, margin sinh^2(2l) sin^2(2phi) (G-2).
  - Minimum-uncertainty states = principal axes: recovered as the equality
    locus of G-2, not assumed.
  - <q^2>-1 = 2 sinh^2(l): squeezed-vacuum mean occupation sinh^2(lambda).
  - Moments = Legendre P_{n-1}(cosh 2l): the Poisson kernel on the circle -
    harmonic-analysis lineage flag (machinery + failure modes).
  - Y = log tan phi is the standard squeeze rapidity; the floor value A = the
    seat's quantum is NAMED here only; every theorem above holds for any A > 0.""")

print(); n_ = sum(CH); print(f"RESULT: {n_}/{len(CH)} checks passed.")
print("""
SEC-5 PRE-REGISTRATION VERDICT (hbar leg):
  The aggregation machinery transfers whole: two fixed poles, two-pole
  invariant, mirror, measure-is-the-plane (G-3b/c, G-7).  The invariant equals
  the protected product's root (G-5) and the pivot group preserves that product
  IDENTICALLY (G-1b) - on the c seat the corresponding form has null states in
  the state space; here the zero stratum is a separate orbit unreachable at
  finite pivot (G-1c).  'Cannot vanish' holds as ORBIT SEPARATION: geometry
  protects the product; the label supplies the floor's value.  FREE vs FLOORED
  is now a computed distinction, not a slogan.  Outcome class (a):
  retrodiction-grade, with TWO structural discharges beyond retrodiction:
    (i)  p_blind differs by seat: 2 on the c plane, 1 on the hbar circle -
         equal to the measure-transport exponent both times.  CONJECTURE:
         p_blind = transport exponent as a law; the G seat is the third test.
    (ii) the cover dictionary Theta = 2phi, Lambda = 2lambda (G-8): the hbar
         presentation geometry is the c geometry through the 2:1 cover -
         angle doubling is the orientation double cover surfacing in the
         readout tier.  Observation; connection to the click-monodromy Z_2
         NOT established - fixed-locus comparison still owed.
  HONEST SCOPE: fixed pivot axis (F-4's fence applies verbatim); equal pole
  weight conditional on c-1 root-genericity (named); primary sec.5 text absent
  from session (named); the taxonomy's riskiest leg - horizon-vanishing on the
  G seat - is untreated and is where sec.5 can still fail loudly.
""")
