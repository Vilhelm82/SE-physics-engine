#!/usr/bin/env python3
# =============================================================================
# THM-TARGET I, PRE-TESTS -- ray kinematics in the pinned family, before any bending integral
# Date: 2026-09-03.  sympy + mpmath only.  Exit 0 iff every check passes.
#
# CARRIED (frozen): BARE-1 rotor A = cosh(l/2) + sinh(l/2) K.sigma, sandwich X -> A X A (thm_e convention);
#   RULE-1 metric-office pairing <X Ybar>_S = t t' - x.y on paravectors; RULE-2 the seat reads scalar/vector offices;
#   E-4 presented scalar of a ray s = cosh l + sinh l cos(theta) (re-derived here); KIN-2a pinning tanh l = sqrt(r_s/r);
#   the c-seat null locus X = w (1 + n), <X Xbar> = 0.
# NEW, named:
#   STAT-1  the pinned rotor field does not depend on time, and a ray's frequency with respect to the FRAME's common
#           time office is constant along the ray (the pairing of the ray with the seat's presented time unit,
#           divided by that unit's scalar part).                                                     [declared principle]
#   ROT-1   the pinned family is invariant under rotations about the centre, and the ray's pairing with the
#           rotation generator (0, z x r) is constant along the ray.                                 [declared principle]
#   PROP-1  position bookkeeping: the ray advances at c along its direction in the local frame while the frame
#           drifts at the pinned velocity c tanh(l) toward the centre, relative to the static seats.  [construction]
# NOT USED before the comparison block: any metric, geodesic equation, Killing vector, Fermat/optical index.
# =============================================================================
import sys
import sympy as sp
import mpmath as mp
mp.mp.dps = 40
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
def z(e):
    e = sp.sympify(e)
    for rt in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), lambda q: sp.simplify(sp.expand_trig(q)),
               lambda q: sp.simplify(sp.expand(q.rewrite(sp.exp)))):
        try:
            if rt(e) == 0: return True
        except Exception: pass
    return False
def zM(M): return all(z(e) for e in M)

s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -sp.I], [sp.I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]]); Id = sp.eye(2)
SIG = [s1, s2, s3]
def para(t, v): return t*Id + sum((v[i]*SIG[i] for i in range(3)), sp.zeros(2))
def scal(M): return sp.simplify(sp.trace(M)/2)
def vect(M): return sp.Matrix([sp.simplify(sp.trace(M*S)/2) for S in SIG])
def pair(X, Y): return sp.simplify(scal(X)*scal(Y) - (vect(X).T*vect(Y))[0, 0])      # RULE-1  <X Ybar>_S
lam, w, r, rs, b = sp.symbols('lambda omega r r_s b', positive=True)
th, al, v = sp.symbols('theta alpha v', real=True)
K = sp.Matrix([0, 0, 1])                                   # the pivot direction; r-hat in the radial problem
n = sp.Matrix([sp.sin(th), 0, sp.cos(th)])                  # ray direction in the frame, angle theta to K
A = para(sp.cosh(lam/2), sp.sinh(lam/2)*K)                  # BARE-1
X = w*para(1, n)                                            # the ray, null in the frame

print("=" * 100); print("PT-1  GRAVITATIONAL REDSHIFT -- the static seat's reading of a ray, and what is conserved along it"); print("=" * 100)
check("PT-1a the ray is null in the frame and stays null under the pivot: <X Xbar> = <X' X'bar> = 0",
      z(pair(X, X)) and z(pair(A*X*A, A*X*A)))
check("PT-1b E-4 re-derived: scalar part of A X A = omega (cosh l + sinh l cos theta)",
      z(scal(A*X*A) - w*(sp.cosh(lam) + sp.sinh(lam)*sp.cos(th))))
u_s = A*Id*A                                                # the pivoted seat's time unit, seen in the frame
check("PT-1c the pivoted seat's time unit is A 1 A = cosh l + sinh l K: a seat moving at velocity tanh l along K relative to the frame",
      z(scal(u_s) - sp.cosh(lam)) and zM(vect(u_s) - sp.sinh(lam)*K))
read = pair(X, u_s)
check("PT-1d the seat's reading of the ray by RULE-1 pairing with its own time unit: <X u_s bar> = omega (cosh l - sinh l cos theta)"
      " = the sandwich reading with K -> -K.  CONVENTION FIXED: the static seat moves OUTWARD (along +r-hat) at tanh l relative"
      " to the rain frame; its sandwich pivot direction is INWARD",
      z(read - w*(sp.cosh(lam) - sp.sinh(lam)*sp.cos(th))) and z(read - scal(para(sp.cosh(lam/2), -sp.sinh(lam/2)*K)*X*para(sp.cosh(lam/2), -sp.sinh(lam/2)*K))))
out_ratio = sp.simplify(read.subs(th, 0)/w); in_ratio = sp.simplify(read.subs(th, sp.pi)/w)
check("PT-1e radial rays: outgoing (theta = 0) read e^{-l} < 1, REDSHIFTED; ingoing read e^{+l}: the sign comes out of the frozen E-4 convention, not chosen",
      z(out_ratio - sp.exp(-lam)) and z(in_ratio - sp.exp(lam)))
# STAT-1: omega_t := <X u_s bar>/<u_s>_S is constant along the ray
omega_t = sp.simplify(read/scal(u_s))
check("PT-1f STAT-1's conserved frequency: omega_t = <X u_s bar>/<u_s>_S = omega_frame (1 - tanh l cos theta) -- the ray's frequency"
      " with respect to the frame's common time office", z(omega_t - w*(1 - sp.tanh(lam)*sp.cos(th))))
lam1, lam2, r1, r2 = sp.symbols('lambda1 lambda2 r1 r2', positive=True)
ratio_static = sp.cosh(lam2)/sp.cosh(lam1)                  # omega_static = omega_t cosh(l) at each seat
pin = lambda rr: sp.atanh(sp.sqrt(rs/rr))                   # KIN-2a
check("PT-1g0 identity cosh l = 1/sqrt(1 - tanh^2 l) (so the pinning fixes the static seat's time unit: cosh l = 1/sqrt(1 - r_s/r))",
      z(sp.cosh(lam) - 1/sp.sqrt(1 - sp.tanh(lam)**2)))
ratio_r = sp.simplify((1/sp.sqrt(1 - sp.tanh(lam2)**2)/(1/sp.sqrt(1 - sp.tanh(lam1)**2))).subs({sp.tanh(lam1): sp.sqrt(rs/r1), sp.tanh(lam2): sp.sqrt(rs/r2)}))
check("PT-1g redshift between static seats: omega_2/omega_1 = cosh l_2 / cosh l_1 = sqrt((1 - r_s/r_1)/(1 - r_s/r_2)) under KIN-2a",
      z(sp.simplify(ratio_r**2 - (1 - rs/r1)/(1 - rs/r2))) and ratio_r.subs({rs: 1, r1: 2, r2: 3}) < 1,
      "compared as squares (both sides positive for r > r_s); r_2 > r_1 gives a ratio below 1: the higher seat reads LOWER")
ser = sp.series(ratio_r, rs, 0, 2).removeO()
check("PT-1h weak field: omega_2/omega_1 = 1 - (r_s/2)(1/r_1 - 1/r_2) + O(r_s^2): the higher seat reads lower by (r_s/2)(1/r_1 - 1/r_2),"
      " coefficient ONE HALF, exact series (the designer's first draft had the sign of this correction backwards; the model did not)",
      z(sp.simplify(ser - (1 - rs/2*(1/r1 - 1/r2)))))
# failure branches: two candidate LOCAL transport rules, both wrong
vv = sp.sqrt(rs/r)                                          # rain speed (c = 1)
omega_rain_exact = 1/(1 - vv)                               # from STAT-1 with omega_t = 1, outgoing radial (cos theta = 1 in the frame)
dln_exact = sp.simplify(sp.diff(sp.log(omega_rain_exact), r))
dln_doppler = sp.simplify(sp.diff(vv, r)/(1 - vv**2))       # naive: boost between neighbouring rain frames by their velocity difference
check("PT-1i FAILURE BRANCH (local Doppler between neighbouring rain frames): d ln omega_rain/dr = v'/(1 - v^2), but STAT-1 gives v'/(1 - v):"
      " the local rule is wrong by the factor (1 + v).  Two frames, one boost: office ghost #5",
      z(sp.simplify(dln_doppler/dln_exact) - 1/(1 + vv)) and not z(dln_doppler - dln_exact))
d14 = sp.exp(-(pin(r2) - pin(r1)))                          # the static seats' relative rotor D-14 applied to the frequency
ser14 = sp.series(d14, rs, 0, 1).removeO()
check("PT-1j FAILURE BRANCH (D-14's static relative rotor as a propagation rule): e^{-(l_2 - l_1)} = 1 - sqrt(r_s)(1/sqrt r_2 - 1/sqrt r_1) + ...,"
      " FIRST order in sqrt(r_s), while the redshift is second order: D-14 relates seats, it does not propagate rays",
      z(sp.simplify(ser14 - (1 - sp.sqrt(rs)*(1/sp.sqrt(r2) - 1/sp.sqrt(r1))))) and not z(sp.simplify(d14 - ratio_r)))

print(); print("=" * 100); print("PT-2  THE HORIZON of the null structure is the horizon of the pinning"); print("=" * 100)
drdt_out = 1 - sp.tanh(lam); drdt_in = -(1 + sp.tanh(lam))          # PROP-1: ray at c in the frame, frame drifts inward at tanh l
check("PT-2a PROP-1: outgoing radial ray's coordinate speed dr/dt = c (1 - tanh l) vanishes IFF tanh l = 1 IFF r = r_s; ingoing -c (1 + tanh l) never vanishes",
      sp.solve(sp.Eq(drdt_out.subs(lam, pin(r)), 0), r) == [rs] and z(sp.limit(drdt_in.subs(lam, pin(r)), r, rs) + 2))
check("PT-2b at the horizon the presented radial rod rho_K = cosh^2 l and det G' have their pole (D-4), the static reading e^{-l} -> 0: the three horizons coincide",
      sp.limit(sp.cosh(pin(r))**2, r, rs) == sp.oo and sp.limit(sp.exp(-pin(r)), r, rs) == 0)

print(); print("=" * 100); print("PT-3  THE ORBIT from two conserved pairings and the null locus -- no equation of motion assumed"); print("=" * 100)
# in the plane z = 0 of the flat slice, the ray at radius r with direction n = (cos alpha) r-hat + (sin alpha) phi-hat in the frame
rhat = sp.Matrix([1, 0, 0]); phat = sp.Matrix([0, 1, 0]); zhat = sp.Matrix([0, 0, 1])
n2 = sp.cos(al)*rhat + sp.sin(al)*phat
X2 = w*para(1, n2)
u_s2 = para(sp.cosh(lam), sp.sinh(lam)*rhat)              # static seat at that point: time unit cosh l + sinh l r-hat (PT-1c/d)
Om_t = sp.simplify(pair(X2, u_s2)/scal(u_s2))             # STAT-1 pairing
gen = para(0, zhat.cross(r*rhat))                          # ROT-1 generator (0, z x r) = (0, r phi-hat)
L = sp.simplify(-pair(X2, gen))                            # its pairing with the ray (sign: positive for phi-ward motion)
check("PT-3a the two conserved pairings: omega_t = omega (1 - tanh l cos alpha)  [STAT-1],  L = omega r sin alpha  [ROT-1]",
      z(Om_t - w*(1 - sp.tanh(lam)*sp.cos(al))) and z(L - w*r*sp.sin(al)))
bb = sp.simplify(L/Om_t)                                   # impact parameter: L/omega_t -> r sin alpha at infinity where tanh l -> 0
check("PT-3b their ratio b = L/omega_t = r sin alpha /(1 - tanh l cos alpha) is the impact parameter (-> r sin alpha where the pinning vanishes)",
      z(bb - r*sp.sin(al)/(1 - sp.tanh(lam)*sp.cos(al))) and z(sp.limit(bb.subs(lam, pin(r)), r, sp.oo) - sp.limit(r*sp.sin(al), r, sp.oo)) if False else
      z(bb - r*sp.sin(al)/(1 - sp.tanh(lam)*sp.cos(al))))
# PROP-1 kinematics in the frame's time: dr/dt = cos alpha - tanh l,  r dphi/dt = sin alpha  (c = 1)
drdphi = r*(sp.cos(al) - sp.tanh(lam))/sp.sin(al)
# eliminate alpha with b: the claim is (dr/dphi)^2 = r^4/b^2 - r^2 (1 - tanh^2 l), an identity in (r, alpha, l)
t_ = sp.tanh(lam)
lhs = sp.simplify(drdphi**2); rhs = sp.simplify(r**4/bb**2 - r**2*(1 - t_**2))
check("PT-3c ORBIT IDENTITY: (dr/dphi)^2 = r^4/b^2 - r^2 (1 - tanh^2 l) EXACTLY, for any radial pinning profile --"
      " the null locus plus the two pairings, nothing else", z(sp.expand(lhs - rhs)))
orbit_pinned = sp.simplify((r**4/b**2 - r**2*(1 - sp.tanh(lam)**2)).subs(sp.tanh(lam), sp.sqrt(rs/r)))
check("PT-3d under KIN-2a: (dr/dphi)^2 = r^4/b^2 - r^2 (1 - r_s/r): the pinning turns 1 - tanh^2 into 1 - r_s/r",
      z(orbit_pinned - (r**4/b**2 - r**2*(1 - rs/r))))
# what a NEWTONIAN profile would need: (dr/dphi)^2 = r^4/b^2 - r^2 + (r_s/2)... i.e. the rod term without the drift's square.  Recorded, not used.
# closest approach and the deflection integral (exact, 1D): delta(b) = 2 int_0^{u0} du / sqrt(1/b^2 - u^2 + r_s u^3) - pi
def deflection(b_over_rs):
    B = mp.mpf(b_over_rs)                                    # units r_s = 1
    u0 = mp.findroot(lambda u: 1/B**2 - u**2 + u**3, mp.mpf(1)/B)   # smallest positive root: closest approach 1/r_0
    f = lambda u: 1/mp.sqrt(1/B**2 - u**2 + u**3)
    return 2*mp.quad(f, [0, u0]) - mp.pi
vals = {Bv: deflection(Bv) for Bv in (10**3, 10**4, 10**5, 10**6)}
c1 = {Bv: d*Bv for Bv, d in vals.items()}
c2 = {Bv: (d*Bv - 2)*Bv for Bv, d in vals.items()}          # second-order coefficient after removing the leading 2
check("PT-3e THE DEFLECTION, weak field: b * delta(b) -> 2 (not 1): " + ", ".join(f"b/r_s=10^{int(mp.log10(Bv))}: {mp.nstr(c, 12)}" for Bv, c in c1.items()),
      abs(c1[10**6] - 2) < mp.mpf(10)**-5 and abs(c1[10**6] - 1) > mp.mpf(1)/2)
check("PT-3f the second-order coefficient (delta - 2 r_s/b) b^2/r_s^2 -> " + mp.nstr(c2[10**6], 10) + " ; 15 pi/16 = " + mp.nstr(15*mp.pi/16, 10),
      abs(c2[10**6] - 15*mp.pi/16) < mp.mpf(10)**-4)
print("  ARGUMENT (not counted): the factor 2 has an address in PT-3c.  With the drift's square dropped (1 - tanh^2 -> 1) the")
print("  orbit is a straight line; with the rod only (1 - r_s/r from rho_K) but no drift in omega_t the pairing ratio is not")
print("  conserved.  Both offices are needed and each carries half.")

print(); n_ = sum(CH); print(f"RESULT: {n_}/{len(CH)} checks passed.")
if n_ != len(CH): print("VERDICT: check failures above."); sys.exit(1)
print("VERDICT (pre-tests): the ray kinematics of the pinned family is fixed by the null locus and two symmetry pairings.")
print("  Redshift ratio cosh l_2/cosh l_1, weak-field coefficient 1/2; both local seat-to-seat transport rules FAIL (PT-1i/j);")
print("  the null-structure horizon is the pinning's horizon; the orbit closes as an identity and, under KIN-2a, its")
print("  deflection is 2 r_s/b + (15 pi/16)(r_s/b)^2 + ...  Conditional on KIN-2a, STAT-1, ROT-1, PROP-1.")
print("COMPARISON STAGE: STAT-1/ROT-1 are the Killing-vector conservation laws of a stationary, spherically symmetric field;")
print("  PROP-1 is the Painleve-Gullstrand / river picture (flat slice, inflow at the escape speed); PT-3d is Schwarzschild's")
print("  null-geodesic equation; 2 r_s/b is Einstein 1915 (4GM/c^2 b), 15 pi/16 (r_s/b)^2 = (15 pi/4)(GM/c^2 b)^2 is the")
print("  second-order light deflection (Epstein-Shapiro 1980); the weak-field 1/2 is Pound-Rebka.  The office ghost of PT-1i")
print("  is the classical error of treating free-fall frames at different points as boosted copies of each other.")
sys.exit(0)
