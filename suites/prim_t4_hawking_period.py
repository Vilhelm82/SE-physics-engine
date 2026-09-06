#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T4 -- Hawking's coefficient is the period of the Wick face times the
#                  scale of the logarithmic twist.  Cl(3) BANNED.   (2026-09-04)
#
# INPUTS: T2's SO(2,1) (two hyperbolic planes through the timelike line c, one
#   compact plane); P6 pivot = imaginary-angle rotation; the pinning
#   tanh(lambda) = sqrt(r_s/r) as the model's own radial profile (THM-K derived it
#   from the thermal declarations; THM-O re-derived its exponent from flat rods;
#   here it is the INPUT whose consequences we read).  Will's ruling: c is timelike.
# BANNED: Cl(3), Pauli, spinor as a given, Bogoliubov, field theory, any Euclidean
#   black-hole metric, any temperature as a given.  SO(2,1)'s double cover SL(2,R)
#   is REACHED below (T4b': the vector rep is the symmetric square of the 2-dim rep).
#
# HELD-OUT PREDICTIONS (written before running):
#   T4a  vector rep: exp(i theta K) has minimal period 2 pi; at theta = pi it is -1
#        on the boost plane (the deck at i^2), +1 off it.
#   T4b  2-dim rep: exp(i theta K_s) has minimal period 4 pi; at 2 pi it is -I.
#   T4b' the 3-dim boost generator IS the symmetric square of the 2-dim one:
#        two one-sided readings compose to one two-sided reading (SECT-1's census).
#   T4c  lambda(r) = ln(sqrt r + sqrt r_s) - (1/2) ln(r - r_s) exactly: logarithmic
#        twist, slope -1/2 in r.
#   T4d  in proper distance rho from the cone, lambda = ln(4 r_s / rho) + o(1):
#        slope -1 (universal), scale 4 r_s (where the physics lives).
#   T4e  kappa two ways: N/rho -> 1/(2 r_s) from the twist's scale, and
#        alpha N -> 1/(2 r_s) from the seat's redshifted acceleration (THM-K's
#        K-1 identity re-derived from lambda alone).  Same number, no fitting.
#   T4f  T = hbar (alpha N) / (2 pi c k_B) with 2 pi from T4a  ->  hbar c / (4 pi k_B r_s).
#        KILL: if the seat reading presented rods (vectors) had period 4 pi, T halves
#        and dE = T dS with S = A/4 breaks by 2 (THM-K's first-law check).
#   T4g  Euclidean regularity: the (tau, rho) plane near the cone is flat polar
#        coordinates iff tau has period 2 pi / kappa.  Same 2 pi, geometric route.
# OUTCOME FORK: (a) all hold -> Hawking's coefficient is DERIVED as (Wick period of
#   the two-sided seat) x (scale of the log twist), the tensoriality theorem is no
#   longer load-bearing (the vector rep has no -1 to hide).  (b) T4b' fails -> the
#   cover is not reached and hbar's 4 pi is an import.  (c) T4e's two kappas differ
#   -> the twist and the acceleration are different objects; record the ratio.
# =============================================================================
import sympy as sp, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand(q)),
              lambda q: sp.simplify(sp.expand(q.rewrite(sp.exp))), sp.cancel):
        try:
            if f(e) == 0: return True
        except Exception: pass
    return False
def zM(M): return all(z(e) for e in M)
I = sp.I
th = sp.Symbol('theta', real=True)

print("=== T4a/b: the Wick face's period on the two-sided and one-sided readings ===")
# SO(2,1), metric diag(1,1,-1), x3 timelike (= c, Will's ruling).  Boost in the (1,3) plane.
Q21 = sp.diag(1, 1, -1)
K = sp.Matrix([[0,0,1],[0,0,0],[1,0,0]])
check("T4a0 K is an so(2,1) generator (K^T Q + Q K = 0) and K^3 = K", zM(K.T*Q21 + Q21*K) and zM(K**3 - K))
# exp(i theta K) from K^3 = K: E = I + (cos th - 1) K^2 + i sin th K.  Verify it IS the exponential.
E = sp.eye(3) + (sp.cos(th) - 1)*K**2 + I*sp.sin(th)*K
check("T4a1 E(theta) solves E' = i K E with E(0) = I (it is exp(i theta K))",
      zM(E.diff(th) - I*K*E) and zM(E.subs(th, 0) - sp.eye(3)))
check("T4a2 E(2 pi) = I: the Wick face closes at 2 pi on the vector reading", zM(E.subs(th, 2*sp.pi) - sp.eye(3)))
check("T4a3 E(pi) = diag(-1, 1, -1): at i^2 the boost plane is inverted, the compact direction untouched (the deck)",
      zM(E.subs(th, sp.pi) - sp.diag(-1, 1, -1)))
# minimal period: E(theta) = I  <=>  cos th = 1 and sin th = 0.  (K^2 and K are independent, nonzero.)
# E - I = (cos th - 1) K^2 + i sin th K, and K^2, K are linearly independent nonzero matrices,
# so E = I iff cos th = 1 and sin th = 0.  Check independence and that the only solutions in
# [0, 4 pi) are 0 and 2 pi (so 2 pi is the minimal positive period, and 4 pi is not new).
indep = sp.Matrix([list(K**2), list(K)]).rank() == 2
sols = sp.solveset(sp.Eq(sp.cos(th), 1), th, sp.Interval(0, 4*sp.pi - sp.Rational(1,10**6)))
check("T4a4 E(theta) = I only at theta in 2 pi Z (minimal period exactly 2 pi)",
      indep and sols == sp.FiniteSet(0, 2*sp.pi), str(sols))
# numeric second path
import mpmath as mp; mp.mp.dps = 30
Kn = mp.matrix([[0,0,1],[0,0,0],[1,0,0]])
En = lambda t: mp.expm(mp.mpc(0,1)*t*Kn)
check("T4a5 mpmath expm agrees: |E(2pi) - I| < 1e-25, |E(pi) - diag(-1,1,-1)| < 1e-25",
      mp.norm(En(2*mp.pi) - mp.eye(3)) < mp.mpf('1e-25') and mp.norm(En(mp.pi) - mp.diag([-1,1,-1])) < mp.mpf('1e-25'))

# The 2-dim real rep: sl(2,R).  A boost generator there is real, symmetric, traceless.
Ks = sp.Rational(1,2)*sp.Matrix([[1,0],[0,-1]])
check("T4b0 K_s is in sl(2,R): real, traceless", Ks.trace() == 0 and all(sp.im(e) == 0 for e in Ks))
Es = (I*th*Ks).exp()
check("T4b1 exp(i theta K_s) = diag(e^{i th/2}, e^{-i th/2})",
      zM(Es - sp.diag(sp.exp(I*th/2), sp.exp(-I*th/2))))
check("T4b2 at theta = 2 pi the one-sided reading returns -I (the sign the two-sided reading cannot see)",
      zM(Es.subs(th, 2*sp.pi) + sp.eye(2)))
check("T4b3 at theta = 4 pi it returns +I: minimal period 4 pi on the cover",
      zM(Es.subs(th, 4*sp.pi) - sp.eye(2)) and not zM(Es.subs(th, 2*sp.pi) - sp.eye(2)))

print("=== T4b': the vector reading is the symmetric square of the one-sided reading (cover REACHED) ===")
# Sym^2 of the 2-dim rep: basis {e1e1, e1e2 sym, e2e2}.  The induced action of a 2x2 A on
# symmetric 2x2 matrices S is S -> A S A^T.  Its generator for A = exp(t Ks) is S -> Ks S + S Ks^T.
# Map symmetric S = [[s11,s12],[s12,s22]] to the vector (s11, s22, s12) up to a normalisation,
# and read the 3x3 generator; then check it is conjugate to K (same eigenvalues, preserves a
# signature-(2,1) form) -- i.e. the two-sided reading of a boost is two one-sided readings.
s11, s12, s22 = sp.symbols('s11 s12 s22', real=True)
S = sp.Matrix([[s11, s12], [s12, s22]])
dS = Ks*S + S*Ks.T
vec = lambda M: sp.Matrix([M[0,0], M[1,1], M[0,1]])
Gsym = sp.Matrix([[sp.diff(vec(dS)[i], v) for v in (s11, s22, s12)] for i in range(3)])
check("T4b'1 the induced generator on Sym^2 is real with eigenvalues {1, -1, 0} = pairwise sums of {+1/2, -1/2}",
      sorted(sp.Matrix(Gsym).eigenvals().keys(), key=lambda x: float(x)) == [-1, 0, 1])
# it preserves the determinant form on symmetric matrices, s11 s22 - s12^2, which has signature (2,1)
# (or (1,2)): that is the (2,1) metric the vector rep carries.  Verify: d/dt det(A S A^T) = 0 <=> generator
# is in the orthogonal algebra of the form s11 s22 - s12^2.
Fdet = sp.Matrix([[0, sp.Rational(1,2), 0],[sp.Rational(1,2), 0, 0],[0, 0, -1]])   # s11 s22 - s12^2 as a quadratic form
check("T4b'2 the Sym^2 generator preserves det S = s11 s22 - s12^2 (a signature-(2,1) form): it is an so(2,1) boost",
      zM(Gsym.T*Fdet + Fdet*Gsym))
check("T4b'3 Sym^2 generator is conjugate to K: same characteristic polynomial",
      z(sp.expand(Gsym.charpoly(th).as_expr() - K.charpoly(th).as_expr())))
check("T4b'4 exp(i theta Gsym) has period 2 pi while exp(i theta Ks) has 4 pi: the two-sided reading of the cover's -1 is +1",
      zM((I*2*sp.pi*Gsym).exp() - sp.eye(3)))

print("=== T4c/d: the twist is logarithmic; slope universal, scale carries the physics ===")
r, rs, rho = sp.symbols('r r_s rho', positive=True)
lam_r = sp.atanh(sp.sqrt(rs/r))                       # the pinning, as the model's own profile
lam_closed = sp.log(sp.sqrt(r) + sp.sqrt(rs)) - sp.Rational(1,2)*sp.log(r - rs)
# same derivative everywhere and same value at one point => same function on r > r_s
d_match = z(sp.simplify(sp.diff(lam_r, r) - sp.diff(lam_closed, r)))
pt = z(sp.simplify((lam_r - lam_closed).subs({r: 4*rs}).rewrite(sp.log)))   # atanh(1/2) = ln 3 / 2 ; closed: ln(3 sqrt r_s) - ln(sqrt(3 r_s)) = ln sqrt 3
num_ok = all(abs(mp.atanh(mp.sqrt(mp.mpf(1)/x)) - (mp.log(mp.sqrt(x)+1) - mp.log(x-1)/2)) < mp.mpf('1e-28') for x in (mp.mpf('1.5'), mp.mpf(3), mp.mpf(50)))
check("T4c1 lambda(r) = ln(sqrt r + sqrt r_s) - (1/2) ln(r - r_s) exactly (derivative + one point + 30-digit spot checks)",
      d_match and pt and num_ok)
eps = sp.Symbol('epsilon', positive=True)             # r = r_s + eps
lam_near = lam_closed.subs(r, rs + eps)
lead = sp.series(lam_near, eps, 0, 1).removeO()
check("T4c2 near the cone: lambda = -(1/2) ln(eps) + ln(2 sqrt r_s) + O(eps): slope -1/2 in r, scale 2 sqrt r_s",
      z(lead - (-sp.Rational(1,2)*sp.log(eps) + sp.log(2*sp.sqrt(rs)))))
# proper radial distance from the cone: d rho = dr / N, N = sech(lambda) = sqrt(1 - r_s/r)
N = 1/sp.cosh(lam_r)
check("T4c3 N = sech(lambda) = sqrt(1 - r_s/r) (the pinning's lapse)", z(sp.simplify(N**2 - (1 - rs/r))))
rho_r = sp.integrate(1/sp.sqrt(1 - rs/r), (r, rs, r))
rho_closed = sp.sqrt(r*(r - rs)) + rs*sp.log((sp.sqrt(r) + sp.sqrt(r - rs))/sp.sqrt(rs))
check("T4c4 rho(r) = sqrt(r(r-r_s)) + r_s ln((sqrt r + sqrt(r-r_s))/sqrt r_s) (proper distance to the cone)",
      z(sp.simplify(sp.diff(rho_closed, r) - 1/sp.sqrt(1 - rs/r))) and z(rho_closed.subs(r, rs)))
lim_ratio = sp.limit(rho_closed.subs(r, rs + eps)/sp.sqrt(eps), eps, 0, '+')
check("T4c5 near the cone rho / sqrt(eps) -> 2 sqrt(r_s): rho = 2 sqrt(r_s eps) + o(sqrt eps)", z(lim_ratio - 2*sp.sqrt(rs)), str(lim_ratio))
# invert: eps = rho^2/(4 r_s); substitute into the log expansion
lam_rho = lead.subs(eps, rho**2/(4*rs))
check("T4d1 lambda = ln(4 r_s / rho) + o(1): slope -1 in proper distance (universal), scale 4 r_s",
      z(sp.expand_log(lam_rho - sp.log(4*rs/rho), force=True)))
# N near the cone from the twist's scale: N = sech(lambda) ~ 2 e^{-lambda} = rho / (2 r_s)
N_from_twist = 2*sp.exp(-lam_rho)
check("T4d2 N ~ 2 e^{-lambda} = rho/(2 r_s): the lapse is linear in proper distance with slope 1/(2 r_s)",
      z(sp.simplify(N_from_twist - rho/(2*rs))))

print("=== T4e: kappa two ways ===")
c = sp.Symbol('c', positive=True)
kappa_twist = sp.limit((N.subs(r, rs + eps)) / (2*sp.sqrt(rs*eps)), eps, 0) * c**2     # N/rho at the cone
check("T4e1 kappa from the twist's scale: N/rho -> 1/(2 r_s)  (times c^2)", z(kappa_twist - c**2/(2*rs)))
# THM-K's route, re-derived from lambda alone: static proper acceleration alpha = (c^2/2)|d beta^2/dr| / N,
# beta = tanh(lambda); the redshifted acceleration alpha N.
beta2 = sp.tanh(lam_r)**2
alphaN = sp.Rational(1,2)*c**2*sp.Abs(sp.diff(beta2, r))
check("T4e2 alpha N = (c^2/2) r_s / r^2 for every r (the K-1 identity from lambda alone)",
      z(sp.simplify(alphaN - c**2*rs/(2*r**2))))
kappa_accel = alphaN.subs(r, rs)
check("T4e3 kappa from the redshifted acceleration: alpha N -> c^2/(2 r_s) at the cone", z(kappa_accel - c**2/(2*rs)))
check("T4e4 THE TWO KAPPAS ARE THE SAME NUMBER, no fitting", z(kappa_twist - kappa_accel))

print("=== T4f/g: the temperature, and the Euclidean route to the same 2 pi ===")
hbar, kB = sp.symbols('hbar k_B', positive=True)
period_vector = 2*sp.pi     # T4a2/T4a4
period_cover  = 4*sp.pi     # T4b3
T_from_vector = hbar*kappa_accel/(period_vector*c*kB)
T_from_cover  = hbar*kappa_accel/(period_cover*c*kB)
check("T4f1 T = hbar kappa / (2 pi c k_B) = hbar c / (4 pi k_B r_s)  [the two-sided seat's period]",
      z(T_from_vector - hbar*c/(4*sp.pi*kB*rs)))
check("T4f2 the cover's period would give exactly half (the live factor-2 error, now impossible on the vector reading)",
      z(T_from_cover - T_from_vector/2))
# first-law consistency (THM-K): with S = k_B A/(4 l_P^2), A = 4 pi r_s^2, r_s = 2GM/c^2, E = M c^2:
G, M = sp.symbols('G M', positive=True)
lP2 = hbar*G/c**3
S_area = kB*4*sp.pi*(2*G*M/c**2)**2/(4*lP2)
dE_dM = c**2
dS_dM = sp.diff(S_area, M)
T_firstlaw = dE_dM/dS_dM
check("T4f3 first law dE = T dS with S = A/4 returns T_from_vector exactly (and NOT T_from_cover)",
      z(sp.simplify(T_firstlaw - T_from_vector.subs(rs, 2*G*M/c**2))) and not z(sp.simplify(T_firstlaw - T_from_cover.subs(rs, 2*G*M/c**2))))
# Euclidean route: near the cone the (t, rho) plane has ds^2 = -N^2 dt^2 + d rho^2 with N = kappa rho / c^2 (T4d2).
# Wick t -> -i tau: ds^2 = (kappa rho / c^2)^2 d tau^2 + d rho^2.  Put u = (kappa/c^2) tau: ds^2 = rho^2 du^2 + d rho^2,
# flat polar coordinates, smooth at rho = 0 iff u has period 2 pi, i.e. tau has period 2 pi c^2 / kappa.
tau, u = sp.symbols('tau u', real=True)
kap = c**2/(2*rs)
metric_E = sp.Matrix([[(kap*rho/c**2)**2, 0],[0, 1]])           # (tau, rho)
polar    = sp.Matrix([[rho**2, 0],[0, 1]])                       # (u, rho)
J = sp.Matrix([[kap/c**2, 0],[0, 1]])                            # d(u,rho)/d(tau,rho)
check("T4g1 the Wick face of the near-cone plane is flat polar coordinates in u = (kappa/c^2) tau",
      zM(J.T*polar*J - metric_E))
period_tau = 2*sp.pi*c**2/kap
T_euclid = hbar/(kB*period_tau/c)     # circle circumference hbar/(k_B T) in time units: tau-period = hbar c/(k_B T) -> T = hbar c/(k_B period)
check("T4g2 regularity (no conical deficit) forces tau-period 2 pi c^2/kappa = 4 pi r_s, and hbar c/(k_B period) = T_from_vector",
      z(sp.simplify(period_tau - 4*sp.pi*rs)) and z(sp.simplify(T_euclid - T_from_vector)))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT (computed from the table):")
print("  T4a/b :", "vector period 2 pi, cover period 4 pi, deck at i^2 -- DERIVED" if all(CH[0:10]) else "period anomaly, read FAILs")
print("  T4b'  :", "the cover is REACHED: vector = Sym^2(2-dim), the -1 squares to +1" if all(CH[10:14]) else "cover NOT reached -> 4 pi is an import")
print("  T4c/d :", "twist logarithmic; slope universal (-1/2 in r, -1 in rho); scale 4 r_s" if all(CH[14:21]) else "twist anomaly")
print("  T4e   :", "kappa from the twist = kappa from the acceleration = c^2/(2 r_s)" if all(CH[21:25]) else "the two kappas DIFFER -> fork (c)")
print("  T4f/g :", "T_H = hbar c/(4 pi k_B r_s) from (2 pi) x (kappa); first law and Euclidean regularity agree" if all(CH[25:30]) else "assembly anomaly")
