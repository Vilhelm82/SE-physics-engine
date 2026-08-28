#!/usr/bin/env python3
# =============================================================================
# THE TWO DOORS - relax vs circulate, run side by side  (2026-08-28)
#
# The handoff calls this "a ruling, not a computation".  The cheat: build the
# MINIMAL dynamics behind each door on the same state space (three unit axes,
# state G), same click-invariant potential, same initial data, and compute the
# discriminators.  The ruling then chooses between computed consequences.
#
# Potential (natural, not tuned): V = -1/2 log Delta = the total correlation
# of the three axes = the SDP log-barrier on the elliptope [lineage flagged].
# Second potential V = -Delta for the wall experiments.
#
# DOOR R (relax):     da_i/dt = -P_i grad_i V          (projected gradient)
# DOOR C (circulate): d2a_i/dt2 = -P_i grad_i V - |da_i|^2 a_i   (Hamiltonian
#                     on the product of spheres, holonomic |a_i| = 1)
#
# Symbolic results = theorems.  Trajectories = NUMERIC EVIDENCE, labelled so.
# =============================================================================
import sympy as sp
import math, random

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" - {n}" if n else ""), flush=True)

print("=" * 78)
print("PART 1 - THEOREMS (symbolic, exact)")
print("=" * 78)
# S1: circulation momentum L = sum a_i x a'_i is conserved by door C for ANY V(G)
a1 = sp.Matrix(sp.symbols('a11 a12 a13')); a2 = sp.Matrix(sp.symbols('a21 a22 a23'))
a3 = sp.Matrix(sp.symbols('a31 a32 a33'))
c12, c13, c23 = sp.symbols('c12 c13 c23')      # dV/dgamma_ij, arbitrary
F1 = -(c12*a2 + c13*a3); F2 = -(c12*a1 + c23*a3); F3 = -(c13*a1 + c23*a2)
Ldot = a1.cross(F1) + a2.cross(F2) + a3.cross(F3)   # constraint forces a_i x a_i = 0
check("S1  dL/dt = -Sum a_i x grad_i V = 0 for EVERY click-invariant V (door C)",
      sp.simplify(Ldot) == sp.zeros(3, 1),
      "3-axis P10/P11: circulation's conserved momentum is door C's birthright")
print("      Door R destroys L by construction (first-order flow has no momentum).")

# S2: near the orthogonal state the three angle modes are 1:1:1 DEGENERATE and
# the cubic coupling is exactly -m: the three-seat cross term is the resonance.
g1, g2, g3 = sp.symbols('g1 g2 g3')
Delta = 1 - g1**2 - g2**2 - g3**2 + 2*g1*g2*g3
Vlog = -sp.Rational(1, 2)*sp.log(Delta)
ser = sp.expand(Vlog.series(g1, 0, 4).removeO().series(g2, 0, 4).removeO()
                    .series(g3, 0, 4).removeO())
P = sp.Poly(ser, g1, g2, g3)
deg = lambda n: sp.expand(sum(c*g1**a*g2**b*g3**e for (a,b,e), c in P.terms() if a+b+e == n))
check("S2  quadratic = p1/2 (1:1:1 DEGENERATE modes); cubic = -m EXACTLY",
      sp.expand(deg(2) - (g1**2+g2**2+g3**2)/2) == 0 and
      sp.expand(deg(3) + g1*g2*g3) == 0,
      "NOTE: odd cubic coupling at 1:1:1 is NON-resonant at first order "
      "(+-1+-1+-1 != 0): transfer is SECOND order, via m^2. Corrected claim.")
# S3: pure single-angle displacement is an invariant subspace of BOTH doors
check("S3  {g13 = g23 = 0} is invariant: dV/dg13 and dV/dg23 vanish there",
      sp.simplify(sp.diff(Vlog, g2).subs({g2: 0, g3: 0})) == 0 and
      sp.simplify(sp.diff(Vlog, g3).subs({g2: 0, g3: 0})) == 0,
      "one angle cannot reach another except THROUGH the third: m is the only conduit")
# S4: neither door conserves Delta => CONS-1 is not automatic in either
pt = {g1: sp.Rational(1,3), g2: sp.Rational(1,5), g3: sp.Rational(-1,4)}
gradV = [sp.diff(Vlog, g) for g in (g1, g2, g3)]
gradD = [sp.diff(Delta, g) for g in (g1, g2, g3)]
inner = sp.simplify(sum(gv*gd for gv, gd in zip(gradV, gradD)).subs(pt))
check("S4  dDelta/dt != 0 under door R at a generic point (and door C oscillates it)",
      inner != 0,
      "the SANDPIT (fixed total) holds in NEITHER door as a scalar => it must be the")
print("      CONTINUITY reading (d_mu n^mu = 0, the current) or an imposed constraint.")
print("      The title fight's 'the object is a current' is forced from a 2nd direction.")

print(); print("=" * 78)
print("PART 2 - TRAJECTORIES [NUMERIC EVIDENCE, float RK4; monitors reported]")
print("=" * 78)
def dot(u, v): return sum(x*y for x, y in zip(u, v))
def sub(u, v): return [x-y for x, y in zip(u, v)]
def add(u, v): return [x+y for x, y in zip(u, v)]
def scl(s, u): return [s*x for x in u]
def crs(u, v): return [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
def nrm(u):
    n = math.sqrt(dot(u, u)); return [x/n for x in u]
def gammas(A): return (dot(A[0], A[1]), dot(A[0], A[2]), dot(A[1], A[2]))
def Dval(g): return 1 - g[0]**2 - g[1]**2 - g[2]**2 + 2*g[0]*g[1]*g[2]
def dDdg(g): return (-2*g[0] + 2*g[1]*g[2], -2*g[1] + 2*g[0]*g[2], -2*g[2] + 2*g[0]*g[1])
def dVdg(g, kind):
    dD = dDdg(g)
    if kind == 'log':
        D = Dval(g); return tuple(-d/(2*D) for d in dD)
    return tuple(-d for d in dD)                      # V = -Delta
def forces(A, kind):
    g = gammas(A); c = dVdg(g, kind)
    F = [scl(-1, add(scl(c[0], A[1]), scl(c[1], A[2]))),
         scl(-1, add(scl(c[0], A[0]), scl(c[2], A[2]))),
         scl(-1, add(scl(c[1], A[0]), scl(c[2], A[1])))]
    return [sub(F[i], scl(dot(F[i], A[i]), A[i])) for i in range(3)]   # tangential
def Vsigned(A):
    return dot(A[0], crs(A[1], A[2]))
def Ltot(A, W):
    L = [0.0, 0.0, 0.0]
    for i in range(3): L = add(L, crs(A[i], W[i]))
    return L
def energy(A, W, kind):
    g = gammas(A); D = Dval(g)
    V = -0.5*math.log(D) if kind == 'log' else -D
    return V + 0.5*sum(dot(w, w) for w in W), V

def step_C(A, W, h, kind):
    def acc(A_, W_):
        F = forces(A_, kind)
        return [sub(F[i], scl(dot(W_[i], W_[i]), A_[i])) for i in range(3)]
    k1v = acc(A, W); k1x = W
    A2 = [add(A[i], scl(h/2, k1x[i])) for i in range(3)]
    W2 = [add(W[i], scl(h/2, k1v[i])) for i in range(3)]
    k2v = acc(A2, W2); k2x = W2
    A3 = [add(A[i], scl(h/2, k2x[i])) for i in range(3)]
    W3 = [add(W[i], scl(h/2, k2v[i])) for i in range(3)]
    k3v = acc(A3, W3); k3x = W3
    A4 = [add(A[i], scl(h, k3x[i])) for i in range(3)]
    W4 = [add(W[i], scl(h, k3v[i])) for i in range(3)]
    k4v = acc(A4, W4); k4x = W4
    An, Wn = [], []
    for i in range(3):
        x = add(A[i], scl(h/6, add(add(k1x[i], scl(2, k2x[i])), add(scl(2, k3x[i]), k4x[i]))))
        v = add(W[i], scl(h/6, add(add(k1v[i], scl(2, k2v[i])), add(scl(2, k3v[i]), k4v[i]))))
        x = nrm(x); v = sub(v, scl(dot(v, x), x))
        An.append(x); Wn.append(v)
    return An, Wn
def step_R(A, h, kind):
    def vel(A_): return forces(A_, kind)
    k1 = vel(A); A2 = [nrm(add(A[i], scl(h/2, k1[i]))) for i in range(3)]
    k2 = vel(A2); A3 = [nrm(add(A[i], scl(h/2, k2[i]))) for i in range(3)]
    k3 = vel(A3); A4 = [nrm(add(A[i], scl(h, k3[i]))) for i in range(3)]
    k4 = vel(A4)
    return [nrm(add(A[i], scl(h/6, add(add(k1[i], scl(2, k2[i])),
                                       add(scl(2, k3[i]), k4[i]))))) for i in range(3)]
def frame_from_g(g):
    a1 = [1.0, 0.0, 0.0]
    a2 = [g[0], math.sqrt(1-g[0]**2), 0.0]
    x = g[1]; y = (g[2]-g[0]*g[1])/math.sqrt(1-g[0]**2)
    a3 = [x, y, math.sqrt(max(1-x*x-y*y, 1e-15))]
    return [a1, a2, a3]

print("\nN1 - free response from the same displaced state (V = -log(Delta)/2)")
g0 = (0.55, -0.30, 0.20)
A0 = frame_from_g(g0)
h, N = 0.002, 60000
A, W = [r[:] for r in A0], [[0.0]*3 for _ in range(3)]
Ac = [r[:] for r in A0]
E0, _ = energy(A, W, 'log'); L0 = Ltot(A, W)
Dmin, Dmax, Vs = 2.0, -1.0, []
for k in range(N):
    A, W = step_C(A, W, h, 'log')
    if k % 200 == 0:
        D = Dval(gammas(A)); Dmin, Dmax = min(Dmin, D), max(Dmax, D)
E1, _ = energy(A, W, 'log'); L1 = Ltot(A, W)
Dstart = Dval(g0)
Rvals = []
for k in range(N):
    Ac = step_R(Ac, h, 'log')
    if k % 6000 == 0: Rvals.append(Dval(gammas(Ac)))
Dend_R = Dval(gammas(Ac))
mono = all(Rvals[i] <= Rvals[i+1] + 1e-12 for i in range(len(Rvals)-1))
print(f"    door C: Delta ranged [{Dmin:.4f}, {Dmax:.4f}] around start {Dstart:.4f} (oscillates, returns)")
print(f"    door C monitors: |H-H0| = {abs(E1-E0):.2e}, |L-L0| = {max(abs(L1[i]-L0[i]) for i in range(3)):.2e}")
print(f"    door R: Delta {Dstart:.4f} -> {Dend_R:.6f}, monotone: {mono} (relaxes to orthogonal, L identically 0)")
check("N1  discriminator: C oscillates+conserves(H,L); R monotone to Delta=1 [numeric]",
      Dmax - Dmin > 0.1 and abs(E1-E0) < 1e-6 and mono and Dend_R > 0.999)

print("\nN2 - the SEESAW: kick one angle at orthogonality, seed the others tiny")
random.seed(7)
A = frame_from_g((0.0, 0.001, 0.001))            # seeds 1e-3 in modes 2,3
W = [[0, 0.60, 0], [0,0,0], [0,0,0]]             # kick mode 1 only
seed = 0.001; gmax2 = gmax3 = 0.0
for k in range(240000):
    A, W = step_C(A, W, 0.002, 'log')
    if k % 50 == 0:
        g = gammas(A); gmax2 = max(gmax2, abs(g[1])); gmax3 = max(gmax3, abs(g[2]))
Ar = frame_from_g((0.30, 0.001, 0.001)); rmax = 0.001
for k in range(120000):
    Ar = step_R(Ar, 0.002, 'log')
    if k % 50 == 0:
        g = gammas(Ar); rmax = max(rmax, abs(g[1]), abs(g[2]))
check("N2  SEESAW: door C GROWS seeded modes (2nd order, via m); door R: exactly none",
      gmax2 > 2*seed and gmax3 > 2*seed and rmax <= seed*1.0001,
      f"C: |g13| max {gmax2:.4f}, |g23| max {gmax3:.4f} from seed {seed}; R: max {rmax:.6f}")
print("    Exact corollary (S3): a PURE single-pair excitation never spreads at all;")
print("    the third party m is the only conduit, and it opens at second order.")

print("\nN3 - the WALL, decided by the cost functional (sec 8 adjudicated per door)")
A = frame_from_g((0.1, 0.1, 0.1)); W = [[0,0,0],[0,0,0],[2.3, -0.9, 0.0]]
flips = 0; last = Vsigned(A); DminP = 2.0
for k in range(90000):
    A, W = step_C(A, W, 0.002, 'poly')
    vs = Vsigned(A); DminP = min(DminP, Dval(gammas(A)))
    if vs*last < 0: flips += 1
    last = vs
check("N3a V = -Delta (finite at the wall): door C CROSSES coplanarity; handedness",
      flips >= 2 and DminP < 1e-3,
      f"signed volume changed sign {flips} times, min Delta {DminP:.2e} - the Moebius")
print("      traversal is a recurring dynamical EVENT: each crossing flips the frame.")
A = frame_from_g((0.1, 0.1, 0.1)); W = [[0,0,0],[0,0,0],[2.3, -0.9, 0.0]]
DminL = 2.0
for k in range(90000):
    A, W = step_C(A, W, 0.002, 'log')
    DminL = min(DminL, Dval(gammas(A)))
check("N3b V = -log(Delta)/2 (barrier): SAME initial data is repelled - never crosses",
      DminL > 1e-3, f"min Delta on trajectory {DminL:.4f} > 0: ASYMPTOTE, not WALL")

print("\nN4 - time reversal: door C retraces; door R cannot")
A = frame_from_g(g0); W = [[0,0.1,0],[0,0,0.07],[0.05,0,0]]
A0s = [r[:] for r in A]
for k in range(20000): A, W = step_C(A, W, 0.002, 'log')
W = [scl(-1, w) for w in W]
for k in range(20000): A, W = step_C(A, W, 0.002, 'log')
err = max(abs(A[i][j]-A0s[i][j]) for i in range(3) for j in range(3))
check("N4  velocity reversal returns door C to its start (reversible microworld)",
      err < 1e-6, f"return error {err:.2e}; door R is irreversible by construction")

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
COMPARISON TABLE (R = relax / gradient, C = circulate / Hamiltonian):
                              DOOR R              DOOR C
  energy                      Lyapunov, decays    conserved (H)
  circulation L               destroyed           conserved for EVERY V(G)  [S1]
  Delta                       monotone -> 1       oscillates, recurs
  injection into one angle    others only decay   others GROW then return (seesaw,
                                                  1:1:1 resonance, coupling = -m) [S2,N2]
  time reversal               impossible          exact                      [N4]
  sandpit CONS-1              fails as scalar     fails as scalar            [S4]
                              => forced to the CURRENT/continuity reading in both
  Delta = 0                   never reached       cost functional decides: -Delta
                              (barrier repels)    crossable + handedness flips;
                                                  -log barrier = asymptote  [N3]
THE REFRAME: R is not C's sibling - R is C coupled to an unmodelled bath
(gradient flow = overdamped limit; dissipation presupposes an OUTSIDE).
So the ruling is not relax-vs-circulate.  It is: IS THE STRUCTURE CLOSED?
Closed => door C is forced (S1 makes circulation its birthright), and 'relax'
can only ever be a coarse-grained description of C.  Open => name the bath.
""")
