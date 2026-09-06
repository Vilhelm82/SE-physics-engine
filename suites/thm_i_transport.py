#!/usr/bin/env python3
# =============================================================================
# THM-TARGET I, PART A -- the local transport rule: seat rotors versus the gradient of the pivot field
# Date: 2026-09-03.  sympy + mpmath only.  Exit 0 iff every check passes.
#
# QUESTION: is there a LOCAL rule X(x+dx) = Lambda X Lambda^dagger, with Lambda built at x, that reproduces the
#   two conserved pairings of thm_i_pre (STAT-1 omega_t, ROT-1 L) and hence the exact orbit -- or is propagation
#   intrinsically global?  Three candidates, all infinitesimal Lorentz maps applied to the ray's null paravector:
#   RULE N  naive product of the seat rotors, A(lambda(x+dx)) A(lambda(x))^-1, both boosts against one background
#           (the D1 non-collinear composition, with its Wigner rotation)
#   RULE F  the frame-transported relative rotor of thm_h2_d1 D-14, R(d theta) A(d lambda, K)
#   RULE G  the boost by the GRADIENT of the pinned velocity field along the ray's co-moving displacement,
#           d lambda = (e . grad) v,  e = n dt  (the ray's own step in the local frame).  Not a product of seat rotors.
# CARRIED: BARE-1, RULE-1 pairing, E-4, KIN-2a, STAT-1, ROT-1, PROP-1 (thm_i_pre.py, 19/19).
# NOT USED before the comparison block: spin connections, Christoffel symbols, geodesic equations, tetrads.
# =============================================================================
import sys
import sympy as sp
import mpmath as mp
mp.mp.dps = 30
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
def z(e):
    e = sp.sympify(e)
    for rt in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), lambda q: sp.simplify(sp.expand_trig(q))):
        try:
            if rt(e) == 0: return True
        except Exception: pass
    return False

s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -sp.I], [sp.I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]]); Id = sp.eye(2)
SIG = [s1, s2, s3]
def para(t, v): return t*Id + sum((v[i]*SIG[i] for i in range(3)), sp.zeros(2))
def scal(M): return sp.trace(M)/2
def vect(M): return sp.Matrix([sp.trace(M*S)/2 for S in SIG])
def dag(M): return M.T.applyfunc(lambda e: e.subs(sp.I, -sp.I))   # Hermitian conjugate; every symbol and function here is real
def act(Lam, X): return Lam*X*dag(Lam)                      # SL(2,C) action on paravectors: boosts A X A, rotations R X R^dagger

# ---------------------------------------------------------------------------
# SYMBOLIC: the plane z = 0, a radial pinned field v(r) with v(r) = tanh(lambda(r)), inflow  vvec = -v r-hat
# ---------------------------------------------------------------------------
x, y, eps, w = sp.symbols('x y epsilon omega', real=True)
psi = sp.Symbol('psi', real=True)                             # ray direction angle in the frame's flat triad
r = sp.sqrt(x**2 + y**2)
vf = sp.Function('v', real=True)                              # general radial profile
vvec = -vf(r)*sp.Matrix([x, y, 0])/r                          # inflow at speed v(r)
n = sp.Matrix([sp.cos(psi), sp.sin(psi), 0])
X = w*para(1, n)                                              # the ray, null
J = sp.Matrix(3, 3, lambda i, j: sp.diff(vvec[i], [x, y, sp.Symbol('zz')][j]) if j < 2 else 0)   # J_ij = d_j v_i
check("A-0  the pinned field is irrotational: its gradient is symmetric (d_j v_i = d_i v_j)", z(J[0, 1] - J[1, 0]))

def conserved(Xp, xp, yp, field):
    """(omega_t, L) of a ray paravector at position (xp, yp): STAT-1 and ROT-1 pairings."""
    om = scal(Xp); nv = vect(Xp)/om
    vv = field.subs({x: xp, y: yp}, simultaneous=True)
    return om*(1 + (vv.T*nv)[0, 0]), om*(xp*nv[1] - yp*nv[0])
def first_order_drift(Lam_eps, field, simp=True):
    """d/d eps at 0 of (omega_t, L) after the step x -> x + dx with the ray mapped by Lam_eps (a function of eps)."""
    dxv = (n + field)*eps                                     # PROP-1: the ray advances at c along n in the frame; the frame drifts
    xp, yp = x + dxv[0], y + dxv[1]
    Xp = act(Lam_eps, X)
    Om1, L1 = conserved(Xp, xp, yp, field)
    d1, d2 = sp.diff(Om1, eps).subs(eps, 0), sp.diff(L1, eps).subs(eps, 0)
    return (sp.simplify(d1), sp.simplify(d2)) if simp else (d1, d2)

# RULE G: boost by the gradient along the co-moving step e = n eps;  passive map = A(-d lambda)
dlamG = (J*n)*eps                                             # (e . grad) v
LamG = Id - para(0, dlamG)/2                                  # first order; exact below in the numerics
dOmG, dLG = first_order_drift(LamG, vvec)
check("A-1  RULE G conserves BOTH pairings to first order, for EVERY radial profile v(r) and every ray direction: d omega_t = d L = 0",
      z(dOmG) and z(dLG), "identically in x, y, psi, omega, and the function v")

# the concrete KIN-2a profile for the two seat-rotor rules (r_s = 1): v = r^(-1/2)
vk = lambda R_: 1/sp.sqrt(R_)
vvec_k = -vk(r)*sp.Matrix([x, y, 0])/r
dxk = (n + vvec_k)*eps; xpk, ypk = x + dxk[0], y + dxk[1]
pt = {x: sp.Rational(3, 2), y: sp.Rational(4, 5), psi: sp.Rational(1, 3), w: 1}
def num(e): return complex(sp.N(e.subs(pt), 20))   # numeric at a rational point, no simplification
def fmt(c): return f'{c.real:.6g}' + (f' + {c.imag:.2g}i' if abs(c.imag) > 1e-12 else '')
# RULE N: the naive product of seat rotors.  rain frame at x = background boosted by rapidity vector lam_vec(x) = atanh(v) (-r-hat)
lamvec = lambda X_, Y_: -sp.atanh(vk(sp.sqrt(X_**2 + Y_**2)))*sp.Matrix([X_, Y_, 0])/sp.sqrt(X_**2 + Y_**2)
def boostvec(lv):
    l = sp.sqrt((lv.T*lv)[0, 0]); return sp.cosh(l/2)*Id + sp.sinh(l/2)*para(0, lv)/l
LamN = boostvec(lamvec(xpk, ypk))*boostvec(-lamvec(x, y))     # A(x + dx) A(x)^-1, exact in eps; differentiated below
dOmN, dLN = first_order_drift(LamN, vvec_k, simp=False)
check("A-2  RULE N (naive product of seat rotors) FAILS: d omega_t != 0 at a generic point under KIN-2a",
      abs(num(dOmN)) > 1e-6, f"d omega_t/dt = {fmt(num(dOmN))}, d L/dt = {fmt(num(dLN))}")

# RULE F: D-14's frame-transported relative rotor: rotation of K by the angle swept, then a collinear boost by the scalar rapidity change
Khat = -sp.Matrix([x, y, 0])/r
dlam_scalar = sp.atanh(vk(sp.sqrt(xpk**2 + ypk**2))) - sp.atanh(vk(r))
dtheta = sp.atan2(ypk, xpk) - sp.atan2(y, x)
Rz = sp.cos(dtheta/2)*Id - sp.I*sp.sin(dtheta/2)*s3
LamF = Rz*(sp.cosh(dlam_scalar/2)*Id + sp.sinh(dlam_scalar/2)*para(0, Khat))
dOmF, dLF = first_order_drift(LamF, vvec_k, simp=False)
check("A-3  RULE F (D-14 frame-transported relative rotor) FAILS: d omega_t != 0 at the same point",
      abs(num(dOmF)) > 1e-6, f"d omega_t/dt = {fmt(num(dOmF))}, d L/dt = {fmt(num(dLF))}")
dOmGk, dLGk = first_order_drift(Id - para(0, (sp.Matrix(3, 3, lambda i, j: sp.diff(vvec_k[i], [x, y, sp.Symbol('zz')][j]) if j < 2 else 0)*n)*eps)/2, vvec_k)
check("A-4  ...while RULE G at that point gives exactly 0 (the identity A-1, instantiated on KIN-2a)", z(dOmGk) and z(dLGk))

# ---------------------------------------------------------------------------
# NUMERICS: integrate each rule as an ODE for a ray at impact parameter b = 100 r_s, from r0 = 2000 r_s in and back out;
# compare the swept polar angle with the exact quadrature of the orbit identity between the same radii.
# ---------------------------------------------------------------------------
print(); print("=" * 100); print("NUMERICS -- each rule as an ODE, against the exact orbit between the same radii"); print("=" * 100)
import math
def generator_rhs(Lam_eps, field):
    """d ln omega/dt and d psi/dt from the first-order action of Lam_eps on the ray (symbolic -> lambdified to floats)."""
    Xp = act(Lam_eps, X); dX = sp.diff(Xp, eps).subs(eps, 0)
    dom = scal(dX); dn = (vect(dX) - n*dom)/w
    dpsi = -sp.sin(psi)*dn[0] + sp.cos(psi)*dn[1]
    import numpy as np
    f1 = sp.lambdify((x, y, psi), (dom/w).subs(w, 1), 'numpy', cse=True); f2 = sp.lambdify((x, y, psi), dpsi.subs(w, 1), 'numpy', cse=True)
    return (lambda X_, Y_, P_: float(np.real(f1(X_, Y_, P_)))), (lambda X_, Y_, P_: float(np.real(f2(X_, Y_, P_))))
Jk = sp.Matrix(3, 3, lambda i, j: sp.diff(vvec_k[i], [x, y, sp.Symbol('zz')][j]) if j < 2 else 0)
RULES = {'G': Id - para(0, (Jk*n)*eps)/2, 'N': LamN, 'F': LamF}
vxf = sp.lambdify((x, y), vvec_k[0], 'math'); vyf = sp.lambdify((x, y), vvec_k[1], 'math')
def integrate(rule, b=100.0, r0=1000.0, h=0.1):
    fom, fpsi = generator_rhs(RULES[rule], vvec_k)
    def rhs(s):
        X_, Y_, P_, LO_ = s
        return (math.cos(P_) + vxf(X_, Y_), math.sin(P_) + vyf(X_, Y_), fpsi(X_, Y_, P_), fom(X_, Y_, P_))
    s = [-math.sqrt(r0**2 - b**2), b, 0.0, 0.0]
    def cons(s):
        X_, Y_, P_, LO_ = s; om = math.exp(LO_); nx, ny = math.cos(P_), math.sin(P_)
        return om*(1 + vxf(X_, Y_)*nx + vyf(X_, Y_)*ny), om*(X_*ny - Y_*nx)
    c0 = cons(s); phi0 = math.atan2(s[1], s[0]); passed = False
    while True:
        k1 = rhs(s); s2 = [s[i] + h/2*k1[i] for i in range(4)]; k2 = rhs(s2)
        s3 = [s[i] + h/2*k2[i] for i in range(4)]; k3 = rhs(s3); s4 = [s[i] + h*k3[i] for i in range(4)]; k4 = rhs(s4)
        s = [s[i] + h/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) for i in range(4)]
        rr = math.hypot(s[0], s[1])
        if rr < 0.9*r0: passed = True
        if passed and rr >= r0: break
    c1 = cons(s); phi1 = math.atan2(s[1], s[0])
    return (phi0 - phi1) % (2*math.pi), c0, c1, s
def exact_sweep(b, r0):
    """swept angle between r0 (in) and r0 (out) from (dr/dphi)^2 = r^4/b^2 - r^2 (1 - r_s/r), r_s = 1, in u = 1/r."""
    B = mp.mpf(b); u0 = 1/mp.mpf(r0)
    P = lambda u: 1/B**2 - u**2 + u**3
    utp = mp.findroot(P, 1/B)
    # substitution u = u0 + (utp - u0) sin^2 th removes the turning-point singularity
    f = lambda th: 2*(utp - u0)*mp.sin(th)*mp.cos(th)/mp.sqrt(abs(P(u0 + (utp - u0)*mp.sin(th)**2)))
    return mp.re(2*mp.quad(f, [0, mp.pi/2]))
res = {}
for rule in ('G', 'N', 'F'):
    sw, c0, c1, s_end = integrate(rule)
    beff = abs(c0[1])/c0[0]                                   # impact parameter L/omega_t of the actual initial ray
    ex = exact_sweep(beff, 1000.0)
    res[rule] = dict(sweep=sw, exact=float(ex), dOm=(c1[0] - c0[0])/c0[0], dL=(c1[1] - c0[1])/c0[1])
    print(f"  rule {rule}: swept angle {sw:.9f}  exact {float(ex):.9f}  |diff| {abs(sw - float(ex)):.2e}   drift omega_t {res[rule]['dOm']:.2e}  L {res[rule]['dL']:.2e}")
check("A-5  RULE G integrates to the exact orbit: swept angle agrees with the quadrature of the orbit identity to 1e-6 (RK4, h = 0.1 r_s, r_0 = 1000 r_s)"
      " and both pairings are conserved along the path to 1e-6", abs(res['G']['sweep'] - res['G']['exact']) < 1e-6 and abs(res['G']['dOm']) < 1e-6 and abs(res['G']['dL']) < 1e-6)
defl = lambda k: res[k]['sweep'] - math.pi
check("A-6  RULES N and F do NOT: their deflections differ from the exact one at order one (exact " + f"{res['G']['exact'] - math.pi:.6f}" +
      f", N {defl('N'):.6f}, F {defl('F'):.6f}) and their pairings drift", abs(defl('N') - defl('G')) > 0.1*abs(defl('G')) and abs(defl('F') - defl('G')) > 0.1*abs(defl('G')))
print("  READING: propagation is LOCAL, but the local law is the derivative of the pivot field along the ray's own step, not a")
print("  product of seat rotors.  Seat rotors relate seats to the frame; linking two of them against a background imports the")
print("  background -- the stopwatch -- and the number comes out wrong.  The gradient rule carries no background.")

# ---------------------------------------------------------------------------
# MASSIVE BODIES: the same two pairings with a timelike paravector, and the perihelion coefficient
# ---------------------------------------------------------------------------
print(); print("=" * 100); print("MASSIVE BODIES -- the pairings with a timelike paravector; the perihelion coefficient"); print("=" * 100)
zeta, al, lamb, rr, m = sp.symbols('zeta alpha lambda r m', positive=True)
Xm = m*para(sp.cosh(zeta), sp.sinh(zeta)*sp.Matrix([sp.cos(al), sp.sin(al), 0]))    # velocity tanh(zeta) at angle alpha to r-hat
u_s = para(sp.cosh(lamb), sp.sinh(lamb)*sp.Matrix([1, 0, 0]))
E = sp.simplify((scal(Xm)*scal(u_s) - (vect(Xm).T*vect(u_s))[0, 0])/scal(u_s))
L = sp.simplify(-(scal(Xm)*0 - (vect(Xm).T*vect(para(0, sp.Matrix([0, rr, 0]))))[0, 0]))
check("M-1  pairings for a massive body: E = m (cosh zeta - sinh zeta tanh l cos alpha) [STAT-1],  L = m r sinh zeta sin alpha [ROT-1]",
      z(E - m*(sp.cosh(zeta) - sp.sinh(zeta)*sp.tanh(lamb)*sp.cos(al))) and z(L - m*rr*sp.sinh(zeta)*sp.sin(al)))
T = sp.tanh(lamb); uu = sp.tanh(zeta)
drdphi = rr*(uu*sp.cos(al) - T)/(uu*sp.sin(al))                # PROP-1: dr/dt = u cos alpha - tanh l, r dphi/dt = u sin alpha
Et, Lt = E/m, L/m
rhs_m = rr**4/Lt**2*(Et**2 - (1 - T**2)*(1 + Lt**2/rr**2))
check("M-2  TIMELIKE ORBIT IDENTITY: (dr/dphi)^2 = (r^4/L~^2)[E~^2 - (1 - tanh^2 l)(1 + L~^2/r^2)] EXACTLY, any radial profile",
      z(sp.expand(sp.simplify(drdphi**2 - rhs_m))))
# rule G conserves E and L for the massive body too: co-moving step e = u n dt
zt = sp.Symbol('zeta', positive=True)
Xm2 = para(sp.cosh(zt), sp.sinh(zt)*n)                          # at the symbolic point (x, y), direction psi, general profile
def conserved_m(Xp, xp, yp, field):
    E_ = scal(Xp); pv = vect(Xp); vv = field.subs({x: xp, y: yp}, simultaneous=True)
    return E_ + (vv.T*pv)[0, 0], xp*pv[1] - yp*pv[0]
def drift_m(Lam_eps, field):
    uvec = sp.tanh(zt)*n
    dxv = (uvec + field)*eps; xp, yp = x + dxv[0], y + dxv[1]
    Xp = act(Lam_eps, Xm2); E1, L1 = conserved_m(Xp, xp, yp, field)
    return sp.simplify(sp.diff(E1, eps).subs(eps, 0)), sp.simplify(sp.diff(L1, eps).subs(eps, 0))
dEm, dLm = drift_m(Id - para(0, (J*(sp.tanh(zt)*n))*eps)/2, vvec)
check("M-3  RULE G with the body's own co-moving step e = u dt conserves E and L identically (general profile, any velocity, any direction)",
      z(dEm) and z(dLm))
# perihelion precession for the KIN-2a profile: P(u) = L~^2 u^3 - L~^2 u^2 + u + (E~^2 - 1), roots u_a < u_p < u_3, sum of roots = 1
def precession(rp, ra):
    up, ua = 1/mp.mpf(rp), 1/mp.mpf(ra)
    # E~^2, L~^2 from P(u_p) = P(u_a) = 0
    A_ = mp.matrix([[1, -(1 - up)*up**2], [1, -(1 - ua)*ua**2]]); rhs_ = mp.matrix([(1 - up), (1 - ua)])
    sol = mp.lu_solve(A_, rhs_); E2, L2 = sol[0], sol[1]
    u3 = 1 - up - ua
    f = lambda th: 1/mp.sqrt(u3 - ua - (up - ua)*mp.sin(th)**2)
    return 4*mp.quad(f, [0, mp.pi/2]) - 2*mp.pi, E2, L2
rp, ra = 20000, 30000; ell = 2*rp*ra/(rp + ra)
dphi, E2, L2 = precession(rp, ra)
coef = dphi*ell/mp.pi
check("M-4  perihelion advance per orbit for r_p = 2e4 r_s, r_a = 3e4 r_s: Delta phi * ell / pi = " + mp.nstr(coef, 8) + " (leading coefficient 3;"
      " E~^2 = " + mp.nstr(E2, 10) + ", bound orbit)", abs(coef - 3) < mp.mpf(10)**-3 and E2 < 1)
print("  ARGUMENT (not counted): with gamma = 1 forced by the boost (rod = 1/lapse^2), the precession coefficient (2 + 2 gamma - beta)/3 x 3")
print("  reads (4 - beta) and 3 means beta = 1: the profile's second-order term is exactly Schwarzschild's.  Part B runs the alternatives.")

print(); n_ = sum(CH); print(f"RESULT: {n_}/{len(CH)} checks passed.")
if n_ != len(CH): print("VERDICT: check failures above."); sys.exit(1)
print("VERDICT (Part A): the local transport law exists and is unique among the candidates: the momentum paravector is boosted by")
print("  the gradient of the pinned velocity field along its own co-moving step.  It conserves both symmetry pairings identically,")
print("  integrates to the exact orbit, and carries massive bodies (timelike identity M-2, perihelion coefficient 3).  Both")
print("  seat-rotor rules fail.  Conditional on KIN-2a, STAT-1, ROT-1, PROP-1.")
print("COMPARISON STAGE: rule G is the spin connection of the rain tetrad in Painleve-Gullstrand form (boost part = velocity")
print("  gradient on the co-moving displacement, rotation part zero for an irrotational flow); M-2 is the Schwarzschild timelike")
print("  geodesic equation; the coefficient 3 is 6 pi GM/(c^2 a (1 - e^2)) -- Mercury's 43 arcsec per century; rule N's error is")
print("  the river-model error of boosting neighbouring elements against a background; rule F is the static-observer relation.")
sys.exit(0)
