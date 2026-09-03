#!/usr/bin/env python3
# =============================================================================
# THM-TARGET O -- the strain law: Gauss and Codazzi on the seat's flat rods, aimed at THM-N's gate.
# Date: 2026-09-03 (late) / 09-04.  sympy exact throughout; numpy only for the independent numeric path.  Exit 0 iff every
# check passes.
#
# THE CLAIM UNDER TEST (handoff 09-03, next-session item 2): the nonlinear closure THM-N demands is e_2(strain) = 0 -- the
#   second elementary symmetric function of the symmetric strain of the pivot field along the co-moving step -- sourced by
#   the blind-mass measure.  Guardrail: e_2 must be DERIVED as the unique admissible invariant from a model principle
#   before it is used, or it is ADM in a hat.
# THE PRINCIPLE: one boost per seat => the rain seat's rods are FLAT (the presented Gram of the free-fall seat is G itself;
#   the static seat's G' = G + sinh^2(l) k k^T is that flat frame seen through one boost).  A flat 3-slice carried by a
#   stationary shift v under a unit lapse has second fundamental form K = sym(grad v), and Gauss's theorema egregium fixes
#   the ambient's normal-normal Einstein component to be e_2(K) EXACTLY; Codazzi fixes the normal-tangential components to
#   (1/2)[curl curl v].  "No source through the seat's normal" is therefore e_2(sym grad v) = 0 and curl curl v = 0.  Not a
#   menu of nine generators: the Gauss invariant of flat rods is e_2 and nothing else.                             [O-1]
# THE PREREGISTERED FORK (handoff, item 2): at O(eps^2) with the rotating shell, does the closure
#   (a) cancel THM-N's P2/r term AND return Q = -J^2/(M c^2) -> Einstein at second order from the cell;
#   (b) neither -> flat frame rods fail at O(a^2) (Garat-Price), the Ricci of the presented Gram G' is the next tier;
#   (c) one of two -> name what is missing.
# NOT USED before the comparison block: the Kerr metric, Einstein's equations by name, ADM, any cancellation term, any
#   value of the quadrupole.
# CONVENTIONS: c = G = 1 inside the algebra (r_s = 2M; swirl amplitude A = 2J so that omega = A/r^3 = 2GJ/(c^2 r^3)).
#   River velocity v = -beta(r) r-hat + omega(r) (z x r) + eps^2 grad psi (thm_l / thm_m signs).  The presentation's
#   4-geometry is the eta-sandwich of the boosted flat frame: ds^2 = -dt^2 + |dx - v dt|^2 (unit lapse, flat rods).
# =============================================================================
import sys, time
import sympy as sp
import numpy as np
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
def z(e):
    e = sp.sympify(e)
    for rt in (lambda q: q, sp.simplify, lambda q: sp.simplify(sp.expand(q)), lambda q: sp.simplify(sp.powsimp(sp.expand(q), force=True))):
        try:
            if rt(e) == 0: return True
        except Exception: pass
    return False
def rnd_eval_zero(e, syms, trials=4, seed=7):
    """Exact evaluation at random rational points (independent of simplify): True iff e vanishes at every point."""
    import random
    rng = random.Random(seed)
    for _ in range(trials):
        sub = {s: sp.Rational(rng.randint(1, 40), rng.randint(1, 13)) for s in syms}
        if sp.nsimplify(sp.N(e.subs(sub), 40)) != 0 and abs(sp.N(e.subs(sub), 40)) > 1e-30: return False
    return True
T0 = time.time()
def lap(): return f"{time.time()-T0:6.1f}s"

# -----------------------------------------------------------------------------
# machinery: 4D Einstein tensor from a metric (any coordinates), and the 3+1 objects of a flat slice under a shift
# -----------------------------------------------------------------------------
def christoffel(g, ginv, X):
    n = len(X)
    return [[[sp.together(sum(ginv[a, d]*(sp.diff(g[d, c], X[b]) + sp.diff(g[d, b], X[c]) - sp.diff(g[b, c], X[d])) for d in range(n))/2)
              for c in range(n)] for b in range(n)] for a in range(n)]
def ricci(Gam, X):
    n = len(X)
    def R(b, d):
        return sum(sp.diff(Gam[a][b][d], X[a]) - sp.diff(Gam[a][b][a], X[d])
                   + sum(Gam[a][a][e]*Gam[e][b][d] - Gam[a][d][e]*Gam[e][b][a] for e in range(n)) for a in range(n))
    return sp.Matrix(n, n, lambda b, d: R(b, d))
def einstein(g, ginv, X):
    Gam = christoffel(g, ginv, X); Ric = ricci(Gam, X)
    Rs = sum(ginv[a, b]*Ric[a, b] for a in range(len(X)) for b in range(len(X)))
    return Ric - g*Rs/2, Ric, Rs

# spherical coordinates (t, r, th, ph); v given by PHYSICAL components (vr, vth, vph) in the orthonormal frame (r-hat, th-hat, ph-hat)
t, th, ph = sp.symbols('t theta phi', real=True); r = sp.symbols('r', positive=True)
X4 = [t, r, th, ph]
def pg_metric(vr, vth, vph):
    """ds^2 = -dt^2 + (dr - vr dt)^2 + (r dth - vth dt)^2 + (r sin th dph - vph dt)^2 : unit lapse, flat rods, river v."""
    # coordinate components of v: v^r = vr, v^th = vth/r, v^ph = vph/(r sin th);  h = diag(1, r^2, r^2 sin^2 th)
    h = sp.diag(1, r**2, r**2*sp.sin(th)**2)
    vup = sp.Matrix([vr, vth/r, vph/(r*sp.sin(th))])
    vdn = h*vup
    g = sp.zeros(4, 4)
    g[0, 0] = -1 + (vdn.T*vup)[0, 0]
    for i in range(3):
        g[0, i+1] = g[i+1, 0] = -vdn[i]
        for j in range(3): g[i+1, j+1] = h[i, j]
    # inverse: g^tt = -1, g^ti = -v^i, g^ij = h^ij - v^i v^j  (unit lapse, shift -v)
    hinv = h.inv()
    ginv = sp.zeros(4, 4); ginv[0, 0] = -1
    for i in range(3):
        ginv[0, i+1] = ginv[i+1, 0] = -vup[i]
        for j in range(3): ginv[i+1, j+1] = hinv[i, j] - vup[i]*vup[j]
    n_up = sp.Matrix([1, vup[0], vup[1], vup[2]])          # unit normal to t = const, future-pointing
    return g, ginv, h, hinv, vup, vdn, n_up
X3 = [r, th, ph]
def flat_slice_objects(h, hinv, vdn):
    """K_ij = (1/2)(D_i b_j + D_j b_i) for the shift one-form b_i = -v_i on the flat slice h; e_2(K); momentum vector."""
    Gam3 = christoffel(h, hinv, X3)
    b = -vdn
    Db = sp.Matrix(3, 3, lambda i, j: sp.diff(b[j], X3[i]) - sum(Gam3[k][i][j]*b[k] for k in range(3)))
    K = (Db + Db.T)/2
    Kmix = hinv*K                                            # K^i_j
    trK = sp.trace(Kmix)
    e2 = sp.together((trK**2 - sp.trace(Kmix*Kmix))/2)
    # momentum: M_i = D_j K^j_i - D_i K
    def DjKji(i):
        return sum(sp.diff(Kmix[j, i], X3[j]) + sum(Gam3[j][j][k]*Kmix[k, i] for k in range(3)) - sum(Gam3[k][j][i]*Kmix[j, k] for k in range(3)) for j in range(3))
    Mom = sp.Matrix([sp.together(DjKji(i) - sp.diff(trK, X3[i])) for i in range(3)])
    return K, Kmix, trK, e2, Mom

print("=" * 110); print("O-1  Gauss and Codazzi on flat rods: e_2(sym grad v) IS the normal-normal Einstein component, curl curl v the rest"); print("=" * 110)
VR, VTH, VPH = (sp.Function(s, real=True) for s in ('vr', 'vth', 'vph'))
vr_g, vth_g, vph_g = VR(r, th), VTH(r, th), VPH(r, th)        # general stationary axisymmetric river
g, ginv, h, hinv, vup, vdn, n_up = pg_metric(vr_g, vth_g, vph_g)
check("O-1a the presentation's 4-geometry has unit lapse and flat rods: n.n = -1, det g = -det h, g^{-1} as stated",
      z((n_up.T*g*n_up)[0, 0] + 1) and z(sp.simplify(g.det() + h.det())) and all(z(sp.simplify((g*ginv)[i, j] - (1 if i == j else 0))) for i in range(4) for j in range(4)))
print(f"  ... Einstein tensor of the general axisymmetric river (three arbitrary functions of r, theta) [{lap()}]"); sys.stdout.flush()
G4, Ric4, Rs4 = einstein(g, ginv, X4)
K, Kmix, trK, e2, Mom = flat_slice_objects(h, hinv, vdn)
Gnn = sp.together((n_up.T*G4*n_up)[0, 0])
print(f"  ... assembled [{lap()}]"); sys.stdout.flush()
d_gauss = sp.simplify(sp.expand(Gnn - e2))
check("O-1b GAUSS: G_nn = e_2(K) identically, K = sym(grad v) the second fundamental form of the flat rods (three arbitrary"
      " functions vr, vth, vph of r and theta; theorema egregium with R(h) = 0)", d_gauss == 0, f"residual {d_gauss}")
# Codazzi: G_{n i} = n^mu G_{mu i}; the momentum vector of the slice with our sign of K
Gni = sp.Matrix([sp.together(sum(n_up[m]*G4[m, i+1] for m in range(4))) for i in range(3)])
d_cod = [sp.simplify(sp.expand(Gni[i] + Mom[i])) for i in range(3)]
check("O-1c CODAZZI: G_{n i} = -(D_j K^j_i - D_i K) identically (three components; the sign is the orientation of n)", all(d == 0 for d in d_cod), f"residuals {d_cod}")
# on flat rods the momentum vector is half the double curl of v.  Verify in Cartesian coordinates (independent route).
x, y, zc, eps = sp.symbols('x y z epsilon', real=True)
Xc = sp.Matrix([x, y, zc]); rr = sp.sqrt(x**2 + y**2 + zc**2); ez = sp.Matrix([0, 0, 1])
def cart_strain(v):
    J = sp.Matrix(3, 3, lambda i, j: sp.diff(v[i], [x, y, zc][j]))
    S = (J + J.T)/2
    e2c = (sp.trace(S)**2 - sp.trace(S*S))/2
    return S, e2c
def curl(v):
    return sp.Matrix([sp.diff(v[2], y) - sp.diff(v[1], zc), sp.diff(v[0], zc) - sp.diff(v[2], x), sp.diff(v[1], x) - sp.diff(v[0], y)])
def momentum_flat(v):                                       # d_j(S_ij - delta_ij tr S) for S = sym grad v
    S, _ = cart_strain(v); trS = sp.trace(S)
    return sp.Matrix([sum(sp.diff(S[i, j], [x, y, zc][j]) for j in range(3)) - sp.diff(trS, [x, y, zc][i]) for i in range(3)])
F1, F2, F3 = (sp.Function(s, real=True) for s in ('f1', 'f2', 'f3'))
v_gen = sp.Matrix([F1(x, y, zc), F2(x, y, zc), F3(x, y, zc)])
cc = curl(curl(v_gen)); mf = momentum_flat(v_gen)
check("O-1d on flat rods d_j(S_ij - delta_ij tr S) = -(1/2) curl curl v for ANY v(x,y,z) (with K = -S this is G_{n i} = -(1/2)[curl curl v]_i): Codazzi in"
      " vacuum says the vorticity of the river is curl-free", all(z(sp.expand(mf[i] + cc[i]/2)) for i in range(3)))
print("  READING: the closure is not chosen.  Flat rods + one shift is the whole of what one boost per seat can present, and")
print("  Gauss/Codazzi say the ambient Einstein tensor seen THROUGH the seat's normal is (e_2(K); (1/2) curl curl v).  Vacuum")
print("  through the normal = e_2(sym grad v) = 0 and curl curl v = 0.  The tangential-tangential components G_ij are NOT")
print("  fixed by the seat's rods -- they are the debt this suite measures at O(eps^2) (O-6).")
print()
print("=" * 110); print("O-2  the invariant menu of the six-number cell, and the pinning"); print("=" * 110)
bet = sp.Function('beta', real=True); w = sp.Function('w', real=True)
v0 = -bet(rr)*Xc/rr                                          # the pinned inflow, any profile
S0, e2_0 = cart_strain(v0)
e1_0 = sp.simplify(sp.trace(S0)); e3_0 = sp.simplify(S0.det()); e2_0s = sp.simplify(e2_0)
rho = sp.symbols('rho', positive=True)                       # r as a plain symbol for the ODEs
def radial(e): return sp.simplify(e.subs({x: rho, y: 0, zc: 0}))
e1r, e2r, e3r = radial(e1_0), radial(e2_0s), radial(e3_0)
check("O-2a e_2(sym grad v) for the inflow v = -beta(r) r-hat is (1/r^2) d(r beta^2)/dr exactly",
      z(e2r - sp.diff(rho*bet(rho)**2, rho)/rho**2))
sol1 = sp.dsolve(sp.Eq(e1r, 0), bet(rho)); sol2 = sp.dsolve(sp.Eq(sp.diff(rho*bet(rho)**2, rho), 0), bet(rho))
check("O-2b e_1 = 0 (divergence-free river) forces beta = C/r^2 -- the WRONG exponent (clock 1/r^4): dead at redshift",
      z(sol1.rhs*rho**2 - sp.Symbol('C1')) , f"{sol1}")
check("O-2c e_3 = 0 (degenerate strain) forces beta' = 0 -- no field at all", z(e3r + sp.diff(bet(rho), rho)*bet(rho)**2/rho**2))
check("O-2d e_2 = 0 forces r beta^2 = const: the pinning beta^2 = r_s/r, i.e. THM-K's derived law, from the Gauss invariant alone",
      all(not sp.simplify(s.rhs**2*rho).has(rho) for s in (sol2 if isinstance(sol2, list) else [sol2])), f"{sol2}")
rs, c2 = sp.symbols('r_s c_2', positive=True)
f_try = rs/rho + c2*rs**2/rho**2
check("O-2e a second-order profile f = r_s/r + c_2 (r_s/r)^2 gives e_2 = -c_2 r_s^2/r^4: c_2 = 0 FORCED, Mercury becomes a check",
      z(sp.diff(rho*f_try, rho)/rho**2 + c2*rs**2/rho**4))
print("  READING: of the three elementary symmetric functions of the strain cell, only e_2 = 0 returns the pinning; e_1 = 0 and e_3 = 0")
print("  are dead on arrival.  O-1 says e_2 is not one of three but THE Gauss invariant of flat rods, so O-2 is a consistency")
print("  receipt (retrodiction against THM-K), not the derivation; the derivation is O-1b.")
print()

print("=" * 110); print("O-3  the swirl at O(eps): Codazzi gives Lense-Thirring's exponent, Gauss is untouched"); print("=" * 110)
v1 = w(rr)*ez.cross(Xc)
vfull = v0 + eps*v1
ode_target = (rho*sp.diff(w(rho), rho, 2) + 4*sp.diff(w(rho), rho))
wp_rr, wpp_rr = sp.diff(w(rho), rho).subs(rho, rr), sp.diff(w(rho), rho, 2).subs(rho, rr)
ccv = curl(curl(vfull))
check("O-3a curl curl v for v = -beta r-hat + eps w(r)(z x r) is -eps (z x r)(r w'' + 4 w')/r exactly, for ANY beta and w (Cartesian)",
      all(z(sp.simplify(ccv[i] + eps*(ez.cross(Xc))[i]*(rr*wpp_rr + 4*wp_rr)/rr)) for i in range(3)))
# the same statement through the spherical 3+1 route: momentum vector of the slice, physical components (-beta, 0, eps w r sin th)
_, _, h_s, hinv_s, _, vdn_s, _ = pg_metric(-bet(r), 0, eps*w(r)*r*sp.sin(th))
_, _, _, e2_s, Mom_s = flat_slice_objects(h_s, hinv_s, vdn_s)
check("O-3a' spherical 3+1 route (independent code): the momentum vector has only a phi-component, -(eps/2) r sin^2(theta)(r w'' + 4 w'),"
      " and e_2 = (1/r^2)(r beta^2)' - eps^2 (r w')^2 sin^2(theta)/4 exactly",
      z(sp.simplify(Mom_s[0])) and z(sp.simplify(Mom_s[1])) and z(sp.simplify(Mom_s[2] + eps*r*sp.sin(th)**2*(r*sp.diff(w(r), r, 2) + 4*sp.diff(w(r), r))/2))
      and z(sp.simplify(e2_s - sp.diff(r*bet(r)**2, r)/r**2 + eps**2*(r*sp.diff(w(r), r))**2*sp.sin(th)**2/4)))
solw = sp.dsolve(sp.Eq(ode_target, 0), w(rho))
check("O-3b Codazzi in vacuum, r w'' + 4 w' = 0, has the two-dimensional solution space w = C1 + C2/r^3; the constant is a rigid"
      " rotation of the frame at infinity (excluded), so w = A/r^3: the 1/r^3 of frame dragging from the closure, no superposition",
      z(sp.simplify(solw.rhs - sp.Symbol('C1') - sp.Symbol('C2')/rho**3)) or z(sp.simplify(solw.rhs - sp.Symbol('C1')/rho**3 - sp.Symbol('C2'))), f"{solw}")
Sf, e2_f = cart_strain(vfull); S1, e2_1 = cart_strain(v1)
check("O-3c e_2(S_0 + eps S_1) = e_2(S_0) + eps^2 e_2(S_1) exactly: the O(eps) cross term tr S_0 tr S_1 - tr(S_0 S_1) vanishes for ANY"
      " beta(r), w(r), so the swirl leaves the pinning alone at first order", z(sp.expand(e2_f - e2_0 - eps**2*e2_1)))
check("O-3d the swirl's own Gauss content is e_2(S_1) = -(r w')^2 sin^2(theta)/4; for w = A/r^3 it is -9 A^2 sin^2(theta)/(4 r^6):"
      " the strain-squared of the drag, the SOURCE of the second-order field",
      z(sp.simplify(e2_1 + wp_rr**2*(x**2 + y**2)/4)))
print("  READING: at first order the closure reproduces THM-M's result without K-6 (linear superposition): the drag profile is the")
print("  curl-free-vorticity solution of Codazzi.  At second order the drag's strain-squared enters Gauss as a source with a")
print("  monopole and a P2 part -- this is the nonlinearity THM-N asked for, and it is quadratic in the swirl, i.e. in J.")
print()

print("=" * 110); print("O-4  second order: Codazzi makes the correction a gradient; Gauss is an ANISOTROPIC operator with half-integer exponents"); print("=" * 110)
# Codazzi at O(eps^2): an axisymmetric poloidal correction with curl-free vorticity and a regular axis is a gradient
Fw = sp.Function('F', real=True)
om_phi = Fw(r, th)                                           # the phi-component of curl v_2 (v_2 poloidal => curl v_2 toroidal)
cond_r = sp.diff(sp.sin(th)*om_phi, th)/(r*sp.sin(th)); cond_th = -sp.diff(r*om_phi, r)/r     # curl(om_phi phi-hat) components
G_th = sp.Function('G', real=True)
gen = sp.dsolve(sp.Eq(sp.diff(sp.sin(th)*G_th(th), th), 0), G_th(th))
check("O-4a curl curl v_2 = 0 for poloidal axisymmetric v_2 forces r*omega_phi = G(theta) and (sin theta G)' = 0, i.e. omega_phi = C/(r sin theta):"
      " singular on the axis unless C = 0.  A regular second-order correction is CURL-FREE: v_2 = grad psi",
      z(sp.simplify(gen.rhs*sp.sin(th) - sp.Symbol('C1'))), f"{gen}")
# Gauss at O(eps^2): polarisation of e_2 around the pinned inflow, for psi = f(r) P_l(cos theta)
fpsi = sp.Function('f', real=True); ell = sp.symbols('ell', integer=True, nonnegative=True)
psi_gen = sp.Function('psi', real=True)(x, y, zc)
v2_gen = sp.Matrix([sp.diff(psi_gen, s) for s in (x, y, zc)])
S2g, _ = cart_strain(v2_gen)
polar = sp.trace(S0)*sp.trace(S2g) - sp.trace(S0*S2g)
lap_psi = sum(sp.diff(psi_gen, s, 2) for s in (x, y, zc))
psi_rr = (Xc.T*S2g*Xc)[0, 0]/rr**2
bp = sp.diff(bet(rho), rho).subs(rho, rr)
check("O-4b for ANY inflow beta(r) and ANY psi(x): tr S_0 tr S_2 - tr(S_0 S_2) = -(beta' + beta/r) lap psi + (beta' - beta/r) psi_rr",
      z(sp.expand(polar - (-(bp + bet(rr)/rr)*lap_psi + (bp - bet(rr)/rr)*psi_rr))))
# specialise to the pinning beta = sqrt(r_s/r) and psi = f(r) P_l(cos theta) for l = 0..3  (spherical 3+1 route: fast and exact)
ops = {}
for L in range(4):
    psiL = fpsi(r)*sp.legendre(L, sp.cos(th))
    _, _, hL, hinvL, _, vdnL, _ = pg_metric(-sp.sqrt(rs/r) + eps**2*sp.diff(psiL, r), eps**2*sp.diff(psiL, th)/r, 0)
    _, _, _, e2L, _ = flat_slice_objects(hL, hinvL, vdnL)
    polL = sp.expand(e2L).coeff(eps, 2)
    target = -(sp.sqrt(rs/r)/(2*r))*(4*sp.diff(fpsi(r), r, 2) + 2*sp.diff(fpsi(r), r)/r - L*(L+1)*fpsi(r)/r**2)*sp.legendre(L, sp.cos(th))
    ops[L] = z(sp.simplify(polL - target))
check("O-4c with the pinning, the polarised Gauss operator on psi = f(r) P_l is -(beta/2r)[4 f'' + 2 f'/r - l(l+1) f/r^2] P_l for l = 0,1,2,3:"
      " NOT the Laplacian -- the radial direction is weighted four to one (spherical 3+1 route)", all(ops.values()), f"{ops}")
p = sp.symbols('p')
indicial = sp.solve(4*p*(p - 1) + 2*p - ell*(ell + 1), p)
check("O-4d the indicial exponents of that operator are p = (l+1)/2 and p = -l/2 -- HALF the Laplace exponents l and -(l+1)",
      set(sp.simplify(q) for q in indicial) == {sp.Rational(1, 2)*(ell + 1), -ell/2}, f"{indicial}")
# what a Newtonian clock multipole would need: |v|^2 = beta^2 + 2 eps^2 v_0.grad psi, and 2 v_0.grad psi ~ P_l r^{-(l+1)} needs psi ~ r^{1/2 - l}
newt = sp.Rational(1, 2) - ell
coincide = [L for L in range(0, 8) if any(sp.simplify(q.subs(ell, L) - newt.subs(ell, L)) == 0 for q in indicial)]
check("O-4e a Newtonian clock multipole P_l/r^{l+1} needs psi ~ r^{1/2 - l}; that exponent is in the flat-rod spectrum for l = 0 and l = 1 ONLY."
      " Flat rods carry a mass and a boost and NO mass multipole of any higher order", coincide == [0, 1], f"coincidences at l = {coincide}")
check("O-4f in particular NO psi ~ P2 r^{1/2} mode exists, so no P2/r term in the clock can arise at O(eps^2): THM-N's Gate 1 is passed"
      " by the STRUCTURE of the closure, with nothing imported", all(sp.simplify(q.subs(ell, 2) - sp.Rational(1, 2)) != 0 for q in indicial))
print("  READING: the second-order correction obeys an anisotropic second-order equation whose exterior solutions are r^{-l/2} P_l.")
print("  The pinning's own l = 0 mode is r^{1/2} (the mass), the l = 1 mode is r (a uniform velocity, a boost of the frame).  For")
print("  l >= 2 the exponent is half-integer and NOT Newton's.  This is a theorem about flat rods, not about e_2: the Gauss")
print("  invariant is doing exactly what Gauss says, and what it says is that a flat slice cannot be bent into a quadrupole.")
print()

print("=" * 110); print("O-5  the exterior at O(eps^2), and the invariant that reads a quadrupole"); print("=" * 110)
A, a0, b2, M = sp.symbols('A a_0 b_2 M', real=True)
c0 = -A**2/(10*sp.sqrt(rs)); c2s = A**2/(8*sp.sqrt(rs))
P2 = sp.legendre(2, sp.cos(th))
psi2 = a0*sp.sqrt(r) + c0*r**sp.Rational(-5, 2) + (b2/r + c2s*r**sp.Rational(-5, 2))*P2
vr_e = -sp.sqrt(rs/r) + eps**2*sp.diff(psi2, r); vth_e = eps**2*sp.diff(psi2, th)/r; vph_e = eps*(A/r**3)*r*sp.sin(th)
gE, ginvE, hE, hinvE, vupE, vdnE, nE = pg_metric(vr_e, vth_e, vph_e)
KE, KmixE, trKE, e2E, MomE = flat_slice_objects(hE, hinvE, vdnE)
e2_ser = [sp.simplify(sp.diff(e2E, eps, k).subs(eps, 0)/sp.factorial(k)) for k in range(3)]
mom_ser = [[sp.simplify(sp.diff(MomE[i], eps, k).subs(eps, 0)/sp.factorial(k)) for k in range(3)] for i in range(3)]
check("O-5a with c_0 = -A^2/(10 sqrt r_s) and c_2 = A^2/(8 sqrt r_s), Gauss e_2 = 0 and Codazzi = 0 hold through O(eps^2) for EVERY a_0, b_2:"
      " the two homogeneous constants (mass renormalisation, l = 2 mode) are exactly the freedom the closure leaves",
      all(q == 0 for q in e2_ser) and all(q == 0 for row in mom_ser for q in row), f"e2 by order {e2_ser}; momentum {mom_ser}")
psi2w = a0*sp.sqrt(r) + sp.Rational(11, 10)*c0*r**sp.Rational(-5, 2) + (b2/r + sp.Rational(11, 10)*c2s*r**sp.Rational(-5, 2))*P2
_, _, hW, hinvW, _, vdnW, _ = pg_metric(-sp.sqrt(rs/r) + eps**2*sp.diff(psi2w, r), eps**2*sp.diff(psi2w, th)/r, vph_e)
_, _, _, e2W, _ = flat_slice_objects(hW, hinvW, vdnW)
wrong2 = sp.simplify(sp.diff(e2W, eps, 2).subs(eps, 0)/2)
check("O-5b perturbing the particular constants by 10% leaves e_2 = O(eps^2) != 0: c_0 and c_2 are FORCED by the swirl's strain-squared",
      wrong2 != 0, f"residual {sp.factor(wrong2)}")
# the clock: |v|^2, the static seat's h_00 = 1 - (d tau/dt)^2, through O(eps^2); powers of r made integer with r = s^2
s_ = sp.symbols('s', positive=True); mu = sp.symbols('mu', real=True)
V2 = sp.simplify(sp.diff((vdnE.T*vupE)[0, 0], eps, 2).subs(eps, 0)/2)
V2s = sp.expand(sp.simplify(V2.subs(sp.cos(th), mu).subs(sp.sin(th)**2, 1 - mu**2)).subs(r, s_**2))
coef = {k: sp.factor(V2s.coeff(s_, -k)) for k in (2, 5, 6, 8)}
check("O-5c the O(eps^2) clock: -sqrt(r_s) a_0/r (mass), 2 sqrt(r_s) b_2 P2 r^{-5/2} (the free l = 2 mode), [A^2/6 - A^2 P2/24]/r^4 (forced);"
      " NO P2/r (THM-N's ghost) and NO P2/r^3 (a Newtonian quadrupole): the P2 content sits at r^{-5/2} and r^{-4}",
      z(coef[2] + sp.sqrt(rs)*a0) and z(coef[5] - 2*sp.sqrt(rs)*b2*sp.legendre(2, mu)) and coef[6] == 0
      and z(coef[8] - (A**2/6 - A**2*sp.legendre(2, mu)/24)), f"{coef}")
print("  h_00 at O(eps^2) =", sp.collect(V2s.subs(s_, sp.sqrt(r)), [a0, b2]))
# the invariant: equatorial circular orbits of the presentation (transport = geodesic, THM-I/THM-L), frequency at infinity vs
# proper circumference.  Flat rods: circumference = 2 pi r exactly, so R_c = r.
Veq = sp.simplify((vdnE.T*vupE)[0, 0].subs(th, sp.pi/2))
gtt = -1 + Veq; gtp = sp.simplify(gE[0, 3].subs(th, sp.pi/2)); gpp = sp.simplify(gE[3, 3].subs(th, sp.pi/2))
Om = sp.symbols('Omega')
circ = sp.diff(gtt, r) + 2*Om*sp.diff(gtp, r) + Om**2*sp.diff(gpp, r)
roots = sp.solve(circ, Om)
avg = sp.simplify(sum(rt**2 for rt in roots)/2)
avg2 = sp.expand(sp.simplify(sp.diff(avg, eps, 2).subs(eps, 0)/2).subs(r, s_**2))
avg0 = sp.simplify(avg.subs(eps, 0))
inv = {k: sp.factor(avg2.coeff(s_, -k)) for k in (6, 9, 10, 12)}
check("O-5d the invariant (1/2)(Omega_+^2 + Omega_-^2) as a function of R_c: r_s/(2R^3) at zeroth order; at O(eps^2) a mass shift"
      " -sqrt(r_s) a_0/(2R^3), a TAIL -(5/4) sqrt(r_s) b_2 R^{-9/2}, ZERO at R^{-5}, and (7/8) A^2 R^{-6}",
      z(avg0 - rs/(2*r**3)) and z(inv[6] + sp.sqrt(rs)*a0/2) and z(inv[9] + sp.Rational(5, 4)*sp.sqrt(rs)*b2) and inv[10] == 0
      and z(inv[12] - sp.Rational(7, 8)*A**2), f"{inv}")
print("  READING: a mass quadrupole is read by orbits as an R^{-5} term in this invariant (Kepler's law with a 1/R^2 correction).  The")
print("  flat-rod exterior has none: its l = 2 freedom shows up as a half-integer tail R^{-9/2}, which no multipole of a source")
print("  produces.  The strain law on flat rods can be FED any source and will never return a quadrupole.")
print()

print("=" * 110); print("O-6  what the seat's rods do NOT fix: the tangential Einstein components at O(eps^2)"); print("=" * 110)
print(f"  ... Einstein tensor of the explicit O(eps^2) field [{lap()}]"); sys.stdout.flush()
G4E, _, _ = einstein(gE, ginvE, X4)
def order2(e): return sp.simplify(sp.diff(e, eps, 2).subs(eps, 0)/2)
def order1(e): return sp.simplify(sp.diff(e, eps).subs(eps, 0))
def order0(e): return sp.simplify(e.subs(eps, 0))
print(f"  ... assembled [{lap()}]"); sys.stdout.flush()
G0 = [[order0(G4E[i, j]) for j in range(4)] for i in range(4)]
G1 = [[order1(G4E[i, j]) for j in range(4)] for i in range(4)]
check("O-6a at O(eps^0) the presentation is vacuum: G_mu nu = 0 (the pinning is exact, all ten components)", all(q == 0 for row in G0 for q in row))
check("O-6b at O(eps^1) the presentation is vacuum: G_mu nu = 0 (flat rods + curl-free drag A/r^3 is an exact first-order field,"
      " all ten components -- THM-L/M's frame dragging is not just constraint-satisfying, it is vacuum)", all(q == 0 for row in G1 for q in row))
# G_nn = e_2(K) and G_ni = -(momentum vector) are IDENTITIES (O-1b, O-1c), and O-5a showed both vanish through O(eps^2):
check("O-6c at O(eps^2) the components the rods DO fix vanish: G_nn = e_2 = 0 (Gauss) and G_ni = -Mom = 0 (Codazzi), for every a_0, b_2"
      " -- corollary of the identities O-1b/c and the expansion O-5a", all(q == 0 for q in e2_ser) and all(q == 0 for row in mom_ser for q in row))
Gij2 = {(i, j): order2(G4E[i+1, j+1]) for i in range(3) for j in range(i, 3)}
# the tangential components carry a residual that no choice of (a_0, b_2) removes: collect in s = sqrt(r) and cos theta
def free_part(e): return sp.expand(sp.simplify(e.subs(sp.cos(th), mu).subs(sp.sin(th)**2, 1 - mu**2)).subs(r, s_**2))
resid = {k: free_part(v_) for k, v_ in Gij2.items()}
const_terms = {k: sp.factor(v_.subs({a0: 0, b2: 0})) for k, v_ in resid.items()}
# explicit linear algebra (sp.solve is not trusted here): G_ij^(2) = C_ij + a_0 A_ij + b_2 B_ij; a_0 must drop out (a mass shift is
# still vacuum); solve b_2 from G_thth at one point, then show the residual survives at other points and in other components
Acoef = {k: sp.simplify(sp.diff(v_, a0)) for k, v_ in Gij2.items()}
diag = [(0, 0), (1, 1), (2, 2)]
pt1, pt2, pt3 = {r: sp.Rational(7, 2), th: sp.pi/3}, {r: sp.Rational(23, 4), th: sp.pi/5}, {r: 9, th: 2*sp.pi/7}
eq1 = sp.simplify(Gij2[(1, 1)].subs(pt1))
b2_star = sp.solve(eq1, b2)
resid_after = [sp.simplify(Gij2[k].subs(pt).subs(b2, b2_star[0])) for k in diag for pt in (pt2, pt3)] if b2_star else None
check("O-6d at O(eps^2) the tangential components G_ij are NOT zero: a_0 drops out of all of them (a mass shift is still vacuum), and the"
      " b_2 that kills G_thth at one point leaves G_rr, G_thth, G_phph nonzero at two others.  Flat rods cannot be vacuum at second"
      " order in the spin, for ANY value of the two free constants",
      all(q == 0 for q in Acoef.values()) and bool(b2_star) and resid_after is not None and any(q != 0 for q in resid_after),
      f"a_0-coefficients {set(Acoef.values())}; b_2* from G_thth at (7/2, pi/3) = {b2_star}; residuals elsewhere {[sp.N(q, 4) for q in resid_after]}")
print("  G_rr, G_thth, G_phph at O(eps^2) (a_0 = b_2 = 0):")
for k in ((0, 0), (1, 1), (2, 2)): print(f"    {k}: {sp.collect(const_terms[k].subs(s_, sp.sqrt(r)), A)}")
print("  READING: Gauss and Codazzi are the WHOLE of what a flat slice can say about its ambient.  They are satisfied.  The six")
print("  tangential components are the ambient's OWN equations -- the ones the seat's rods do not see -- and at O(eps^2) they")
print("  fail for every value of the two free constants.  That is fork (b): the rods must bend at second order in J, and the")
print("  bending is a rank-two object, which one boost cannot present.")
print("  Scaffold register: NORMAL-1 (declared; Einstein's constraint renamed -- removal route: the blind-mass theorem, sourced Gauss")
print("  from MASS-1); K-6 linearity (declared; the integer-multipole kill of the tail rests on it, its GROUND is J_2's a^{-7/2});")
print("  drag amplitude 2 = 1 + 1 (THM-M superposition; removal route: Codazzi sourced by the shell's momentum density, not run here).")
print()

print("=" * 110); print("O-7  independent numeric path: the assembled field, finite differences at 40 digits, no symbolic derivative reused"); print("=" * 110)
import mpmath as mp, random
mp.mp.dps = 40
def field_factory(rsn, An, a0n, b2n, scale_c=1):
    c0n = -An**2/(10*mp.sqrt(rsn))*scale_c; c2n = An**2/(8*mp.sqrt(rsn))*scale_c
    def psi(P):
        xn, yn, zn = P; rn = mp.sqrt(xn**2 + yn**2 + zn**2); ct = zn/rn; P2n = (3*ct**2 - 1)/2
        return a0n*mp.sqrt(rn) + c0n*rn**mp.mpf('-2.5') + (b2n/rn + c2n*rn**mp.mpf('-2.5'))*P2n
    def grad(f, P, hh):
        out = []
        for i in range(3):
            Pp = list(P); Pm = list(P); Pp[i] += hh; Pm[i] -= hh
            out.append((f(Pp) - f(Pm))/(2*hh))
        return out
    def v(P, e):
        xn, yn, zn = P; rn = mp.sqrt(xn**2 + yn**2 + zn**2); b = mp.sqrt(rsn/rn)
        v0 = [-b*xn/rn, -b*yn/rn, -b*zn/rn]
        wn = An/rn**3; v1 = [-wn*yn, wn*xn, mp.mpf(0)]
        g = grad(psi, P, mp.mpf('1e-9'))
        return [v0[i] + e*v1[i] + e**2*g[i] for i in range(3)]
    return v, grad
def jac(vf, P, hh):
    J = [[None]*3 for _ in range(3)]
    for j in range(3):
        Pp = list(P); Pm = list(P); Pp[j] += hh; Pm[j] -= hh
        vp, vm = vf(Pp), vf(Pm)
        for i in range(3): J[i][j] = (vp[i] - vm[i])/(2*hh)
    return J
def e2_num(vf, P, hh):
    J = jac(vf, P, hh); S = [[(J[i][j] + J[j][i])/2 for j in range(3)] for i in range(3)]
    tr = S[0][0] + S[1][1] + S[2][2]; tr2 = sum(S[i][j]*S[j][i] for i in range(3) for j in range(3))
    return (tr**2 - tr2)/2
def curl_of(vf, hh):
    def c(P):
        J = jac(vf, P, hh)
        return [J[2][1] - J[1][2], J[0][2] - J[2][0], J[1][0] - J[0][1]]
    return c
rng = random.Random(2026)
pts = [[mp.mpf(rng.uniform(-6, 6)) for _ in range(3)] for _ in range(3)]
pts = [P for P in pts if mp.sqrt(sum(q**2 for q in P)) > 3] or [[mp.mpf(4), mp.mpf(1), mp.mpf(2)]]
rsn, An, a0n, b2n = mp.mpf(1), mp.mpf('0.7'), mp.mpf('0.3'), mp.mpf('-0.5')
vF, _ = field_factory(rsn, An, a0n, b2n)
vW, _ = field_factory(rsn, An, a0n, b2n, scale_c=mp.mpf(2))
ratios, ratios_w, ccs = [], [], []
for P in pts:
    e_big, e_small = mp.mpf('0.2'), mp.mpf('0.1')
    q1 = e2_num(lambda Q: vF(Q, e_big), P, mp.mpf('1e-8')); q2 = e2_num(lambda Q: vF(Q, e_small), P, mp.mpf('1e-8'))
    ratios.append(q1/q2)
    w1 = e2_num(lambda Q: vW(Q, mp.mpf('0.04')), P, mp.mpf('1e-8')); w2 = e2_num(lambda Q: vW(Q, mp.mpf('0.02')), P, mp.mpf('1e-8'))
    ratios_w.append(w1/w2)
    cc_f = curl_of(curl_of(lambda Q: vF(Q, e_big), mp.mpf('1e-7')), mp.mpf('1e-6'))
    ccs.append(max(abs(q) for q in cc_f(P)))
check("O-7a numeric (mpmath, nested central differences, 40 digits): halving eps divides e_2 by 16 at three random exterior points --"
      " the Gauss residual is O(eps^4): the O(eps^2) closure holds with the forced constants", all(abs(q - 16) < mp.mpf('0.5') for q in ratios),
      "ratios " + ", ".join(mp.nstr(q, 6) for q in ratios))
check("O-7b with the particular constants doubled the ratio (eps 0.04 vs 0.02) drops to ~4: e_2 = O(eps^2), the constants are not decorative",
      all(abs(q - 4) < mp.mpf('0.6') for q in ratios_w), "ratios " + ", ".join(mp.nstr(q, 6) for q in ratios_w))
check("O-7c curl curl v vanishes to finite-difference precision at every order (Codazzi is linear in v; both the drag and the gradient"
      " correction are curl-curl-free exactly)", all(q < mp.mpf('1e-9') for q in ccs), "max |curl curl v| " + ", ".join(mp.nstr(q, 3) for q in ccs))
print()

print("=" * 110); print("COMPARISON BLOCK (names allowed from here on)"); print("=" * 110)
a_k = sp.symbols('a', real=True); Rc = sp.symbols('R', positive=True); Mk = sp.symbols('M', positive=True); Jk = sp.symbols('J', real=True)
rb = sp.symbols('r_BL', positive=True)
Om_p = sp.sqrt(Mk)/(rb**sp.Rational(3, 2) + a_k*sp.sqrt(Mk)); Om_m = -sp.sqrt(Mk)/(rb**sp.Rational(3, 2) - a_k*sp.sqrt(Mk))
avgK = sp.simplify((Om_p**2 + Om_m**2)/2)
# proper equatorial circumference of Kerr: 2 pi sqrt(g_phiphi) with g_phiphi = r^2 + a^2 + 2 M a^2/r; invert R(r) to O(a^2)
u2 = sp.symbols('u2')
rb_of_R = Rc + a_k**2*u2
eqR = sp.expand(sp.series((rb_of_R**2 + a_k**2 + 2*Mk*a_k**2/rb_of_R - Rc**2), a_k, 0, 3).removeO())
u2_sol = sp.solve(eqR.coeff(a_k, 2), u2)[0]
avgK_R = sp.series(avgK.subs(rb, Rc + a_k**2*u2_sol), a_k, 0, 3).removeO()
avgK_R = sp.expand(sp.simplify(avgK_R))
kerr_coeffs = {k: sp.factor(avgK_R.subs(a_k, Jk/Mk).coeff(Rc, -k)) for k in (3, 5, 6)}
check("C-1 Kerr, same invariant, same variable (frequency at infinity vs proper circumference): M/R^3 + (3 J^2/2M) R^{-5} + 6 J^2 R^{-6} + O(a^4)."
      " The quadrupole Q = -J^2/M is the R^{-5} coefficient, -3Q/(2M) per R^2 on Kepler",
      z(kerr_coeffs[3] - Mk) and z(kerr_coeffs[5] - sp.Rational(3, 2)*Jk**2/Mk) and z(kerr_coeffs[6] - 6*Jk**2), f"{kerr_coeffs}")
model_coeffs = {k: v_.subs({rs: 2*Mk, A: 2*Jk}) for k, v_ in inv.items()}
print("  model (flat rods, r_s = 2M, A = 2J, eps absorbed):  R^-3:", sp.simplify(Mk + model_coeffs[6]), " R^-9/2:", model_coeffs[9],
      " R^-5:", model_coeffs[10], " R^-6:", model_coeffs[12])
print("  Kerr:                                                R^-3:", kerr_coeffs[3], " R^-9/2: 0  R^-5:", kerr_coeffs[5], " R^-6:", kerr_coeffs[6])
check("C-2 GATE 1 (THM-N): no P2/r mass aspect -- PASSED by the closure (O-4f, O-5c).  GATE 2 (Q = -J^2/M): NOT REACHABLE -- the model's"
      " R^{-5} coefficient is identically zero and its l = 2 freedom is a half-integer tail Kerr does not have.  Fork (b).",
      model_coeffs[10] == 0 and kerr_coeffs[5] != 0)
print("  BAND, named loudly: O-1 on flat rods with unit lapse IS the ADM Hamiltonian and momentum constraint pair in Painleve-Gullstrand")
print("  gauge; O-2 is Schwarzschild in PG form; O-3 and O-6b are the Doran / Hamilton-Lisle river at first order in J.  That is the")
print("  band this suite reproduces.  BAND EDGE (the model's own statement, no GR name): flat rods carry a mass and a boost and no")
print("  other multipole, for ANY source -- O-4e.  Garat & Price 2000 (PRD 61, 124011; no conformally flat axisymmetric slice of Kerr,")
print("  obstruction at O(a^2)) and Valiente Kroon 2004 are the shadow of that edge, narrower and about Kerr only; Natario 2009's")
print("  PG-like Kerr coordinates have curved slices.  THM-N's Gate 2 as written -- Q = -J^2/M FROM A SHELL -- is Kerr's")
print("  no-hair number, not a rigidly rotating shell's: de la Cruz & Israel 1968 and Pfister & Braun 1985 need an oblate, stressed")
print("  shell to source Kerr at O(a^2).  The gate to re-declare is: (i) a hole -- Q fixed by regularity at the wall (THM-H at O(J^2));")
print("  (ii) a shell -- Q fixed by the source, compared with GR's shell, not with Kerr.  On flat rods neither question arises.")
print("  COROLLARY (same operator, static source): a static oblate body's l = 2 clock would fall as P2 r^{-5/2} on flat rods, not")
print("  P2 r^{-3}; the J_2 nodal regression would scale as a^{-3} instead of Newton's a^{-7/2} (LAGEOS/GPS/lunar: a^{-7/2}).  This")
print("  kills FLAT RODS for any non-spherical source, not the model: THM-M's static rods are a sum of rank-ones and are not flat.")
print("  It locates the multipoles: they live in the bending of the rods, never in the strain of a flat slice.")
print()
print("=" * 110); print("VERDICT (THM-O)"); print("=" * 110)
print("  e_2(strain) = 0 is DERIVED: it is Gauss's theorema egregium for the seat's flat rods (G_nn = e_2(K) identically, O-1b), and")
print("  Codazzi supplies curl curl v = 0 (O-1c/d).  Conditional on: one boost per seat (flat rain rods, unit lapse), the presentation")
print("  metric = eta-sandwich of the boosted frame (BARE-1/M-2), and 'no source through the seat's normal' in vacuum (declared: NORMAL-1).")
print("  Retrodictions: the pinning (O-2d, = THM-K), c_2 = 0 (O-2e), Lense-Thirring's 1/r^3 (O-3b, = THM-M without K-6), and -- new --")
print("  the first-order drag on flat rods is exact vacuum, all ten components (O-6b).")
print("  Second order (proved): Codazzi makes the correction a gradient; Gauss is the anisotropic operator with exponents (l+1)/2, -l/2.")
print("  Gate 1 PASSED by structure: no P2/r mode exists (O-4f, O-5c).  Gate 2 NOT REACHABLE: no R^{-5} slot in the orbital invariant")
print("  (O-5d, C-1/C-2); the l = 2 freedom is a half-integer tail, killed on the GROUND by J_2's scaling and, provisionally, by K-6's")
print("  linearity (a scaffold).  The tangential Einstein components fail at O(eps^2) for every")
print("  (a_0, b_2) (O-6d).  FORK (b): flat rods fail at second order in the spin.  The multipoles of any source live in the bending")
print("  of the rods; the bending is a rank-two object; one boost presents rank one.  Next tier: what bends the rods -- the same")
print("  rank-two debt the TT wave carries (handoff item 4).  Kill condition for THIS suite: any of O-1b/c, O-4d, O-5d, O-6d failing")
print("  on re-run, or a flat-rod field with a nonzero R^{-5} coefficient.")
n_pass = sum(CH); print(f"\n{n_pass}/{len(CH)} checks passed [{lap()}]")
sys.exit(0 if all(CH) else 1)
