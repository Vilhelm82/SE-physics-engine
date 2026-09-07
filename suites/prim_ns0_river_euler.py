# prim_ns0_river_euler.py -- NS-0: is the river a Navier-Stokes flow, and does its wave sector give the light cone?
# INPUTS: DYN-1's river v = -c sqrt(r_s/r) (the free-fall rapidity). Euler/NS for a fluid; continuity; linearised sound.
# BANNED: any metric as input. The acoustic line element is BUILT from (rho, v, c_s) and only then compared.
import sympy as sp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")
r, rs, c, t, G, M = sp.symbols('r r_s c t G M', positive=True)
v = -c*sp.sqrt(rs/r)
print("=== NS-0a: the river's own inertia ===")
conv = sp.simplify(v*sp.diff(v, r))                             # (v.grad) v, radial
print("  (v.grad)v =", conv, "  with r_s = 2GM/c^2:", sp.simplify(conv.subs(rs, 2*G*M/c**2)))
check("a1", sp.simplify(conv.subs(rs, 2*G*M/c**2) + G*M/r**2)==0, "the river's convective acceleration IS Newton's -GM/r^2: Euler with p = 0 and no body force. Gravity is the flow's own inertia")
phi = sp.integrate(v, r)
check("a2", sp.simplify(sp.diff(phi, r) - v)==0, f"the river is potential flow, v = grad(phi), phi = {sp.simplify(phi)}: irrotational, so viscosity only renormalises pressure (grad of div v)")
bern = sp.simplify(v**2/2 - c**2*rs/(2*r))
check("a3", bern==0, "Bernoulli: v^2/2 = GM/r exactly -- the river is the energy integral of free fall")
divv = sp.simplify(sp.diff(r**2*v, r)/r**2)
print("  div v =", divv, " (compressive: the medium is not incompressible; continuity fixes rho(r) ~ 1/(r^2 v))")
print("=== NS-0b: the wave sector -- linearised Euler + continuity on the river ===")
# Unruh/Visser: sound in a barotropic irrotational flow sees ds^2 = (rho/c_s) [ -(c_s^2 - v^2) dt^2 - 2 v dr dt + dr^2 + r^2 dOmega^2 ]
cs = sp.symbols('c_s', positive=True)
g = sp.Matrix([[-(cs**2 - v**2), -v],[-v, 1]])                   # (t, r) block, conformal factor dropped (does not affect the cone)
# null directions of the acoustic cone: dr/dt = v +- c_s
disc = sp.solve(sp.Eq(g[0,0] + 2*g[0,1]*sp.Symbol('u') + g[1,1]*sp.Symbol('u')**2, 0), sp.Symbol('u'))
print("  acoustic null slopes dr/dt =", [sp.simplify(d) for d in disc])
check("b1", set(sp.simplify(d - v) for d in disc) == {cs, -cs}, "sound moves at +-c_s relative to the river: the cone is carried by the flow (SHEET-1's theta_+ = 0 at v = -c_s is this cone stalling)")
# set c_s = c and transform t -> T = t + int v/(c^2 - v^2) dr : is the (t,r) block Schwarzschild?
gc = g.subs(cs, c)
f = sp.simplify(v/(c**2 - v**2))
T_of_t = sp.Symbol('T'); dt_dT = 1          # dt = dT - f dr
# substitute dt = dT - f dr into ds^2 = g00 dt^2 + 2 g01 dt dr + g11 dr^2
dT, dr_ = sp.symbols('dT dr')
ds2 = gc[0,0]*(dT - f*dr_)**2 + 2*gc[0,1]*(dT - f*dr_)*dr_ + gc[1,1]*dr_**2
ds2 = sp.expand(ds2)
gTT = sp.simplify(ds2.coeff(dT, 2)); gTr = sp.simplify(ds2.coeff(dT, 1).coeff(dr_, 1)); grr = sp.simplify(ds2.coeff(dr_, 2))
print("  after the time shift: g_TT =", gTT, "  g_Tr =", gTr, "  g_rr =", grr)
check("b2", sp.simplify(gTT + c**2*(1 - rs/r))==0 and gTr==0 and sp.simplify(grr - 1/(1 - rs/r))==0,
      "with c_s = c the acoustic metric of the river IS the potential ratio's metric: g_TT = -c^2(1 - r_s/r), g_rr = 1/(1 - r_s/r). The AC sector is linearised Euler on the river")
print("  the characteristic ratio sqrt(L/C) of this sector is rho c_s: the medium's acoustic characteristic ratio. THM-K's scalar half is its VALUE.")
print("=== NS-0c: where the O(J^2) lives ===")
th, K = sp.symbols('theta K', positive=True)
vphi = K*sp.sin(th)/r**2                                       # DYN-3's swirl
centrif = sp.simplify(vphi**2/r)                               # the swirl's own inertia: (v.grad v) radial part from v_phi
print("  centrifugal term of the swirl v_phi^2/r =", centrif, " ~ K^2 sin^2(theta)/r^5 : the O(J^2) radial force the DC medium ignored")
check("c1", sp.degree(sp.denom(centrif), r)==5, "the swirl's inertia enters at 1/r^5 in force (1/r^4 in potential): a DIFFERENT radial law from Kerr's mass quadrupole (1/r^3). Kerr's Q = -J^2/M is not produced by the swirl's inertia alone")
n=sum(CH); print(f"\n=== NS-0: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
