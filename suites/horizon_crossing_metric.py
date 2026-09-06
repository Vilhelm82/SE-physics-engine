#!/usr/bin/env python3
"""Exact native Lorentzian horizon-crossing derivations.

Scope:
  * T7's seat-derived Gram form, with delta > 0.
  * A regular radial/null frame obtained from the dual normal.
  * A conditional spacetime completion with explicit clock/radial/transverse
    coframe functions; the fourth positive direction is the T8 branch.
  * Horizon observables in an explicitly stated spherical specialization.

No legacy Euclidean/Cl(3) runner, field equation, or radial mass profile is used.
Run: python -B horizon_crossing_metric.py
The JSON records exact symbolic identities, not a fitted dynamical solution.
"""
import json
import sympy as s

checks = {}


def zero(expr):
    return s.factor(s.cancel(s.simplify(expr)))


def check(name, actual, expected=0):
    residual = zero(actual - expected)
    if residual != 0:
        raise AssertionError((name, residual))
    checks[name] = {"verified": True}


def check_matrix(name, actual, expected):
    residual = actual - expected
    for item in residual:
        if zero(item) != 0:
            raise AssertionError((name, residual))
    checks[name] = {"verified": True}


# 1. A complete coframe that needs no inverse ruler-plane minor.
a, b, gamma = s.symbols("a b gamma", real=True)
G = s.Matrix([[-1, a, b], [a, 1, gamma], [b, gamma, 1]])
eps = 1 - gamma**2
delta = eps + a*a + b*b - 2*a*b*gamma
W = delta - eps
F_native = eps / delta
c, h, g = (s.eye(3)[:, i] for i in range(3))
q = lambda u, v: (u.T * G * v)[0]
nu = G.inv() * c
k_native = -nu
radial_raw = -(nu + c)
transverse_raw = b*h - a*g

check("native_determinant", G.det(), -delta)
check("dual_normal_norm", q(nu, nu), -F_native)
check_matrix("dual_normal_pairings", G*nu, c)
check("future_normal_clock_pairing", q(k_native, c), -1)
check("radial_raw_norm", q(radial_raw, radial_raw), W/delta)
check("radial_raw_seat_orthogonality", q(c, radial_raw))
check("normal_radial_pairing", q(k_native, radial_raw), W/delta)
check("transverse_raw_norm", q(transverse_raw, transverse_raw), W)
check("transverse_seat_orthogonality", q(c, transverse_raw))
check("transverse_normal_orthogonality", q(k_native, transverse_raw))

# Coframe square-root identity with rational diagonal weights.  On delta > 0,
# both positive entries of D can be square-rooted without crossing eps = 0.
A2 = 1 + a*a
P = s.Matrix([[1, -a, -b], [0, 1, (gamma+a*b)/A2], [0, 0, 1]])
D = s.diag(-1, A2, delta/A2)
check_matrix("regular_coframe_factorization", P.T*D*P, G)

# 2. Native normalized radial block and unique finite null partner in span(c,k).
# The preceding rational identities justify this orthonormal representation.
w = s.symbols("w", positive=True)
Q = s.diag(-1, 1, 1)
c0 = s.Matrix([1, 0, 0])
er = s.Matrix([0, 1, 0])
et = s.Matrix([0, 0, 1])
k0 = c0 + w*er
ell0 = (er-c0)/(1+w)
R = s.Matrix.hstack(k0, er, et)
N = s.Matrix.hstack(k0, ell0, et)
check_matrix("river_frame", R.T*Q*R,
             s.Matrix([[w*w-1, w, 0], [w, 1, 0], [0, 0, 1]]))
check_matrix("null_frame", N.T*Q*N,
             s.Matrix([[w*w-1, 1, 0], [1, 0, 0], [0, 0, 1]]))
check("null_frame_determinant", N.det(), 1)
y = s.symbols("y", real=True)
F_w = 1-w*w
x = -1-F_w*y
partner_norm = -x*x-2*x*y-F_w*y*y
y_regular = 1/(w*(1+w))
y_other = -1/(w*(1-w))
check("regular_partner_root", partner_norm.subs(y, y_regular))
check("other_partner_root", partner_norm.subs(y, y_other))
check("regular_partner_horizon_limit", s.limit(y_regular, w, 1), s.Rational(1, 2))
check("other_partner_has_horizon_pole",
      s.limit((1-w)*y_other, w, 1), -1)

# Horizon values from both selected-pole branches.
for sigma in (1, -1):
    sub = {gamma: sigma}
    depth_gap = a-sigma*b
    check(f"branch_{sigma}_determinant", delta.subs(sub), depth_gap**2)
    check(f"branch_{sigma}_signed_lapse", F_native.subs(sub))
    check(f"branch_{sigma}_drift_squared", (W/delta).subs(sub), 1)

# 3. General null-coordinate completion, allowing two transverse directions,
# transverse drift, time dependence, and a non-round transverse metric.
A, B, F = s.symbols("A B F", nonzero=True, real=True)
u1, u2, h11, h12, h22 = s.symbols("u1 u2 h11 h12 h22", real=True)
H = s.Matrix([[h11, h12], [h12, h22]])
U = s.Matrix([u1, u2])
metric = s.zeros(4)
metric[0, 0] = -A*A*F + (U.T*H*U)[0]
metric[0, 1] = metric[1, 0] = A*B
for i in range(2):
    metric[0, i+2] = metric[i+2, 0] = -(H*U)[i]
    for j in range(2):
        metric[i+2, j+2] = H[i, j]
inverse = s.zeros(4)
inverse[0, 1] = inverse[1, 0] = 1/(A*B)
inverse[1, 1] = F/B**2
for i in range(2):
    inverse[1, i+2] = inverse[i+2, 1] = U[i]/(A*B)
    for j in range(2):
        inverse[i+2, j+2] = H.inv()[i, j]
check("general_metric_determinant", metric.det(), -A*A*B*B*H.det())
check_matrix("general_metric_inverse", metric*inverse, s.eye(4))
check("horizon_metric_determinant", metric.det().subs(F, 0),
      -A*A*B*B*H.det())
fv, fr, f1, f2 = s.symbols("f_v f_r f_1 f_2", real=True)
df = s.Matrix([fv, fr, f1, f2])
angular_df = s.Matrix([f1, f2])
normal_norm = (df.T*inverse*df)[0]
check("moving_boundary_normal",
      normal_norm.subs(F, 0),
      2*fr*(fv+u1*f1+u2*f2)/(A*B) +
      (angular_df.T*H.inv()*angular_df)[0])

# 4. Spherical dynamical geometry, with r explicitly the areal radius.
v, r, theta, phi = s.symbols("v r theta phi", real=True)
Af, Bf, Ff = (s.Function(name)(v, r) for name in ("A", "B", "F"))
base = s.Matrix([[-Af**2*Ff, Af*Bf], [Af*Bf, 0]])
base_inv = s.Matrix([[0, 1/(Af*Bf)], [1/(Af*Bf), Ff/Bf**2]])
coords = (v, r)
Gamma = [[[zero(sum(base_inv[i, d] *
                   (s.diff(base[d, k], coords[j]) +
                    s.diff(base[d, j], coords[k]) -
                    s.diff(base[j, k], coords[d]))/2
                   for d in range(2)))
           for k in range(2)] for j in range(2)] for i in range(2)]
outgoing = s.Matrix([1, Af*Ff/(2*Bf)])
ingoing = s.Matrix([0, -1/(Af*Bf)])
check("spherical_outgoing_null", (outgoing.T*base*outgoing)[0])
check("spherical_ingoing_null", (ingoing.T*base*ingoing)[0])
check("spherical_null_pairing", (outgoing.T*base*ingoing)[0], -1)
exp_out = 2*outgoing[1]/r
exp_in = 2*ingoing[1]/r
check("outgoing_expansion", exp_out, Af*Ff/(Bf*r))
check("ingoing_expansion", exp_in, -2/(Af*Bf*r))

# In the stationary spherical sector, the ingoing unit-mass radial geodesic
# follows from its conserved Killing energy E and its norm.
E, z = s.symbols("E z", positive=True)
root = s.sqrt(E*E-A*A*F)
vdot = 1/(E+root)
rdot = -root/(A*B)
velocity = s.Matrix([vdot, rdot])
stationary_base = s.Matrix([[-A*A*F, A*B], [A*B, 0]])
check("ingoing_geodesic_normalization",
      (velocity.T*stationary_base*velocity)[0], -1)
check("ingoing_geodesic_conserved_energy",
      A*A*F*vdot-A*B*rdot, E)
check("ingoing_geodesic_horizon_vdot", vdot.subs(F,0), 1/(2*E))
check("ingoing_geodesic_horizon_rdot", rdot.subs(F,0), -E/(A*B))

# Freeze F's derivatives before setting F=0 at a point on its zero locus.
jet = {Ff: s.Symbol("F0")}
for total in (1, 2):
    for nv in range(total+1):
        nr = total-nv
        jet[s.diff(Ff, v, nv, r, nr)] = s.Symbol("F_"+"v"*nv+"r"*nr)
def at_horizon(expr):
    return zero(expr.xreplace(jet).subs(jet[Ff], 0))

Fv = jet[s.diff(Ff, v)]
Fr = jet[s.diff(Ff, r)]
check("outer_horizon_derivative", at_horizon(
    ingoing[1]*s.diff(exp_out, r)), -Fr/(Bf**2*r))
check("null_generator_inaffinity_at_fixed_boundary",
      at_horizon(Gamma[0][0][0]),
      s.diff(s.log(Af*Bf), v)+Af*Fr/(2*Bf))

# Direct four-dimensional Ricci_vv, independently of the warped-product formula.
g4 = s.diag(0, 0, r*r, r*r*s.sin(theta)**2)
g4[:2, :2] = base
gi4 = s.diag(0, 0, r**-2, 1/(r*r*s.sin(theta)**2))
gi4[:2, :2] = base_inv
co4 = (v, r, theta, phi)
conn = {}
def connection(i, j, k):
    key = (i, min(j,k), max(j,k))
    if key not in conn:
        conn[key] = zero(sum(gi4[i, d] *
            (s.diff(g4[d, k], co4[j]) + s.diff(g4[d, j], co4[k]) -
             s.diff(g4[j, k], co4[d]))/2 for d in range(4)))
    return conn[key]
Ricci_vv = zero(sum(
    s.diff(connection(i, 0, 0), co4[i]) -
    s.diff(connection(i, 0, i), v) +
    sum(connection(i, 0, 0)*connection(j, i, j) -
        connection(j, 0, i)*connection(i, 0, j) for j in range(4))
    for i in range(4)))
# At F=0, g_vv=0, hence Einstein_vv=Ricci_vv.
check("geometric_horizon_flux_component", at_horizon(Ricci_vv),
      -Af*Fv/(Bf*r))

# Scalar wave operator: direct connection contraction versus divergence form.
psi = s.Function("psi")(*co4)
box_direct = sum(gi4[i,j] * (
    s.diff(psi,co4[i],co4[j]) -
    sum(connection(k,i,j)*s.diff(psi,co4[k]) for k in range(4)))
    for i in range(4) for j in range(4))
box_divergence = (
    s.diff(r*r*s.diff(psi,r),v) +
    s.diff(r*r*s.diff(psi,v)+(Af/Bf)*r*r*Ff*s.diff(psi,r),r)
    )/(Af*Bf*r*r) + (
    s.diff(psi,theta,theta) + s.cot(theta)*s.diff(psi,theta) +
    s.diff(psi,phi,phi)/s.sin(theta)**2)/(r*r)
check("spherical_wave_operator", box_direct, box_divergence)

# 5. The horizon slope depends on gamma's derivative and the native depth gap.
av, bv, gv = (s.Function(name)(v, r) for name in ("a", "b", "gamma"))
dstate = 1-gv*gv+av*av+bv*bv-2*av*bv*gv
Fs = (1-gv*gv)/dstate
state_jet = {}
for f in (av, bv, gv):
    for co in (v, r):
        state_jet[s.diff(f, co)] = s.Symbol(str(f.func)+"_"+str(co))
for sigma in (1, -1):
    dr_F = zero(s.diff(Fs, r).xreplace(state_jet).subs(gv, sigma))
    dv_F = zero(s.diff(Fs, v).xreplace(state_jet).subs(gv, sigma))
    gap = av-sigma*bv
    check(f"branch_{sigma}_radial_horizon_slope", dr_F,
          -2*sigma*state_jet[s.diff(gv,r)]/gap**2)
    check(f"branch_{sigma}_temporal_horizon_slope", dv_F,
          -2*sigma*state_jet[s.diff(gv,v)]/gap**2)
    check(f"branch_{sigma}_horizon_velocity", -dv_F/dr_F,
          -state_jet[s.diff(gv,v)]/state_jet[s.diff(gv,r)])

# 6. Nontrivial crossing state: a=-1, b=0, gamma=1-rho.
# Valid on a two-sided neighbourhood of rho=0; no exterior radial law is imposed.
rho = s.symbols("rho", real=True)
delta_example = delta.subs({a:-1, b:0, gamma:1-rho})
F_example = zero(F_native.subs({a:-1, b:0, gamma:1-rho}))
check("crossing_state_full_rank", delta_example.subs(rho,0), 1)
check("crossing_state_lapse_slope", s.diff(F_example,rho).subs(rho,0), 2)
check("second_crossing_full_rank", delta_example.subs(rho,2), 1)
check("second_crossing_signed_lapse", F_example.subs(rho,2))
check("second_crossing_lapse_slope", s.diff(F_example,rho).subs(rho,2), -2)
for location, direction in ((s.Rational(1,10), 1), (-s.Rational(1,10), -1)):
    assert delta_example.subs(rho,location) > 0
    assert direction*F_example.subs(rho,location) > 0
checks["crossing_state_both_sides"] = {"verified": True}

print(json.dumps({
    "scope": "Exact native frame identities and conditional metric completion; no derived load dynamics or fixed mass profile.",
    "verified_count": len(checks),
    "checks": checks,
    "results": {
        "F": "(1-gamma^2)/(1-gamma^2+a^2+b^2-2*a*b*gamma)",
        "metric": "-A^2 F dv^2 + 2 A B dv dr + h_AB(dx^A-U^A dv)(dx^B-U^B dv)",
        "determinant": "-A^2 B^2 det(h)",
        "stationary_surface_gravity": "A_H F_r,H/(2 B_H) = -sigma A_H gamma_r,H/[B_H(a_H-sigma b_H)^2]",
        "spherical_expansions": ["A F/(B r)", "-2/(A B r)"],
        "horizon_speed": "dr_H/dv = -gamma_v/gamma_r",
        "geometric_horizon_flux": "G_vv|H = -A_H F_v,H/(B_H r_H)",
        "crossing_example": {"a":-1,"b":0,"gamma":"1-rho","F":str(F_example)}
    }
}, indent=2))
