#!/usr/bin/env python3
"""Exact checks for the implicit normal-geometry extension.

The proof and geometric hypotheses are in
docs/2026-09-05-cella-normal-geometry.md.  These finite examples check the
construction, gauge law, mixed elementary invariants and normal curvature;
they do not select the native physical cut or transport law.
"""
import json
from pathlib import Path

import sympy as s

from cella_tensor_valuation import CurvatureNumerators


checks = {}
witnesses = {}


def simp(value):
    return s.factor(s.cancel(value))


def check(name, actual, expected=0):
    if isinstance(actual, s.MatrixBase):
        if not isinstance(expected, s.MatrixBase):
            expected = s.zeros(*actual.shape)
        residual = [simp(v) for v in actual - expected]
    else:
        residual = [simp(actual - expected)]
    if any(v != 0 for v in residual):
        raise AssertionError((name, residual))
    checks[name] = {"verified": True}


def require(name, condition):
    if not condition:
        raise AssertionError(name)
    checks[name] = {"verified": True}


def diagonal_part(matrix):
    return s.diag(*matrix.diagonal())


def geometry(metric, derivative):
    j = metric.inv() * derivative.T
    c = derivative * j
    w = j * c.inv()
    p = s.eye(metric.rows) - w * derivative
    return j, c, w, p


# A codimension-two graph in 1+4 dimensions, with nonzero slopes, genuinely
# mixed Hessians and noncommuting shape operators.  Nothing is spherical or
# simultaneously diagonalized.  All arithmetic is exact.
t, wcoord, x, y, z = s.symbols("t w x y z", real=True)
coords = (t, wcoord, x, y, z)
uv = (x, y, z)
metric = s.diag(-1, 1, 1, 1, 1)
f = x/s.Integer(4) + y/s.Integer(5) + z/s.Integer(6)
f += (x*x + 2*x*y + 3*y*y + 4*x*z + 2*y*z + 5*z*z)/2
f += x*x*y + y*z*z
gfun = x/s.Integer(7) - y/s.Integer(8) + z/s.Integer(9)
gfun += (2*x*x - 2*x*y + y*y + 6*x*z - 4*y*z - 3*z*z)/2
gfun += x*y*z
relation = s.Matrix([t-f, wcoord-gfun])
embedding = s.Matrix([f, gfun, x, y, z])
origin = dict.fromkeys(coords, 0)
d = relation.jacobian(coords).subs(origin)
hess = [s.hessian(a, coords).subs(origin) for a in relation]
tangent = embedding.jacobian(uv).subs(origin)
h = tangent.T * metric * tangent
hinv = h.inv()
j, c, w, p = geometry(metric, d)
check("cross_graph_tangent_kernel", d*tangent)
check("cross_graph_projector_idempotent", p*p, p)
check("cross_graph_projector_metric_self_adjoint", p.T*metric, metric*p)
check("cross_graph_projector_tangent", p*tangent, tangent)
check("cross_graph_normal_right_inverse", d*w, s.eye(2))
check("cross_graph_normal_metric", w.T*metric*w, c.inv())
require("cross_graph_spacelike", all(h[:i, :i].det() > 0 for i in range(1, 4)))
require("cross_graph_lorentzian_normal", c.det() < 0)
restricted = [tangent.T*hh*tangent for hh in hess]
for aa in range(3):
    for bb in range(aa, 3):
        implicit_b = -w*s.Matrix([ee[aa, bb] for ee in restricted])
        graph_b = (s.eye(5)-p)*embedding.diff(uv[aa], uv[bb]).subs(origin)
        check(f"implicit_vs_embedded_B_{aa}{bb}", implicit_b, graph_b)

# An arbitrary normal uses n=J lambda.  Its shape is linear in lambda; the
# ambient projected formula has the same nonzero characteristic factor.
l0, l1, zz, tc, ts = s.symbols("lambda0 lambda1 zz tc ts", real=True)
lam = s.Matrix([l0, l1])
normal = j*lam
shape_basis = [-hinv*ee for ee in restricted]
shape = l0*shape_basis[0] + l1*shape_basis[1]
require("cross_graph_noncommuting_shapes", shape_basis[0]*shape_basis[1] != shape_basis[1]*shape_basis[0])
combined_hess = l0*hess[0] + l1*hess[1]
ambient_shape = -p*metric.inv()*combined_hess*p
check("shape_ambient_restriction", ambient_shape*tangent, tangent*shape)
check("shape_ambient_normal_kernel", ambient_shape*w)
check("shape_metric_self_adjoint", h*shape, shape.T*h)
# Comparing power traces through dimension three verifies all characteristic
# coefficients without expanding a redundant five-dimensional determinant.
for power in (1, 2, 3):
    check(f"ambient_shape_power_trace_{power}", s.trace(ambient_shape**power), s.trace(shape**power))
char = s.Poly((s.eye(3)+zz*shape).det(), zz)
check("shape_sigma1", char.nth(1), s.trace(shape))
check("shape_sigma2", char.nth(2), (s.trace(shape)**2-s.trace(shape*shape))/2)
check("shape_sigma3", char.nth(3),
      (s.trace(shape)**3-3*s.trace(shape)*s.trace(shape*shape)+2*s.trace(shape**3))/6)
mixed_det = s.Poly(shape.det(), l0, l1)
require("all_cubic_normal_monomials_present", all(mixed_det.coeff_monomial(l0**i*l1**(3-i)) != 0 for i in range(4)))
bordered = s.zeros(2).row_join(d).col_join(d.T.row_join(metric-zz*(hess[0]+2*hess[1])))
check("codimension_two_bordered_fingerprint", bordered.det()/(metric.det()*c.det()),
      (s.eye(3)+zz*(shape_basis[0]+2*shape_basis[1])).det())

# Full A(x)F gauge, including normal mixing and nonconstant first derivatives.
gauge = s.Matrix([[2+x+2*wcoord, 1+y+t], [1+z-wcoord, 1+y-x]])
gauged_relation = gauge*relation
a0 = gauge.subs(origin)
dt = gauged_relation.jacobian(coords).subs(origin)
ht = [s.hessian(a, coords).subs(origin) for a in gauged_relation]
jt, ct, wt, pt = geometry(metric, dt)
check("gauge_derivative", dt, a0*d)
check("gauge_gram", ct, a0*c*a0.T)
check("gauge_normal_frame", wt, w*a0.inv())
check("gauge_projector", pt, p)
gauge_term = []
for aa in range(2):
    a_hess = sum((a0[aa, bb]*hess[bb] for bb in range(2)), s.zeros(5))
    expected = s.zeros(5)
    for bb in range(2):
        da = s.Matrix([s.diff(gauge[aa, bb], cc).subs(origin) for cc in coords])
        expected += da*d[bb, :] + d[bb, :].T*da.T
    gauge_term.append(expected)
    check(f"gauge_hessian_{aa}", ht[aa]-a_hess, expected)
    check(f"gauge_hessian_tangent_cancellation_{aa}", tangent.T*expected*tangent)

lambda_tilde = a0.T.inv()*lam
gauged_shape = -hinv*sum((lambda_tilde[aa]*(tangent.T*ht[aa]*tangent)
                         for aa in range(2)), s.zeros(3))
check("gauge_full_shape", gauged_shape, shape)
hs = [diagonal_part(hh) for hh in hess]
hss = [diagonal_part(hh) for hh in ht]
shape_self = -hinv*sum((lam[aa]*(tangent.T*hs[aa]*tangent)
                       for aa in range(2)), s.zeros(3))
shape_coupling = shape-shape_self
shape_self_tilde = -hinv*sum((lambda_tilde[aa]*(tangent.T*hss[aa]*tangent)
                             for aa in range(2)), s.zeros(3))
shape_coupling_tilde = gauged_shape-shape_self_tilde
shift = -hinv*sum((lambda_tilde[aa]*(tangent.T*diagonal_part(gauge_term[aa])*tangent)
                  for aa in range(2)), s.zeros(3))
check("channel_self_gauge_shift", shape_self_tilde, shape_self+shift)
check("channel_coupling_gauge_shift", shape_coupling_tilde, shape_coupling-shift)
require("channel_gauge_shift_nonzero", shift != s.zeros(3))
check("channel_fingerprint_gauge_matrix", tc*shape_coupling_tilde+ts*shape_self_tilde,
      tc*shape_coupling+ts*shape_self+(ts-tc)*shift)
check("channel_trace_shift_cancellation", s.trace(shape_self_tilde+shape_coupling_tilde), s.trace(shape))

# Test the nonconstant normal-frame connection law directly from the gauged
# Hessians, independent of differentiating a claimed transformation rule.
for ii in range(3):
    tangent_i = tangent[:, ii]
    omega = -s.Matrix(2, 2, lambda aa, bb: (tangent_i.T*hess[aa]*w[:, bb])[0])
    omega_tilde = -s.Matrix(2, 2, lambda aa, bb: (tangent_i.T*ht[aa]*wt[:, bb])[0])
    derivative_a = sum((tangent_i[rr]*gauge.diff(coords[rr]).subs(origin)
                        for rr in range(5)), s.zeros(2))
    check(f"normal_connection_gauge_{ii}", omega_tilde,
          a0*omega*a0.inv()-derivative_a*a0.inv())

# Polarization in tangent dimension two with entirely symbolic, noncommuting
# self-adjoint operators and arbitrary nonsingular tangent metric.
hp, hq, hr = s.symbols("hp hq hr", real=True)
a, b, cc, dd, ee, ff = s.symbols("a b c d e f", real=True)
hm = s.Matrix([[hp, hq], [hq, hr]])
aa = hm.inv()*s.Matrix([[a, b], [b, cc]])
bb = hm.inv()*s.Matrix([[dd, ee], [ee, ff]])
mixed = s.trace(aa)*s.trace(bb)-s.trace(aa*bb)
check("general_second_polarization", (l0*aa+l1*bb).det(),
      l0*l0*aa.det()+l0*l1*mixed+l1*l1*bb.det())

# Nontrivial normal curvature from a cross-coupled two-surface in flat 1+3.
# The slopes vanish at the origin, making omega vanish there while d omega
# remains nonzero.  Thus pointwise omega=0 cannot be mistaken for flatness.
u, v = s.symbols("u v", real=True)
metric4 = s.diag(-1, 1, 1, 1)
f2 = (2*u*u+2*u*v+3*v*v)/2 + u*u*v
g2 = (u*u+4*u*v-v*v)/2 + u*v*v
emb2 = s.Matrix([f2, g2, u, v])
tan2 = emb2.jacobian((u, v))
der2 = s.Matrix([[1, 0, -s.diff(f2, u), -s.diff(f2, v)],
                 [0, 1, -s.diff(g2, u), -s.diff(g2, v)]])
j2, c2, w2, p2 = geometry(metric4, der2)
omega2 = [(der2*w2.diff(var)).applyfunc(simp) for var in (u, v)]
zero2 = {u: 0, v: 0}
for ii in range(2):
    check(f"curved_normal_connection_vanishes_at_point_{ii}", omega2[ii].subs(zero2))
curvature2 = (omega2[1].diff(u)-omega2[0].diff(v)
              +omega2[0]*omega2[1]-omega2[1]*omega2[0]).subs(zero2).applyfunc(simp)
require("normal_curvature_nonzero", curvature2 != s.zeros(2))
w20 = w2.subs(zero2)
normal_metric2 = w20.T*metric4*w20
check("normal_curvature_metric_skew", normal_metric2*curvature2 + curvature2.T*normal_metric2)
second2 = [[emb2.diff(va, vb).subs(zero2) for vb in (u, v)] for va in (u, v)]
shape2 = [s.Matrix(2, 2, lambda ii, jj: (second2[ii][jj].T*metric4*w20[:, nn])[0])
          for nn in range(2)]
comm = shape2[0]*shape2[1]-shape2[1]*shape2[0]
check("normal_curvature_Ricci_identity", (normal_metric2*curvature2)[1, 0], comm[1, 0])
for ii, var in enumerate((u, v)):
    check(f"normal_connection_metric_compatibility_{ii}",
          c2.inv().diff(var), omega2[ii].T*c2.inv()+c2.inv()*omega2[ii])
witnesses["normal_curvature_at_origin"] = str(curvature2)
witnesses["noncommuting_shapes_at_origin"] = [str(ss) for ss in shape2]

# A genuinely varying, off-diagonal ambient metric: the covariant Hessian is
# necessary.  Compare with differentiating the embedding and adding the
# independently calculated Christoffel acceleration.
coord4 = (t, wcoord, x, y)
f3 = x/s.Integer(4)+y/s.Integer(5)+(x*x+2*x*y+3*y*y)/2
g3 = x/s.Integer(7)-y/s.Integer(8)+(2*x*x-2*x*y+y*y)/2
emb3 = s.Matrix([f3, g3, x, y])
rel3 = s.Matrix([t-f3, wcoord-g3])
met3 = s.diag(-1-2*x, 1+y, 1+t, 1+2*wcoord)
met3[0, 2] = met3[2, 0] = x+y
met3[1, 3] = met3[3, 1] = t-x
zero3 = dict.fromkeys(coord4, 0)
met30 = met3.subs(zero3)
inv30 = met30.inv()
d30 = rel3.jacobian(coord4).subs(zero3)
tan30 = emb3.jacobian((x, y)).subs(zero3)
j30, c30, w30, p30 = geometry(met30, d30)
gamma30 = [[[sum(inv30[kk, ll]*(s.diff(met3[ll, jj], coord4[ii])
               +s.diff(met3[ll, ii], coord4[jj])-s.diff(met3[ii, jj], coord4[ll]))
               for ll in range(4)).subs(zero3)/2
             for jj in range(4)] for ii in range(4)] for kk in range(4)]
raw_hess30 = [s.hessian(a, coord4).subs(zero3) for a in rel3]
cov_hess30 = [s.Matrix(4, 4, lambda ii, jj: raw_hess30[aa][ii, jj]
                       -sum(gamma30[kk][ii][jj]*d30[aa, kk] for kk in range(4)))
              for aa in range(2)]
for ii in range(2):
    for jj in range(ii, 2):
        accel = s.Matrix([sum(gamma30[kk][rr][ss]*tan30[rr, ii]*tan30[ss, jj]
                              for rr in range(4) for ss in range(4)) for kk in range(4)])
        direct = (s.eye(4)-p30)*(emb3.diff((x, y)[ii], (x, y)[jj]).subs(zero3)+accel)
        implied = -w30*s.Matrix([(tan30[:, ii].T*hh*tan30[:, jj])[0] for hh in cov_hess30])
        check(f"varying_metric_embedded_B_{ii}{jj}", implied, direct)
require("varying_metric_covariant_hessian_required", cov_hess30 != raw_hess30)

# All codimension-two causal strata, including the zero and null mean cases.
ell, nu = s.symbols("ell nu", real=True)
normal_metric_null = s.Matrix([[0, -1], [-1, 0]])
alpha = s.Matrix([[ell, nu]])
area_dual = s.Matrix([-nu, -ell])
area_preserving = s.Matrix([-nu, ell])
check("area_dual_covector", area_dual.T*normal_metric_null, alpha)
check("area_dual_norm", (area_dual.T*normal_metric_null*area_dual)[0], -2*ell*nu)
check("area_preserving_kernel", alpha*area_preserving)
check("area_preserving_norm", (area_preserving.T*normal_metric_null*area_preserving)[0], 2*ell*nu)
for name, el, nv, causal in [
    ("future_trapped", -2, -3, -1),
    ("past_trapped", 2, 3, -1),
    ("opposite_expansions", 2, -3, 1),
    ("future_marginal", 0, -3, 0),
    ("past_marginal", 2, 0, 0),
]:
    require(f"causal_stratum_{name}", s.sign((-2*ell*nu).subs({ell: el, nu: nv})) == causal)
check("zero_mean_dual_zero", area_dual.subs({ell: 0, nu: 0}), s.zeros(2, 1))
check("zero_mean_every_normal_area_preserving", alpha.subs({ell: 0, nu: 0})*s.Matrix([l0, l1]))

# Forward native construction from an actual family of implicit cuts.  The
# covariant Hessians come from the metric and phi=(T-T0,r-r0); neither D_area
# nor rho is supplied to their construction.
fwd_T, fwd_T0, fwd_r0, fwd_theta, fwd_phi = s.symbols(
    "T_native T0_native r0_native theta_native phi_native", real=True)
fwd_r = s.symbols("r_native", positive=True)
fwd_coords = (fwd_T, fwd_r, fwd_theta, fwd_phi)
fwd_R = s.Function("R")(fwd_T, fwd_r)
fwd_xi = s.Function("xi")(fwd_T, fwd_r)
fwd_metric = s.diag(-1, 1, fwd_R**2, fwd_R**2*s.sin(fwd_theta)**2)
fwd_inverse = fwd_metric.inv()
fwd_relation = s.Matrix([fwd_T-fwd_T0, fwd_r-fwd_r0])
fwd_d = fwd_relation.jacobian(fwd_coords)
fwd_gamma = [s.Matrix(4, 4, lambda ii, jj: simp(sum(
    fwd_inverse[kk, ll]*(s.diff(fwd_metric[ll, jj], fwd_coords[ii])
                        +s.diff(fwd_metric[ll, ii], fwd_coords[jj])
                        -s.diff(fwd_metric[ii, jj], fwd_coords[ll]))
    for ll in range(4))/2)) for kk in range(4)]
fwd_hessians = [s.hessian(fwd_relation[aa], fwd_coords)
                -sum((fwd_d[aa, kk]*fwd_gamma[kk] for kk in range(4)), s.zeros(4))
                for aa in range(2)]
fwd_j, fwd_c, fwd_w, fwd_p = geometry(fwd_metric, fwd_d)
fwd_tangent = s.eye(4)[:, 2:4]
fwd_h = fwd_tangent.T*fwd_metric*fwd_tangent
fwd_e = [fwd_tangent.T*hh*fwd_tangent for hh in fwd_hessians]
fwd_t = s.Matrix([simp(s.trace(fwd_h.inv()*ee)) for ee in fwd_e])
check("implicit_cut_normal_gram", fwd_c, s.diag(-1, 1))
check("implicit_cut_normal_frame", fwd_w, s.eye(4)[:, 0:2])
check("implicit_cut_hessian_time", fwd_e[0], -s.diff(fwd_R, fwd_T)*fwd_h/fwd_R)
check("implicit_cut_hessian_radial", fwd_e[1], s.diff(fwd_R, fwd_r)*fwd_h/fwd_R)
check("implicit_cut_traced_hessian", fwd_t,
      s.Matrix([-2*s.diff(fwd_R, fwd_T)/fwd_R, 2*s.diff(fwd_R, fwd_r)/fwd_R]))

fwd_K = s.Matrix([1, fwd_xi, 0, 0])
fwd_J = s.Matrix([-1, 1, 0, 0])/(1+fwd_xi)
fwd_L = (1+fwd_xi)*s.Matrix([1, 1, 0, 0])/2
fwd_N = s.Matrix([1, -1, 0, 0])/(1+fwd_xi)
fwd_F = 1-fwd_xi**2
fwd_pair = lambda left, right: simp((left.T*fwd_metric*right)[0])
fwd_expansion = lambda normal_vector: simp((fwd_t.T*fwd_c.inv()*fwd_d*normal_vector)[0])
fwd_D = fwd_expansion(fwd_K)
fwd_rho = fwd_expansion(fwd_J)
fwd_theta_L = fwd_expansion(fwd_L)
fwd_theta_N = fwd_expansion(fwd_N)
check("implicit_cut_rest_normal_norm", fwd_pair(fwd_K, fwd_K), -fwd_F)
check("implicit_cut_rest_null_partner", fwd_pair(fwd_K, fwd_J), 1)
check("implicit_cut_outgoing_null", fwd_pair(fwd_L, fwd_L))
check("implicit_cut_inward_null", fwd_pair(fwd_N, fwd_N))
check("implicit_cut_null_pair", fwd_pair(fwd_L, fwd_N), -1)
check("implicit_cut_outgoing_normal_bridge", fwd_L, fwd_K+fwd_F*fwd_J/2)
check("implicit_cut_forward_area_K", fwd_D,
      2*(s.diff(fwd_R, fwd_T)+fwd_xi*s.diff(fwd_R, fwd_r))/fwd_R)
check("implicit_cut_forward_area_J", fwd_rho,
      2*(s.diff(fwd_R, fwd_r)-s.diff(fwd_R, fwd_T))/(fwd_R*(1+fwd_xi)))
check("implicit_cut_forward_outgoing", fwd_theta_L,
      (1+fwd_xi)*(s.diff(fwd_R, fwd_T)+s.diff(fwd_R, fwd_r))/fwd_R)
check("implicit_cut_forward_inward", fwd_theta_N,
      2*(s.diff(fwd_R, fwd_T)-s.diff(fwd_R, fwd_r))/(fwd_R*(1+fwd_xi)))
check("implicit_cut_area_bridge", fwd_theta_L, fwd_D+fwd_rho*fwd_F/2)
fwd_shape_K = -fwd_h.inv()*sum(((fwd_c.inv()*fwd_d*fwd_K)[aa]*fwd_e[aa]
                               for aa in range(2)), s.zeros(2))
check("implicit_cut_full_rest_shape", fwd_shape_K,
      -(s.diff(fwd_R, fwd_T)+fwd_xi*s.diff(fwd_R, fwd_r))*s.eye(2)/fwd_R)

# The resolved constitutive law supplies xi(T), while R(T,r) remains the
# separately chosen area/metric profile.  The latter is differentiated forward.
fwd_k = s.symbols("k_native", positive=True)
fwd_witness_R = fwd_r-fwd_T-fwd_k*fwd_T**2/2
fwd_witness_xi = 1+fwd_k*fwd_T
fwd_resolved_response = s.Matrix([[1+fwd_k**2, fwd_k], [fwd_k, 1]])
fwd_pushforward = s.Matrix([[fwd_T, fwd_witness_xi], [0, 1]])
fwd_original_effort = s.Matrix([0, 1])
check("implicit_cut_resolved_response_positive_determinant", fwd_resolved_response.det(), 1)
check("implicit_cut_resolved_load_generates_clock_path",
      fwd_resolved_response*fwd_pushforward.T*fwd_original_effort,
      s.Matrix([s.diff(fwd_witness_xi, fwd_T), 1]))


def fwd_specialize(expression):
    substituted = expression.subs(
        {fwd_R: fwd_witness_R, fwd_xi: fwd_witness_xi}, simultaneous=True).doit()
    return substituted.applyfunc(simp) if isinstance(substituted, s.MatrixBase) else simp(substituted)


fwd_witness_D = fwd_specialize(fwd_D)
fwd_witness_rho = fwd_specialize(fwd_rho)
fwd_witness_theta_L = fwd_specialize(fwd_theta_L)
fwd_witness_theta_N = fwd_specialize(fwd_theta_N)
fwd_witness_N = fwd_specialize(fwd_N)
fwd_witness_L = fwd_specialize(fwd_L)
fwd_witness_F = 1-fwd_witness_xi**2
check("implicit_cut_witness_area_K", fwd_witness_D)
check("implicit_cut_witness_area_J", fwd_witness_rho, 2/fwd_witness_R)
check("implicit_cut_witness_inward", fwd_witness_theta_N, -2/fwd_witness_R)
check("implicit_cut_witness_outgoing", fwd_witness_theta_L, fwd_witness_F/fwd_witness_R)
check("implicit_cut_witness_initial_marginality", fwd_witness_theta_L.subs(fwd_T, 0))
fwd_outer = simp(sum(fwd_witness_N[ii]*s.diff(fwd_witness_theta_L, fwd_coords[ii])
                     for ii in range(4)))
check("implicit_cut_witness_outer_derivative", fwd_outer.subs(fwd_T, 0), -fwd_k/fwd_r)
check("implicit_cut_witness_regular_crossing", s.diff(fwd_witness_F, fwd_T).subs(fwd_T, 0), -2*fwd_k)

# The tensor-calculus extension independently supplies the full Ricci tensor
# of precisely the metric used by the implicit-cut calculation.
fwd_witness_metric = fwd_specialize(fwd_metric)
fwd_curvature = CurvatureNumerators(fwd_witness_metric, fwd_coords)
fwd_ricci = fwd_curvature.ricci.applyfunc(s.trigsimp)
fwd_angular_ricci = fwd_witness_xi**2-fwd_k*fwd_witness_R
fwd_expected_ricci = s.diag(2*fwd_k/fwd_witness_R, 0, fwd_angular_ricci,
                           fwd_angular_ricci*s.sin(fwd_theta)**2)
check("implicit_cut_tensor_full_Ricci", fwd_ricci, fwd_expected_ricci)
check("implicit_cut_tensor_null_focusing", (fwd_witness_L.T*fwd_ricci*fwd_witness_L)[0],
      fwd_k*(1+fwd_witness_xi)**2/(2*fwd_witness_R))
check("implicit_cut_tensor_scalar_curvature", s.trace(fwd_witness_metric.inv()*fwd_ricci),
      2*fwd_witness_xi**2/fwd_witness_R**2-4*fwd_k/fwd_witness_R)
witnesses["native_implicit_cut"] = {
    "cut": "phi=(T-T0,r-r0)",
    "area_K_general": str(fwd_D),
    "area_profile_chosen": str(fwd_witness_R),
    "xi_from_constant_resolved_response": str(fwd_witness_xi),
    "area_K_witness": str(fwd_witness_D),
    "rho_witness": str(fwd_witness_rho),
    "theta_L_witness": str(fwd_witness_theta_L),
    "outer_derivative_at_T0": str(simp(fwd_outer.subs(fwd_T, 0))),
    "Ricci_diagonal": [str(simp(fwd_ricci[ii, ii])) for ii in range(4)],
}

# A trace-free, nonzero second fundamental form can have zero mean vector.
zero_mean_b0 = s.diag(1, -1)
zero_mean_b1 = s.Matrix([[0, 2], [2, 0]])
check("zero_mean_nontrivial_B_time_trace", s.trace(zero_mean_b0))
check("zero_mean_nontrivial_B_space_trace", s.trace(zero_mean_b1))
require("zero_mean_does_not_force_totally_geodesic", zero_mean_b0 != s.zeros(2) and zero_mean_b1 != s.zeros(2))

report = {
    "date": "2026-09-05",
    "scope": "Conditional nondegenerate implicit submanifold geometry; native cut and transport selection remain inputs.",
    "checks_passed": len(checks),
    "checks": checks,
    "witnesses": witnesses,
    "source_provenance": {
        "GFL:thm:mean_curvature_decomp": "ace98fbfd85a92cca4da4f348c410c81fe38a0409d57840915d98fa3359f9033",
        "DBP:thm:three_channel_kg_ext": "dcb8e175e3c43f9b708afda98e8cbf02c55f910e2d67e28b300edf676ab7c2fd",
    },
}
target = Path(__file__).resolve().parent / "docs" / "cella-normal-geometry-checks.json"
target.write_text(json.dumps(report, indent=2)+"\n")
print(json.dumps({"checks_passed": len(checks), "report": str(target), "witnesses": witnesses}, indent=2))
