#!/usr/bin/env python3
"""Exact checks supporting the rest-frame / trapping equivalence theorem.

The proof and hypotheses are in docs/results/2026-09-05/2026-09-05-rest-frame-trapping-equivalence.md.
The native scalar F is preserved. No area-transport law is imposed silently.
"""
import json
import sympy as s

checks = {}


def check(name, expression, expected=0):
    residual = s.factor(s.cancel(s.simplify(expression-expected)))
    if residual != 0:
        raise AssertionError((name, residual))
    checks[name] = {"verified": True}


def matrix_check(name, actual, expected):
    for entry in actual-expected:
        if s.simplify(entry) != 0:
            raise AssertionError((name, entry))
    checks[name] = {"verified": True}


# General normal-plane theorem; expansion is linear over smooth coefficients
# for vectors normal to the surface.
F, D = s.symbols("F D", real=True)
rho, boost = s.symbols("rho boost", positive=True)
normal_metric = s.Matrix([[-F, 1], [1, 0]])
K = s.Matrix([1, 0])
J = s.Matrix([0, 1])
L = K+F*J/2
N = -J
expansion = lambda V: D*V[0]+rho*V[1]
pair = lambda V,W: (V.T*normal_metric*W)[0]
check("native_normal_norm", pair(K,K), -F)
check("outgoing_null", pair(L,L))
check("ingoing_null", pair(N,N))
check("null_pair_normalization", pair(L,N), -1)
check("expansion_bridge", expansion(L), D+rho*F/2)
check("inward_expansion", expansion(N), -rho)

Ftrap = F+2*D/rho
Karea = K-D*J/rho
check("area_normal_expansion", expansion(Karea))
check("area_normal_pairing", pair(Karea,J), 1)
check("area_normal_norm", pair(Karea,Karea), -Ftrap)
Z = rho*L+expansion(L)*N
matrix_check("dual_mean_curvature_vector", Z, rho*K-D*J)
check("dual_mean_curvature_norm", pair(Z,Z), -rho**2*Ftrap)
check("positive_boost_covariance",
      boost**2*F+2*(boost*D)/(rho/boost), boost**2*Ftrap)
Xi = s.symbols("Xi", real=True)
check("regular_equivalence_factor", Ftrap.subs(D,F*Xi), F*(1+2*Xi/rho))

# Arbitrary transverse metric and angular shift: verify the area formula by
# tracing the full spacetime Lie derivative, rather than assuming a radial ansatz.
v,r,x,y = s.symbols("v r x y", real=True)
coords = (v,r,x,y)
A,B,Ff,U1,U2 = (s.Function(n)(*coords) for n in ("A","B","F","U1","U2"))
h11,h12,h22 = (s.Function(n)(*coords) for n in ("h11","h12","h22"))
h = s.Matrix([[h11,h12],[h12,h22]])
hinv = h.inv()
U = s.Matrix([U1,U2])
metric = s.zeros(4)
metric[0,0] = -A*A*Ff+(U.T*h*U)[0]
metric[0,1] = metric[1,0] = A*B
for i in range(2):
    metric[0,i+2] = metric[i+2,0] = -(h*U)[i]
    for j in range(2):
        metric[i+2,j+2] = h[i,j]
Kf = s.Matrix([1/A,0,U1/A,U2/A])
Jf = s.Matrix([0,1/B,0,0])


def screen_expansion(V):
    result = 0
    for i in range(2):
        for j in range(2):
            a,b = i+2,j+2
            lie = sum(V[k]*s.diff(metric[a,b],coords[k]) +
                      metric[k,b]*s.diff(V[k],coords[a]) +
                      metric[a,k]*s.diff(V[k],coords[b])
                      for k in range(4))
            result += hinv[i,j]*lie/2
    return s.factor(result)


det_h = h.det()
Dexpected = (s.diff(det_h,v)/(2*det_h) +
             U1*s.diff(det_h,x)/(2*det_h) +
             U2*s.diff(det_h,y)/(2*det_h) +
             s.diff(U1,x)+s.diff(U2,y))/A
rho_expected = s.diff(det_h,r)/(2*B*det_h)
check("general_metric_native_normal", (Kf.T*metric*Kf)[0], -Ff)
check("general_metric_normal_pairing", (Kf.T*metric*Jf)[0], 1)
check("full_metric_area_drift", screen_expansion(Kf), Dexpected)
check("full_metric_inward_contraction", screen_expansion(Jf), rho_expected)
check("full_metric_expansion_bridge", screen_expansion(Kf+Ff*Jf/2),
      Dexpected+Ff*rho_expected/2)

# Transverse conformal rescaling preserves the normal metric while changing its
# first derivative area data. Determinant h -> exp(4 chi) determinant h.
chi = s.Function("chi")(*coords)
newdet = s.exp(4*chi)*det_h
newD = (s.diff(newdet,v)/(2*newdet) +
        U1*s.diff(newdet,x)/(2*newdet) +
        U2*s.diff(newdet,y)/(2*newdet) +
        s.diff(U1,x)+s.diff(U2,y))/A
newrho = s.diff(newdet,r)/(2*B*newdet)
Kchi = (s.diff(chi,v)+U1*s.diff(chi,x)+U2*s.diff(chi,y))/A
check("arbitrary_area_first_derivative", newD-Dexpected, 2*Kchi)
check("conformal_inward_contraction", newrho-rho_expected, 2*s.diff(chi,r)/B)

# An actual native state family admits a flat-spacetime soldering.
gamma, radius, theta = s.symbols("gamma radius theta", real=True)
delta = 2-gamma**2
G = s.Matrix([[-1,-1,0,0],[-1,1,gamma,0],
              [0,gamma,1,0],[0,0,0,1]])
c = s.Matrix([1,0,0,0])
er = s.Matrix([-1,1,-gamma,0])/s.sqrt(delta)
ep = s.Matrix([0,0,1,0])
e4 = s.Matrix([0,0,0,1])
flat_solder = s.Matrix.hstack(c,er,radius*ep,radius*s.sin(theta)*e4)
flat_metric = s.diag(-1,1,radius**2,radius**2*s.sin(theta)**2)
matrix_check("native_flat_soldering", flat_solder.T*G*flat_solder, flat_metric)
check("native_flat_soldering_determinant", G.det(), -delta)
Fnative = (1-gamma**2)/delta
native_k = -G.inv()*c
matrix_check("native_flat_soldering_normal",
             flat_solder*s.Matrix([1,1/s.sqrt(delta),0,0]), native_k)
check("flat_native_horizon", Fnative.subs(gamma,1))
check("flat_native_horizon_full_rank", delta.subs(gamma,1), 1)

w = s.symbols("w", positive=True)
gm = s.diag(-1,1)
Km = s.Matrix([1,w])
Jm = s.Matrix([-1,1])/(1+w)
Lm = Km+(1-w*w)*Jm/2
Nm = -Jm
Dm = 2*w/radius
rhom = 2/(radius*(1+w))
check("flat_outgoing_expansion", 2*Lm[1]/radius, (1+w)/radius)
check("flat_ingoing_expansion", 2*Nm[1]/radius, -rhom)
check("flat_geometric_trapping_invariant", 1-w*w+2*Dm/rhom, (1+w)**2)
matrix_check("flat_area_rest_vector", Km-Dm*Jm/rhom, s.Matrix([1+w,0]))
check("flat_horizon_still_expanding", ((1+w)/radius).subs(w,1), 2/radius)

# A null-aligned boundary can still have nonzero outgoing expansion.
lam = s.symbols("lambda", real=True)
R = s.exp(lam*v)*r
drift_R = 2*s.diff(R,v)/R
contract_R = 2*s.diff(R,r)/R
check("null_aligned_example_area_drift", drift_R, 2*lam)
check("null_aligned_example_inward_contraction", contract_R, 2/r)
check("null_aligned_example_horizon_expansion",
      (drift_R+F*contract_R/2).subs(F,0), 2*lam)
check("null_aligned_example_trapping_shift",
      F+2*drift_R/contract_R, F+2*lam*r)

# A broad matching class: time-dependent shape and angular drift can preserve
# transverse area. This is a sufficient class, not an assumed native law.
theta,phi = s.symbols("theta phi",real=True)
shape = s.Function("s")(v,r,theta)
Rr = s.Function("R")(r)
Omega = s.Function("Omega")(v,r,theta)
hs = s.diag(Rr**2*s.exp(2*shape),
            Rr**2*s.exp(-2*shape)*s.sin(theta)**2)
check("distorted_screen_area", hs.det(), Rr**4*s.sin(theta)**2)
check("distorted_screen_area_drift",
      s.diff(hs.det(),v)/(2*hs.det())+s.diff(Omega,phi))
check("distorted_screen_radial_expansion",
      s.diff(hs.det(),r)/(2*hs.det()), 2*s.diff(Rr,r)/Rr)

# Native gradient in the two horizon branches, and the outer derivative law.
a,b,gam = s.symbols("a b gamma",real=True)
delta_general = 1-gam**2+a*a+b*b-2*a*b*gam
Fs = (1-gam**2)/delta_general
da,db,dg = s.symbols("da db dgamma",real=True)
dFs = s.diff(Fs,a)*da+s.diff(Fs,b)*db+s.diff(Fs,gam)*dg
for sigma in (1,-1):
    check(f"native_gradient_branch_{sigma}", dFs.subs(gam,sigma),
          -2*sigma*dg/(a-sigma*b)**2)
Lambda,NF = s.symbols("Lambda NF",real=True)
check("outer_derivative_coefficient",
      (rho/2+Xi)*NF, rho*(1+2*Xi/rho)*NF/2)

print(json.dumps({
    "scope": "Exact bridge identities, general transverse area formula, explicit non-equivalence examples, and a sufficient matching class. Native area-transport dynamics remains a theorem target.",
    "verified_count": len(checks),
    "checks": checks,
    "results": {
        "expansion_bridge": "theta_L = D + rho F/2; rho = -theta_N > 0",
        "geometric_trapping_invariant": "F_trap = F + 2 D/rho",
        "area_rest_vector": "K_area = K - (D/rho) J",
        "regular_equivalence_iff": "D = F Xi and 1 + 2 Xi/rho > 0 on the regular boundary",
        "strong_matching_condition": "D = theta_K = 0",
        "null_generator_alignment": "dF proportional to K_flat on F=0",
        "flat_counterexample": "F crosses zero while F_trap = (1+w)^2 > 0",
        "null_aligned_counterexample": "R = exp(lambda v) r gives theta_L|F=0 = 2 lambda"
    }
},indent=2))
