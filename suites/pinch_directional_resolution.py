#!/usr/bin/env python3
"""Exact directional continuation of PINCH-1's finite normal limits.

This is a proposed extension of the state space, not a selected dynamical law.
See docs/results/2026-09-05/2026-09-05-pinch-directional-resolution.md for the proof and scope.
"""
import json
import sympy as s

checks = {}


def check(name, actual, expected=0):
    difference = actual - expected
    entries = list(difference) if isinstance(difference, s.MatrixBase) else [difference]
    for entry in entries:
        residual = s.factor(s.cancel(s.expand(entry)))
        if residual != 0:
            residual = s.simplify(residual)
        if residual != 0:
            raise AssertionError((name, residual))
    checks[name] = {"verified": True}


a, b, gamma, d, xi = s.symbols("a b gamma d xi", real=True)
p = 1 + a**2
eps = 1 - gamma**2
delta = eps + a**2 + b**2 - 2*a*b*gamma
check("directional_identity", delta - (b-a*gamma)**2, p*eps)

F = (1-xi**2)/p
resolved = s.Matrix([[-1,a,xi], [a,1,0], [xi,0,F]])
check("resolved_determinant", resolved.det(), -1)
q = s.diag(1,1,-1)
c = s.Matrix([0,0,1])
h = s.Matrix([s.sqrt(p),0,-a])
v = s.Matrix([a*xi/s.sqrt(p),1/s.sqrt(p),-xi])
frame = s.Matrix.hstack(c,h,v)
check("resolved_oriented_volume", frame.det(), 1)
check("resolved_ambient_gram", frame.T*q*frame, resolved)
check("normalized_difference_ruler_pairing", (h.T*q*v)[0])
check("normalized_difference_clock_pairing", (c.T*q*v)[0], xi)
check("normalized_difference_norm", (v.T*q*v)[0], F)

nu = s.Matrix([a/s.sqrt(p),xi/s.sqrt(p),-1])
nu_coeff = s.Matrix([-F,a*F,xi])
check("resolved_dual_normal_coefficients", resolved*nu_coeff, s.Matrix([1,0,0]))
check("resolved_dual_normal_vector", frame*nu_coeff, nu)
check("resolved_dual_normal_pairings", frame.T*q*nu, s.Matrix([1,0,0]))
check("resolved_dual_normal_norm", (nu.T*q*nu)[0], -F)
check("tilt_identity", (a*a+xi*xi)/(1-xi*xi), 1/F-1)

for sigma in (1,-1):
    gam = sigma*s.sqrt(1-F*d*d)
    bb = a*gam+xi*d
    g = gam*h+d*v
    original_frame = s.Matrix.hstack(c,h,g)
    gram = s.Matrix([[-1,a,bb], [a,1,gam], [bb,gam,1]])
    prefix = f"branch_{sigma}"
    check(prefix+"_original_gram", original_frame.T*q*original_frame, gram)
    check(prefix+"_original_volume", original_frame.det(), d)
    check(prefix+"_original_determinant", gram.det(), -d*d)
    check(prefix+"_original_lapse", 1-gam*gam, F*d*d)
    check(prefix+"_direction_coordinate", bb-a*gam, xi*d)
    check(prefix+"_normalized_difference", g-gam*h, d*v)
    check(prefix+"_pinch_gamma", gam.subs(d,0), sigma)
    check(prefix+"_pinch_depth", bb.subs(d,0), sigma*a)
    check(prefix+"_pinch_merged_rulers", g.subs(d,0), sigma*h)

center = (-c+a*h)/p
check("pinch1_transverse_limit", nu.subs(xi,0), center)
check("pinch1_bisector", (nu.subs(xi,1)+nu.subs(xi,-1))/2, center)
check("pinch1_center_norm", (center.T*q*center)[0], -1/p)
for tau in (1,-1):
    check(f"null_direction_{tau}", F.subs(xi,tau))
    check(f"regular_horizon_gradient_{tau}", s.diff(F,xi).subs(xi,tau), -2*tau/p)

# For fixed finite a, every finite slope xi has an exact real chart for small d.
# In particular, these checks include the hyperbolic side, not only rest limits.
check("compact_center", F.subs({a:0,xi:0}), 1)
check("compact_intermediate", F.subs({a:0,xi:s.sqrt(s.Rational(1,3))}), s.Rational(2,3))
check("hyperbolic_example", F.subs({a:0,xi:2}), -3)

# The stationary radial completion has a finite surface-gravity coefficient in
# these variables, conditional on smooth A,B and radial state fields.
A,B,ar,xir = s.symbols("A B a_r xi_r", real=True)
dFdr = s.diff(F,a)*ar+s.diff(F,xi)*xir
for tau in (1,-1):
    check(f"stationary_surface_gravity_{tau}", (A*dFdr/(2*B)).subs(xi,tau),
          -tau*A*xir/(B*p))

print(json.dumps({
    "scope": "Exact finite-direction resolution of the native pinch. State-space extension is proposed; direction selection and area-transport dynamics remain open.",
    "verified_count": len(checks),
    "checks": checks,
    "results": {
        "retained_direction": "xi = (b-a gamma)/d, d = oriented det[c,h,g], d^2 = delta",
        "resolved_lapse": "F = (1-xi^2)/(1+a^2)",
        "replacement_ruler": "v = (g-gamma h)/d",
        "resolved_gram": "[[-1,a,xi],[a,1,0],[xi,0,F]], determinant -1",
        "dual_normal": "nu = -F c + a F h + xi v",
        "pinch_fiber": "finite xi is retained at d=0",
        "horizons": "xi=+1 or xi=-1, regular including d=0",
        "next_native_target": "derive the retained direction and area transport from the native dynamics"
    }
}, indent=2))
