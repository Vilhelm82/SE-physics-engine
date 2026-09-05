"""Exact arithmetic receipts for the foundation audit, 2026-09-05.

These check counterexamples derived in the audit, independently of the repository
runners. They do not validate a physical model or establish literature novelty.
Run: python3 audit-counterexamples.py
"""
import json
import sympy as s

out = {}

# A finite, realizable compact-sector path has both zero lapse and rank loss.
t = s.symbols('t', positive=True)
g = 1-t**4
G = s.Matrix([[-1,t,0],[t,1,g],[0,g,1]])
eta = s.cancel(t**2/(1-g**2))
N2 = s.cancel(1/(1+eta))
assert s.expand(G.det()-t**2*(t**6-2*t**2-1)) == 0
assert s.limit(eta,t,0,dir='+') == s.oo
assert s.limit(N2,t,0,dir='+') == 0
assert G.subs(t,0).rank() == 2
out['joint_lapse_and_rank_loss'] = {
    'domain': '0 < t < 1: compact ruler Gram positive definite; full inertia (2,1)',
    'determinant': str(s.factor(G.det())), 'eta': str(eta), 'lapse_squared': str(N2),
    'limit_lapse_squared': 0, 'endpoint_rank': 2}

# Pair-plane type is not a count of coordinate boosts in an orthogonal basis.
G3 = s.Matrix([[1,-2,-2],[-2,1,-2],[-2,-2,1]])
assert G3.eigenvals() == {3:2,-3:1}
assert all(G3.extract(p,p).det() == -3 for p in [(0,1),(0,2),(1,2)])
out['three_hyperbolic_pair_spans'] = {'eigenvalues': ['3','3','-3'], 'pair_determinants': [-3,-3,-3]}

# Trine frames: SO(3) action remains free and mechanical curvature finite at V=0.
axes = [s.Matrix([1,0,0]), s.Matrix([-s.Rational(1,2),s.sqrt(3)/2,0]),
        s.Matrix([-s.Rational(1,2),-s.sqrt(3)/2,0])]
ez = s.Matrix([0,0,1])
X = [ez/s.sqrt(3) for _ in axes]
y = [1/s.sqrt(2),-1/s.sqrt(2),0]
Y = [yi*ez.cross(ai) for yi,ai in zip(y,axes)]
sumv = lambda vs: sum(vs,s.zeros(3,1))
inertia = 3*s.eye(3)-sum((ai*ai.T for ai in axes),s.zeros(3))
assert inertia == s.diag(s.Rational(3,2),s.Rational(3,2),3)
assert sumv(ai.cross(xi) for ai,xi in zip(axes,X)) == s.zeros(3,1)
assert sumv(ai.cross(yi) for ai,yi in zip(axes,Y)) == s.zeros(3,1)
assert s.simplify(sum(xi.dot(xi) for xi in X)) == 1
assert s.simplify(sum(yi.dot(yi) for yi in Y)) == 1
assert sum(xi.dot(yi) for xi,yi in zip(X,Y)) == 0
F = s.simplify(2*inertia.inv()*sumv(xi.cross(yi) for xi,yi in zip(X,Y)))
Fnorm2 = s.simplify(F.dot(F))
vert2 = s.simplify((F.T*inertia*F)[0])
Ktotal = s.simplify(sum(xi.dot(xi)*yi.dot(yi)-xi.dot(yi)**2 for xi,yi in zip(X,Y)))
Kbase = s.simplify(Ktotal+s.Rational(3,4)*vert2)
assert Fnorm2 == s.Rational(8,9)
assert Kbase == s.Rational(4,3)
out['regular_coplanar_quotient'] = {'inertia': str(inertia), 'F_norm_squared': str(Fnorm2),
    'product_sectional_curvature': str(Ktotal), 'quotient_sectional_curvature': str(Kbase)}

# Swapping central algebra blocks need not be a fourth Clifford vector.
J = s.Matrix([[0,1],[-1,0]])
H = s.Matrix([[0,1],[1,0]])
K = s.diag(1,-1)
GC,GH,GG = s.diag(J,-J),s.diag(H,H),s.diag(K,K)
omega = GC*GH*GG
E = s.Matrix.vstack(s.Matrix.hstack(s.zeros(2),s.eye(2)),s.Matrix.hstack(s.eye(2),s.zeros(2)))
assert E*omega+omega*E == s.zeros(4)
assert E*GH-GH*E == s.zeros(4)
assert E*GH+GH*E != s.zeros(4)
out['block_swap_is_not_necessarily_a_vector'] = {'swaps_blocks': True, 'anticommutes_with_all_vectors': False}

# A real symplectic bilinear does not define a positive expectation rule.
B = s.diag(J,J)
A = s.Matrix.vstack(s.Matrix.hstack(s.zeros(2),J),s.Matrix.hstack(-J,s.zeros(2)))
u = s.Matrix(s.symbols('u0:4',real=True))
assert B.inv()*A.T*B == A
assert s.expand((u.T*B*u)[0]) == 0
assert s.expand((u.T*B*A*u)[0]) == 0
out['symplectic_adjoint_without_probabilities'] = {'A_is_B_self_adjoint': True,
    'u_transpose_B_u': 0, 'u_transpose_B_A_u': 0}

# Both inverse branches of cos(2 phi)=u: adding pi/2 does not give the second.
u0 = s.Rational(1,3)
bad = s.acos(u0)/2+s.pi/2
good = s.pi-s.acos(u0)/2
assert s.simplify(s.cos(2*bad)) == -u0
assert s.simplify(s.cos(2*good)) == u0
out['cover_inverse_branch'] = {'chosen_u': str(u0), 'old_branch_1_image': str(-u0), 'correct_branch_1_image': str(u0)}

print(json.dumps(out,indent=2))
