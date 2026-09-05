"""Exact checks supporting the corpus/Seated Root comparison (2026-09-05).

Run with Python 3 and SymPy. These examples check stated mathematical
boundaries, not physical viability or historical novelty.
"""
import json
import sympy as s

out = {}

# A reusable valuation identity in the compact-ruler sector.
a, b, gamma, eps, d = s.symbols('a b gamma eps d', real=True)
G = s.Matrix([[-1,a,b],[a,1,gamma],[b,gamma,1]])
assert s.expand(-G.det() - ((1+b*b)*(1-gamma*gamma)+(a-b*gamma)**2)) == 0
t = s.symbols('t', positive=True)
limits = {}
for p,q in [(4,1),(2,1),(2,2)]:
    N2 = t**p/(5*t**p+9*t**(2*q))
    limits[f'p={p},q={q}'] = str(s.limit(N2,t,0,dir='+'))
assert limits == {'p=4,q=1':'0','p=2,q=1':'1/14','p=2,q=2':'1/5'}
assert s.limit((t**4/(5*t**4+9*t**2))/t**2,t,0) == s.Rational(1,9)
out['lapse_valuation_bridge'] = {'identity':'-det G=(1+b^2)(1-gamma^2)+(a-b gamma)^2',
    'sample_limits_b2_d3':limits}

# The Lorentzian coupling form has the claimed algebraic signature.
Q = 2*s.eye(3)-s.ones(3)
assert Q.eigenvals() == {2:2,-1:1}
out['coupling_form'] = {'eigenvalues':[2,2,-1]}

# Constraint-gradient contributions may cancel in the coupling graph.
J = s.Matrix([[1,1],[1,-1]])
assert J.det() != 0 and (J.T*J)[0,1] == 0
out['lifted_graph_cancellation'] = {'constraints':['x+y','x-y'],
    'spatial_Hessian':str(2*J.T*J),'mixed_entry':0}

# The squared coupling carrier is locally invertible, not globally injective.
A,B,C = s.symbols('A B C', real=True)
carrier = s.Matrix([-B**2/9,-(A-B)**2/9,-(C-B)**2/9]) # a=b=1, q0=3
jet = (2,1,3)
image = carrier.subs(dict(zip((A,B,C),jet)))
assert image == carrier.subs(dict(zip((A,B,C),(-2,-1,-3))))
# At a=b=1 the only active S3 elements preserving this first jet are
# identity and input swap. Neither relates the displayed second jets.
assert (-2,-1,-3) not in [jet,(jet[2],jet[1],jet[0])]
det = s.factor(carrier.jacobian((A,B,C)).det())
assert s.simplify(det-8*B*(A-B)*(C-B)/729) == 0
equal = carrier.subs({A:2,B:1,C:2})
anisotropy = sum((equal[i]-equal[j])**2 for i,j in [(0,1),(1,2),(2,0)])
assert anisotropy == 0 and det.subs({A:2,B:1,C:2}) != 0
out['squared_carrier'] = {'distinct_jets':[jet,(-2,-1,-3)],'same_carrier':str(image),
    'zero_anisotropy_full_rank_det':str(det.subs({A:2,B:1,C:2}))}

# A smooth positive unit creates a pole after leading monomial cancellation.
# Compute Ricci from Christoffels independently of the corpus normal-form code.
x,y,z = s.symbols('x y z', positive=True)
coords = [x,y,z]
metric = s.diag(x*x,1+x,1)
inv = metric.inv()
Christoffel = [[[s.simplify(sum(inv[k,l]*(s.diff(metric[l,j],coords[i])+
    s.diff(metric[l,i],coords[j])-s.diff(metric[i,j],coords[l]))/2 for l in range(3)))
    for j in range(3)] for i in range(3)] for k in range(3)]
Ricci = s.zeros(3)
for i in range(3):
    for j in range(3):
        Ricci[i,j] = s.simplify(sum(s.diff(Christoffel[k][i][j],coords[k])-
            s.diff(Christoffel[k][i][k],coords[j]) + sum(
            Christoffel[k][i][j]*Christoffel[l][k][l]-
            Christoffel[l][i][k]*Christoffel[k][j][l] for l in range(3)) for k in range(3)))
R = s.factor(sum(inv[i,j]*Ricci[i,j] for i in range(3) for j in range(3)))
assert s.simplify(R-(2+3*x)/(2*x**3*(1+x)**2)) == 0
assert s.limit(x**3*R,x,0) == 1
out['corner_unit_after_cancellation'] = {'metric':'diag(x^2,1+x,1)',
    'scalar_curvature':str(R),'limit_x3_R':1,'unit_removed_metric':'diag(x^2,1,1), flat for x>0'}

# Vanishing relative step alone does not control the derivative of its error.
f = s.symbols('f', positive=True)
eta = f*(2+s.sin(1/f**2))
transfer = f*(2+eta) # g(f)=f^2, step delta_f=f*eta
slope = s.simplify(f*s.diff(transfer,f)/transfer)
on_sequence = s.simplify(slope.subs({s.sin(f**(-2)):0,s.cos(f**(-2)):1}))
assert s.limit(on_sequence,f,0) == -s.oo
out['transfer_slope_step_condition'] = {'observable':'g(f)=f^2',
    'relative_step':str(eta),'transfer':str(transfer),
    'slope_on_f_equals_inverse_sqrt_2pi_n':str(on_sequence),
    'slope_subsequence_limit':'-infinity',
    'amplitude_ratio_T_over_2f_limit':str(s.limit(transfer/(2*f),f,0))}
assert out['transfer_slope_step_condition']['amplitude_ratio_T_over_2f_limit'] == '1'

print(json.dumps(out,indent=2))
