"""Fresh Lorentzian derivations for the Seated Root upgrade assessment.

Inputs are a real form of signature (2,1), unit selected poles, an explicit
common-pivot observation protocol, and separately labelled response proposals.
No retired Euclidean runners or Clifford representations are imported.
"""
import json
import sympy as s
results = {}
def record(name, expr, target=0):
    residual=s.simplify(expr-target)
    assert residual==0, (name,residual)
    results[name]={'verified':True,'expression':str(s.simplify(expr))}

# A. The unnormalised dual normal stays regular at the ruler-plane horizon.
a,b,c = s.symbols('a b gamma', real=True)
G=s.Matrix([[-1,a,b],[a,1,c],[b,c,1]])
eps=1-c*c
delta=eps+a*a+b*b-2*a*b*c
ell=s.Matrix([1,0,0])
nu=G.inv()*ell
record('gram_determinant', G.det(), -delta)
for i in range(3): record('dual_normal_'+str(i),(G*nu-ell)[i])
record('normal_norm', (nu.T*G*nu)[0], -eps/delta)
boundary={a:-1,b:0,c:1}
assert G.subs(boundary).det()==-1
assert nu.subs(boundary)==s.Matrix([0,-1,1])
record('regular_null_normal_at_boundary',(nu.T*G*nu)[0].subs(boundary))

# B. A proposed positive kinetic metric constructed from the seated norm,
# equal ruler weights, and elimination of the common azimuthal velocity.
u,v,phi=s.symbols('u v phi',real=True)
Q=s.diag(1,1,-1)
seat=s.Matrix([0,0,1])
positive_seat_form=Q+2*(Q*seat)*(Q*seat).T
assert positive_seat_form==s.eye(3)
H=s.Matrix([s.cosh(u)*s.cos(phi),s.cosh(u)*s.sin(phi),s.sinh(u)])
record('ruler_depth_metric',s.trigsimp(H.diff(u).dot(H.diff(u))),s.cosh(2*u))
record('ruler_azimuth_metric',s.trigsimp(H.diff(phi).dot(H.diff(phi))),s.cosh(u)**2)
record('ruler_mixed_metric',s.trigsimp(H.diff(u).dot(H.diff(phi))))
Ch,Cg,dt,omega=s.symbols('Ch Cg dt omega',positive=True)
T=Ch*omega**2+Cg*(omega+dt)**2
opt=-Cg*dt/(Ch+Cg)
record('eliminate_common_azimuth',T.subs(omega,opt),Ch*Cg*dt**2/(Ch+Cg))
kinetic=s.diag(s.cosh(2*u),s.cosh(2*v),s.cosh(u)**2*s.cosh(v)**2/(s.cosh(u)**2+s.cosh(v)**2))
metric_boundary=kinetic.subs({u:s.asinh(1),v:0}).applyfunc(lambda e:s.simplify(e.rewrite(s.exp)))
metric_neutral=s.simplify(kinetic.subs({u:0,v:0}))
assert metric_boundary==s.diag(3,1,s.Rational(2,3))
assert metric_neutral==s.diag(1,1,s.Rational(1,2))
results['kinetic_metric_boundary_fixture']={'verified':True,'coordinates':'u=asinh(1), v=0, t=pi/4','diagonal':[str(e) for e in metric_boundary.diagonal()],'status':'positive and finite; candidate metric, not a derived physical law'}
results['kinetic_metric_neutral_fixture']={'verified':True,'coordinates':'u=v=0, t=pi/2','diagonal':[str(e) for e in metric_neutral.diagonal()]}

# C. Active readout inversion under an explicit common-pivot protocol.
C=s.symbols('C',real=True)
p=s.sqrt(1-C*C)
z=s.symbols('z',real=True)
Q=s.diag(1,1,-1)
Ah=s.sqrt(1+a*a); Ag=s.sqrt(1+b*b)
h=s.Matrix([Ah,0,-a]); ruler_g=s.Matrix([Ag*C,Ag*p,-b])
gam=(h.T*Q*ruler_g)[0]
def observed(n):
    ah=(h.T*Q*n)[0]; bg=(ruler_g.T*Q*n)[0]
    return (gam+ah*bg)/s.sqrt((1+ah*ah)*(1+bg*bg))
n1=s.Matrix([z,0,s.sqrt(1+z*z)])
n2=s.Matrix([0,z,s.sqrt(1+z*z)])
m1=s.diff(observed(n1),z).subs(z,0)
m2=s.diff(observed(n2),z).subs(z,0)
record('pivot_readout_x',m1,(1-C*C)*b/Ag)
record('pivot_readout_y',m2,p*(a/Ah-C*b/Ag))
rec_b=m1/(1-C*C)
rec_a=m2/p+C*rec_b
record('recover_bounded_depth_a',rec_a,a/Ah)
record('recover_bounded_depth_b',rec_b,b/Ag)
example={a:s.Rational(3,4),b:s.Rational(4,3),C:s.Rational(3,5)}
record('readout_x_fixture',m1.subs(example),s.Rational(64,125))
record('readout_y_fixture',m2.subs(example),s.Rational(12,125))
record('depth_a_fixture',(rec_a/s.sqrt(1-rec_a**2)).subs(example),s.Rational(3,4))
record('depth_b_fixture',(rec_b/s.sqrt(1-rec_b**2)).subs(example),s.Rational(4,3))
# Jacobian: angle + two slopes identify all three continuous coordinates.
Jac=s.Matrix([C,(1-C*C)*b/Ag,p*(a/Ah-C*b/Ag)]).jacobian([a,b,C])
record('readout_jacobian',Jac.det(),-(1-C*C)**s.Rational(3,2)/(Ah**3*Ag**3))

# D. Symmetry constrains the tangent response without fixing a single wave cone.
# At a=b=gamma=0, selected-pole reversals and ruler exchange act on
# tangent Gram coordinates da,db,dgamma by these matrices.
R1=s.diag(-1,1,-1)
R2=s.diag(1,-1,-1)
Swap=s.Matrix([[0,1,0],[1,0,0],[0,0,1]])
aa,ab,ac,bb,bc,cc=s.symbols('aa ab ac bb bc cc',real=True)
T=s.Matrix([[aa,ab,ac],[ab,bb,bc],[ac,bc,cc]])
eqs=[]
for R in [R1,R2,Swap]: eqs.extend(list(R.T*T*R-T))
sol=s.solve(eqs,[ab,ac,bb,bc],dict=True)
assert sol==[{ab:0,ac:0,bb:aa,bc:0}]
results['neutral_response_symmetry']={'verified':True,'general_symmetric_form':'diag(A,A,D)','free_components':2,'hypotheses':'neutral response invariant under independent pole reversals and ruler exchange at the orthogonal state'}

# Hypothesis: linearized second-order local dynamics about that state, with
# positive kinetic K and spatial coupling C0; isotropic resolved space.
kd,kt,cd,ct,omega,k=s.symbols('k_d k_t c_d c_t omega k',positive=True)
K=s.diag(kd,kd,kt)
C0=s.diag(cd,cd,ct)
principal=-omega**2*K+k*k*C0
record('native_principal_polynomial',principal.det(),(cd*k*k-kd*omega**2)**2*(ct*k*k-kt*omega**2))
v0=s.symbols('v0',positive=True)
record('single_cone_condition',principal.det().subs({cd:v0*v0*kd,ct:v0*v0*kt}),kd*kd*kt*(v0*v0*k*k-omega**2)**3)
results['wave_speed_classes']={'depth_squared':'c_d/k_d (multiplicity 2)','angle_squared':'c_t/k_t','one_cone_iff':'c_d/k_d = c_t/k_t','match_seat_cone_requires':'both equal the squared speed defined by the seat clock and rulers'}

# E. An elimination law that can connect a load to the preceding response.
A0,B0,C1,f,y=s.symbols('A0 B0 C1 f y',nonzero=True)
action=(A0*f*f+2*B0*f*y+C1*y*y)/2
record('load_elimination',action.subs(y,-B0*f/C1),(A0-B0*B0/C1)*f*f/2)
x=s.symbols('x',positive=True)
L=s.diag(x,x*x); drive=s.Matrix([1,1])
readouts=[s.Matrix([[0,1]]),s.Matrix([[1,0]]),s.Matrix([[1,-x]])]
actual=[s.simplify((B*L.inv()*drive)[0]) for B in readouts]
assert actual==[x**-2,x**-1,0]
results['observed_response_orders']={'verified':True,'operator':'diag(x,x^2)','drive':'(1,1)','readouts':['(0,1)','(1,0)','(1,-x)'],'responses':[str(v) for v in actual]}

# F. Native normalization cover over Q(a,b,gamma), without a real-node group.
radicands=[delta,eps,1+a*a,1+b*b]
primes=[delta,1-c,1+a*a,1+b*b]
variables=(a,b,c)
def multiplicity(f,p):
    n=0
    f=s.Poly(f,*variables,domain=s.QQ)
    p=s.Poly(p,*variables,domain=s.QQ)
    while True:
        quotient,remainder=s.div(f,p)
        if not remainder.is_zero: return n
        n+=1; f=quotient
V=s.Matrix([[multiplicity(f,p)%2 for f in radicands] for p in primes])
assert V==s.eye(4)
# Irreducibility of delta as a quadratic in gamma follows from its
# discriminant 4(1+a^2)(1+b^2), which is not a square in Q(a,b).
record('normalization_discriminant',s.discriminant(delta,c),4*(1+a*a)*(1+b*b))
for p0 in primes:
    _,factors=s.factor_list(p0,*variables)
    assert len(factors)==1 and factors[0][1]==1
results['native_normalization_kummer_module']={'verified':True,'base':'Q(a,b,gamma)','radicands':[str(v) for v in radicands],'parity_matrix':[[int(e) for e in row] for row in V.tolist()],'generic_degree':16,'generic_group':'(C2)^4','scope':'algebraic cover of the selected normalizations; not a count of physical states or a click-group identification'}

# Character kernels describe exactly what each displayed scalar can read.
results['readout_characters']={'order_of_radicals':['sqrt(delta)','sqrt(epsilon)','sqrt(1+a^2)','sqrt(1+b^2)'],'N_squared':[0,0,0,0],'N':[1,1,0,0],'visible_cosine':[0,0,1,1],'oriented_volume':[1,0,0,0],'dual_normal':'rational: no radical sign is required'}

print(json.dumps(results,indent=2))
