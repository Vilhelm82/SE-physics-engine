#!/usr/bin/env python3
"""Native pinch-rate selection, area transport, and passive-load conditions.

Dependencies and proofs: docs/results/2026-09-05/2026-09-05-pinch-load-transport.md.
The illustrative loads are declared candidates, not deductions of a unique load.
Run this file to emit the exact-check record as JSON; failures raise exceptions.
"""
import json
import sympy as s

checks = {}


def check(name, actual, expected=0):
    diff = actual-expected
    entries = list(diff) if isinstance(diff,s.MatrixBase) else [diff]
    for entry in entries:
        residual = s.factor(s.cancel(s.expand(entry)))
        if residual != 0:
            residual = s.simplify(residual)
        if residual != 0:
            raise AssertionError((name,residual))
    checks[name] = {"verified": True}


# Native directions and smooth crossing: the second jet determines the rate
# of the quotient at the crossing; its first jet determines the limiting normal.
t,a,xi,d,x = s.symbols("t a xi d x",real=True)
p0,p1,q0,q1 = s.symbols("p0 p1 q0 q1",nonzero=True,real=True)
quotient = (p0*t+p1*t*t/2)/(q0*t+q1*t*t/2)
check("crossing_direction",s.limit(quotient,t,0),p0/q0)
check("crossing_direction_rate",s.limit(s.diff(quotient,t),t,0),
      (p1*q0-p0*q1)/(2*q0*q0))
F = (1-xi*xi)/(1+a*a)
adot,xidot = s.symbols("a_dot xi_dot",real=True)
check("native_lapse_rate",s.diff(F,a)*adot+s.diff(F,xi)*xidot,
      -2*xi*xidot/(1+a*a)-2*a*F*adot/(1+a*a))

l1,l2,ang,ell = s.symbols("l1 l2 angle ell",real=True)
v1,v2,vt = s.symbols("l1_dot l2_dot angle_dot",real=True)
aa,bb = -s.sinh(l1),-s.sinh(l2)
gam = s.cosh(l1)*s.cosh(l2)*s.cos(ang)-s.sinh(l1)*s.sinh(l2)
xx = bb-aa*gam
dd = s.cosh(l1)*s.cosh(l2)*s.sin(ang)
def native_rate(expr):
    return sum(s.diff(expr,z)*dz for z,dz in ((l1,v1),(l2,v2),(ang,vt)))
for sigma,angle in ((1,0),(-1,s.pi)):
    sub = {l1:ell,l2:sigma*ell,ang:angle}
    P = s.simplify(native_rate(xx).subs(sub))
    Q = s.simplify(native_rate(dd).subs(sub))
    check(f"native_depth_rate_{sigma}",P,s.cosh(ell)*(sigma*v1-v2))
    check(f"native_volume_rate_{sigma}",Q,sigma*s.cosh(ell)**2*vt)
    check(f"native_selected_direction_{sigma}",P/Q,
          (v1-sigma*v2)/(s.cosh(ell)*vt))

# A pinch-preserving vector field induces a projective Riccati equation.
m11,m12,m21,m22 = s.symbols("m11 m12 m21 m22",real=True)
P=m11*x+m12*d
Q=m21*x+m22*d
riccati=m12+(m11-m22)*xi-m21*xi*xi
check("pinch_riccati",((P-xi*Q)/d).subs(x,xi*d),riccati)
check("pinch_volume_rate",Q.subs(x,xi*d),d*(m21*xi+m22))
for tau in (1,-1):
    check(f"null_direction_drift_{tau}",riccati.subs(xi,tau),
          m12-m21+tau*(m11-m22))
check("both_null_directions_invariant",
      riccati.subs({m22:m11,m21:m12}),m12*(1-xi*xi))

# Explicit reciprocal passive relaxation, with declared quadratic storage.
kappa,beta,x0,d0 = s.symbols("kappa beta x0 d0",real=True)
xt=s.exp(-kappa*t)*(x0*s.cosh(beta*t)+d0*s.sinh(beta*t))
dt=s.exp(-kappa*t)*(d0*s.cosh(beta*t)+x0*s.sinh(beta*t))
check("passive_selector_x_solution",s.diff(xt,t),-kappa*xt+beta*dt)
check("passive_selector_d_solution",s.diff(dt,t),beta*xt-kappa*dt)
check("passive_selector_projective_flow",
      ((-kappa*x+beta*d)-xi*(beta*x-kappa*d)).subs(x,xi*d)/d,
      beta*(1-xi*xi))
power=kappa*(x*x+d*d)-2*beta*x*d
check("passive_selector_power",power,
      (kappa-beta)*(x+d)**2/2+(kappa+beta)*(x-d)**2/2)
mobility=s.Matrix([[kappa,-beta],[-beta,kappa]])
check("passive_selector_eigenvalues",mobility.charpoly().as_expr(),
      (s.Symbol("lambda")-kappa+beta)*(s.Symbol("lambda")-kappa-beta))
check("passive_selector_lapse_rate",s.diff(F,xi)*beta*(1-xi*xi),-2*beta*xi*F)

# Same applied effort and same input power admit all three pinch causal types.
u,ex,ed = s.symbols("u e_x e_d",real=True)
Mu=s.Matrix([[1+u*u,u],[u,1]])
eu=s.Matrix([ex,ed])
check("passive_direction_family_determinant",Mu.det(),1)
check("passive_direction_family_power",(eu.T*Mu*eu)[0],ex*ex+(u*ex+ed)**2)
check("passive_direction_family_fixed_effort",Mu*s.Matrix([0,1]),s.Matrix([u,1]))
check("passive_direction_family_fixed_input_power",(s.Matrix([0,1]).T*Mu*s.Matrix([0,1]))[0],1)

# The physical clock C and the native rest normal K are distinct fields.
w,alpha,beta_area = s.symbols("w alpha beta_area",real=True)
D=alpha+w*beta_area
rho=(beta_area-alpha)/(1+w)
theta=(1+w)*(alpha+beta_area)/2
check("clock_area_bridge",theta,D+rho*(1-w*w)/2)
check("inward_area_expansion",(alpha-beta_area)/(1+w),-rho)
Xi=s.symbols("Xi",real=True)
check("clock_area_matching_law",D.subs(alpha,-w*beta_area+(1-w*w)*Xi),
      (1-w*w)*Xi)

# A general two-port-plus-area passive response, including a reactive part.
# Positivity is a theorem hypothesis on R0, not inferred from symbolic entries.
f=s.symbols("F",real=True)
r11,r12,r22,c1,c2,r0,j12,j1,j2=s.symbols("r11 r12 r22 c1 c2 r0 j12 j1 j2",real=True)
R0=s.Matrix([[r11,r12,c1],[r12,r22,c2],[c1,c2,r0]])
J0=s.Matrix([[0,j12,j1],[-j12,0,j2],[-j1,-j2,0]])
S=s.diag(1,1,f)
M=S*(R0+J0)*S
e=s.Matrix(s.symbols("e1 e2 e_area",real=True))
check("passive_response_power",(e.T*M*e)[0],((S*e).T*R0*(S*e))[0])
check("area_response_row",(M*e)[2],
      f*((c1-j1)*e[0]+(c2-j2)*e[1]+f*r0*e[2]))
check("quadratic_area_self_response",M[2,2],f*f*r0)
check("area_boundary_annihilation",M.subs(f,0)[2,:],s.zeros(1,3))
check("unrestricted_orientation_row",(M*e)[2].subs({j1:c1,j2:c2}),f*f*r0*e[2])
rho0=s.symbols("rho0",positive=True)
check("robust_boundary_factor",(1+2*f*r0*e[2]/rho0).subs(f,0),1)

# Full finite-time realization.  C=partial_t remains timelike through t=0.
# Native load: dot d=1, dot x=1+2 k d, dot a=0.
k=s.symbols("k",positive=True)
r,th,ph=s.symbols("r theta phi",real=True)
xi_t=1+k*t
d_t=t
x_t=t+k*t*t
F_t=1-xi_t*xi_t
R=r-t-k*t*t/2
check("crossing_witness_native_x_rate",s.diff(x_t,t),1+2*k*d_t)
check("crossing_witness_native_d_rate",s.diff(d_t,t),1)
check("crossing_witness_direction",s.cancel(x_t/d_t),xi_t)
check("crossing_witness_selected_horizon",xi_t.subs(t,0),1)
check("crossing_witness_selected_direction_rate",s.diff(xi_t,t),k)
area_mobility=s.symbols("area_mobility",positive=True)
load=Mu.subs(u,1+2*k*d).row_join(s.zeros(2,1)).col_join(s.Matrix([[0,0,area_mobility]]))
check("crossing_witness_declared_passive_load",load*s.Matrix([0,1,0]),
      s.Matrix([1+2*k*d,1,0]))
check("crossing_witness_load_nondegenerate",load.det(),area_mobility)
gamma_t=s.sqrt(1-F_t*d_t*d_t)
original_gram=s.Matrix([[-1,0,x_t],[0,1,gamma_t],[x_t,gamma_t,1]])
check("crossing_witness_original_determinant",original_gram.det(),-t*t)
check("crossing_witness_original_lapse",1-gamma_t*gamma_t,F_t*t*t)

# Independent soldering of the resolved frame into this spacetime metric.
Gres=s.Matrix([[-1,0,xi_t,0],[0,1,0,0],
               [xi_t,0,F_t,0],[0,0,0,1]])
E=s.Matrix([[1,-xi_t,0,0],[0,0,R,0],[0,-1,0,0],
            [0,0,0,R*s.sin(th)]])
g=s.diag(-1,1,R*R,R*R*s.sin(th)**2)
coords=(t,r,th,ph)
check("crossing_witness_soldering",E.T*Gres*E,g)
check("crossing_witness_resolved_determinant",Gres.det(),-1)
K=s.Matrix([1,xi_t,0,0])
L=(1+xi_t)*s.Matrix([1,1,0,0])/2
N=s.Matrix([1,-1,0,0])/(1+xi_t)
check("crossing_witness_native_normal",E*K,s.Matrix([F_t,0,-xi_t,0]))
check("crossing_witness_outgoing_null",(L.T*g*L)[0])
check("crossing_witness_inward_null",(N.T*g*N)[0])
check("crossing_witness_null_pairing",(L.T*g*N)[0],-1)

def screen_expansion(V):
    return s.factor(sum((sum(V[c]*s.diff(g[A,A],coords[c])+
                      2*g[c,A]*s.diff(V[c],coords[A]) for c in range(4)))/g[A,A]
                     for A in (2,3))/2)

theta_L=screen_expansion(L)
theta_N=screen_expansion(N)
check("crossing_witness_native_area_transport",screen_expansion(K))
check("crossing_witness_outgoing_expansion",theta_L,F_t/R)
check("crossing_witness_inward_expansion",theta_N,-2/R)
def directional(V,expr):
    return s.factor(sum(V[c]*s.diff(expr,coords[c]) for c in range(4)))
check("crossing_witness_outer_derivative",directional(N,theta_L).subs(t,0),-k/r)
check("crossing_witness_regular_transition",s.diff(F_t,t).subs(t,0),-2*k)

# Compute the complete Ricci tensor from the metric, rather than assigning it.
ginv=g.inv()
Gamma={}
for i in range(4):
    for j in range(4):
        for l in range(4):
            Gamma[i,j,l]=s.factor(sum(ginv[i,m]*(s.diff(g[m,l],coords[j])+
                              s.diff(g[m,j],coords[l])-s.diff(g[j,l],coords[m]))
                              for m in range(4))/2)
Ric=s.zeros(4)
for i in range(4):
    for j in range(4):
        Ric[i,j]=s.simplify(sum(s.diff(Gamma[m,i,j],coords[m])-
                         s.diff(Gamma[m,i,m],coords[j])+
                         sum(Gamma[m,i,j]*Gamma[n,m,n]-Gamma[n,i,m]*Gamma[m,j,n]
                             for n in range(4)) for m in range(4)))
Ric_expected=s.diag(2*k/R,0,xi_t**2-k*R,(xi_t**2-k*R)*s.sin(th)**2)
check("crossing_witness_full_ricci",Ric,Ric_expected)
RLL=(L.T*Ric*L)[0]
check("crossing_witness_null_focusing",RLL,k*(1+xi_t)**2/(2*R))
acceleration=s.Matrix([sum(L[j]*s.diff(L[i],coords[j]) for j in range(4))+
                       sum(Gamma[i,j,l]*L[j]*L[l] for j in range(4) for l in range(4))
                       for i in range(4)])
check("crossing_witness_pregeodesic_normal",acceleration,k*L/2)
vr,angular_norm=s.symbols("v_r angular_norm",real=True)
null_Ric=2*k/R*(vr*vr+R*R*angular_norm)+(xi_t**2-k*R)*angular_norm
check("crossing_witness_all_null_convergence",null_Ric,
      2*k*vr*vr/R+(xi_t**2+k*R)*angular_norm)
check("crossing_witness_raychaudhuri",directional(L,theta_L),
      k*theta_L/2-theta_L**2/2-RLL)
check("crossing_witness_horizon_rate_balance",
      (directional(L,xi_t)-RLL/(2/R)).subs(t,0))
focus,LRXi=s.symbols("focus L_xi",real=True)
for tau in (1,-1):
    predicted=tau*(1+a*a)*focus/(rho0+2*Xi)
    check(f"general_horizon_focusing_balance_{tau}",
          (-tau*(rho0+2*Xi)*LRXi/(1+a*a)+focus).subs(LRXi,predicted))

samples=[]
for tt in (-s.Rational(1,4),s.Integer(0),s.Rational(1,4)):
    sub={t:tt,k:1,r:2}
    vals={"clock":str(tt),"oriented_volume":str(d_t.subs(sub)),
          "xi":str(xi_t.subs(sub)),"F":str(F_t.subs(sub)),
          "theta_L":str(s.factor(theta_L.subs(sub))),
          "theta_N":str(s.factor(theta_N.subs(sub)))}
    samples.append(vals)
check("finite_time_compact_side",F_t.subs({k:1,t:-s.Rational(1,4)}),s.Rational(7,16))
check("finite_time_trapped_side",F_t.subs({k:1,t:s.Rational(1,4)}),-s.Rational(9,16))

print(json.dumps({
    "scope":"Native rate-selection and transport identities; conditional passive-response classification; an explicit declared load with a regular finite-time trapping/pinch crossing. No unique native load is inferred.",
    "verified_count":len(checks),
    "checks":checks,
    "results":{
        "finite_crossing":"xi_0 = dot x_0/dot d_0 when dot d_0 != 0",
        "native_angle_rates":"xi_0 = (dot l1 - sigma dot l2)/(cosh(l1) dot t_angle)",
        "invariant_pinch_flow":"dot xi = m12 + (m11-m22) xi - m21 xi^2 on the resolved pinch",
        "area_clock_law":"theta_C = -sqrt(1-F) theta_R + F Xi",
        "passive_matching_class":"M = diag(I,F) (R0+J0) diag(I,F), R0 symmetric PSD, J0 skew",
        "area_self_response":"M_area,area = F^2 r0 under smooth two-sided passivity and all-effort boundary matching",
        "finite_time_witness":"d=t; x=t+k t^2; xi=1+k t; R=r-t-k t^2/2; theta_L=F/R; theta_N=-2/R",
        "native_remaining_input":"Derive the actual effort, response, and area coupling from the seat/load construction."
    },
    "finite_time_witness_samples":samples
},indent=2))
