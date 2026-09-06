import sympy as sp
from sympy import Rational as R
def cross(u,v): return sp.Matrix([u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]])
def build_M(a):
    M=sp.zeros(9,9); r=0
    for i in range(3):
        for c in range(3): M[r,3*i+c]=a[i][c]
        r+=1
    for k in range(3):
        for i in range(3):
            for c in range(3):
                M[r,3*i+c]=cross(a[i],sp.Matrix([1 if j==c else 0 for j in range(3)]))[k]
        r+=1
    for (i,j) in [(0,1),(0,2),(1,2)]:
        for c in range(3):
            M[r,3*j+c]+=a[i][c]; M[r,3*i+c]+=a[j][c]
        r+=1
    return M
def inertia(a): return 3*sp.eye(3)-sum((v*v.T for v in a),sp.zeros(3,3))
Z3=sp.zeros(3,1)
def C2(u1,u2,u3):
    a0=[sp.Matrix([u1[0],u1[1],0]),sp.Matrix([u2[0],u2[1],0]),sp.Matrix([u3[0],u3[1],0])]
    M0=build_M(a0)
    N=M0.nullspace()[0]; L=M0.T.nullspace()[0]
    M1=build_M([sp.Matrix([0,0,1]),Z3,Z3])
    LM1N=sp.simplify((L.T*M1*N)[0,0]); I0=inertia(a0)
    Maug=M0.col_join(N.T)                                # 10x9, full column rank
    def psolve(rhs):
        return (Maug.T*Maug).solve(Maug.T*rhs.col_join(sp.zeros(1,1)))
    def T(w): return sum((cross(N[3*i:3*i+3,0],w[3*i:3*i+3,0]) for i in range(3)),Z3)
    cs=[];X0=[]
    for p in range(3):
        b=sp.zeros(9,1); b[6+p]=1
        c=sp.simplify((L.T*b)[0,0]/LM1N); cs.append(c)
        X0.append(psolve(b-c*M1*N))
    tot=0
    for p,q in [(0,1),(0,2),(1,2)]:
        S=cs[p]*T(X0[q])-cs[q]*T(X0[p])
        F=2*I0.solve(S); tot+=(F.T*F)[0,0]
    wedge=u2[0]*u3[1]-u2[1]*u3[0]
    return sp.simplify(tot*wedge**2), N, L
A,NA,LA=C2(sp.Matrix([1,0]),sp.Matrix([R(3,5),R(4,5)]),sp.Matrix([R(-5,13),R(-12,13)]))
print(f"config A: C^2 = {A} = {float(A):.10f}   Richardson gave 4.0062618062")
B,_,_=C2(sp.Matrix([1,0]),sp.Matrix([R(4,5),R(3,5)]),sp.Matrix([R(-12,13),R(5,13)]))
print(f"config B: C^2 = {B} = {float(B):.10f}   Richardson gave 6.8538442185")
# structural closed forms for the kernel data
print("\nkernel structure at config A (theory: N_i = sin(th_j - th_k) cyclic, out-of-plane):")
s23=R(3,5)*R(-12,13)-R(4,5)*R(-5,13)      # u2 ^ u3 = sin(th3-th2)
s31=R(-5,13)*0-R(-12,13)*1                # u3 ^ u1 = sin(th1-th3)
s12=1*R(4,5)-0*R(3,5)                     # u1 ^ u2 = sin(th2-th1)
zN=[sp.simplify(NA[2]),sp.simplify(NA[5]),sp.simplify(NA[8])]
pred=[s23,s31,s12]
rat=[sp.simplify(zN[i]/pred[i]) for i in range(3)]
print(f"  N z-components      : {zN}")
print(f"  sin(angle diffs)    : {pred}")
print(f"  ratios (must match) : {rat}")
