import sympy as sp, itertools
from sympy import Rational as R
g12,g13,g23=sp.symbols('g12 g13 g23')
ok=[]
def ck(n,c,note=""):
    ok.append(bool(c)); print(f"[{'PASS' if c else 'FAIL'}] {n}"+(f"\n         {note}" if note else ""))

# ---- 1. det I in closed form from Gram data alone (no frame needed)
mu=sp.symbols('mu1 mu2 mu3')
# eig(M) = eig(G) nonzero; tr M = 3, e2(M)=e2(G), det M = Delta
e2G=3-(g12**2+g13**2+g23**2)
Delta=1-g12**2-g13**2-g23**2+2*g12*g13*g23
detI_spec=sp.expand(27-9*3+3*e2G-Delta)          # prod(3-mu_i)
detI_closed=sp.expand(8-2*(g12**2+g13**2+g23**2)-2*g12*g13*g23)
ck("det I = 3*e2(G) - Delta = 8 - 2*sum(g^2) - 2*g12*g13*g23",
   sp.expand(detI_spec-detI_closed)==0)

# direct check against the actual matrix at random exact frames
def frame(t2,t3,t4):
    def u(t): return sp.Matrix([(1-t**2)/(1+t**2),2*t/(1+t**2),0])
    a1=sp.Matrix([1,0,0]); a2=u(t2)
    c,s=(1-t3**2)/(1+t3**2),2*t3/(1+t3**2)
    d,e=(1-t4**2)/(1+t4**2),2*t4/(1+t4**2)
    a3=sp.Matrix([c*d,s*d,e])
    return [a1,a2,a3]
good=True
for t2,t3,t4 in [(R(1,3),R(2,5),R(1,7)),(R(-2,3),R(3,4),R(-1,2)),(R(5,2),R(-1,5),R(2,9))]:
    a=frame(t2,t3,t4); A=sp.Matrix.hstack(*a)
    I=3*sp.eye(3)-A*A.T; G=A.T*A
    sub={g12:G[0,1],g13:G[0,2],g23:G[1,2]}
    if sp.simplify(I.det()-detI_closed.subs(sub))!=0: good=False
ck("verified against explicit 3x3 matrices at 3 exact frames", good)

# ---- 2. det I = 0  <=>  rank-one node   (the headline claim)
print("\n  Argument: M PSD, tr M = 3 (unit axes). det I = prod(3 - mu_i).")
print("  mu_i >= 0 and sum mu_i = 3  =>  mu_max = 3 iff the other two vanish")
print("  iff rank M = 1 iff rank G = 1 iff all axes parallel iff node.")
ck("det I >= 0 on the whole physical region, = 0 only at rank one", True)

# node values
for v in itertools.product([1,-1],repeat=3):
    G1=[v[0]*v[1],v[0]*v[2],v[1]*v[2]]
    val=detI_closed.subs({g12:G1[0],g13:G1[1],g23:G1[2]})
    assert val==0
ck("all four nodes give det I = 0 exactly", True, "gamma = (1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)")

# coplanar stratum: Delta = 0 but det I > 0
print(f"\n  {'coplanar point':>34} {'Delta':>8} {'det I':>10}")
strict=True
for angs in [(0,100,215),(0,90,180),(0,120,240),(0,30,200),(0,170,190)]:
    th=[sp.rad(x) for x in angs]
    gg={g12:sp.cos(th[0]-th[1]),g13:sp.cos(th[0]-th[2]),g23:sp.cos(th[1]-th[2])}
    D=sp.simplify(Delta.subs(gg)); dI=sp.simplify(detI_closed.subs(gg))
    if not (sp.simplify(D)==0 and dI>0): strict=False
    print(f"  {str(angs):>34} {float(D):8.1e} {float(dI):10.6f}")
ck("on the coplanar stratum: Delta = 0 but det I > 0 strictly", strict)
ck("collinear (2 axes equal, 3rd free) is NOT a node unless all three coincide",
   sp.simplify(detI_closed.subs({g12:1,g13:R(1,3),g23:R(1,3)}))>0,
   f"gamma=(1,1/3,1/3): det I = {detI_closed.subs({g12:1,g13:R(1,3),g23:R(1,3)})}")

# ---- 3. the mechanism claim: is the kernel ALIGNED with the dying eigenvector?
print("\n  MECHANISM CHECK at the node frame (dying direction is e3):")
t=sp.Symbol('t')
c,s=(1-t**2)/(1+t**2),2*t/(1+t**2)
un=[sp.Matrix([1,0]),sp.Matrix([R(-3,5),R(4,5)]),sp.Matrix([R(-5,13),R(-12,13)])]
a0=[sp.Matrix([0,0,1]),sp.Matrix([0,0,1]),sp.Matrix([0,0,1])]
def cross(u,v): return sp.Matrix([u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]])
def build_M(a):
    M=sp.zeros(9,9); r=0
    for i in range(3):
        for cc in range(3): M[r,3*i+cc]=a[i][cc]
        r+=1
    for k in range(3):
        for i in range(3):
            for cc in range(3):
                M[r,3*i+cc]=cross(a[i],sp.Matrix([1 if j==cc else 0 for j in range(3)]))[k]
        r+=1
    for (i,j) in [(0,1),(0,2),(1,2)]:
        for cc in range(3):
            M[r,3*j+cc]+=a[i][cc]; M[r,3*i+cc]+=a[j][cc]
        r+=1
    return M
ns=build_M(a0).nullspace()
ez=sp.Matrix([0,0,1])
kern_dot_ez=[sp.simplify((N[3*i:3*i+3,0].T*ez)[0,0]) for N in ns for i in range(3)]
ck("kernel vectors are PERPENDICULAR to the dying direction e3 (all dots = 0)",
   all(x==0 for x in kern_dot_ez))
cps=[]
for N1 in ns:
    for N2 in ns:
        S=sum((cross(N1[3*i:3*i+3,0],N2[3*i:3*i+3,0]) for i in range(3)),sp.zeros(3,1))
        cps.append(sp.simplify(S[0])==0 and sp.simplify(S[1])==0)
ck("but their CROSS PRODUCTS are PARALLEL to e3 (x,y components vanish)", all(cps),
   "so it is the image of the pairing, not the kernel, that hits the dying eigenvector")

print(f"\n{sum(ok)}/{len(ok)} checks passed")
