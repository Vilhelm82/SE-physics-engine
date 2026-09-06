import sympy as sp
from sympy import Rational as R
e1,e2,e3=sp.symbols('e1 e2 e3')
finf=sp.sympify(open('Cb_symmetric.srepr').read())
def C2_of(g12,g13,g23):
    E={e1:g12+g13+g23,e2:g12*g13+g12*g23+g13*g23,e3:g12*g13*g23}
    return sp.nsimplify(finf.subs(E)/(16*(1-E[e3])**4))
# trine
tr=C2_of(R(-1,2),R(-1,2),R(-1,2))
print(f"trine:  C^2 = {tr} = {float(tr):.8f},  C = {float(sp.sqrt(tr)):.8f}")
# positivity sweep over rational boundary points
import itertools
neg=0; tested=0
for t2 in [R(i,9) for i in (-8,-5,-2,1,3,6,8)]:
    for t3 in [R(i,7) for i in (-6,-4,-1,2,5)]:
        if t2==t3: continue
        c2,s2=(1-t2**2)/(1+t2**2),2*t2/(1+t2**2)
        c3,s3=(1-t3**2)/(1+t3**2),2*t3/(1+t3**2)
        v=C2_of(c2,c3,c2*c3+s2*s3); tested+=1
        if v<0: neg+=1
print(f"positivity: {tested-neg}/{tested} boundary samples nonnegative")
# det inertia identity: 1 - e3 = (1/2) det I_inplane on the locus
c2,s2=(R(3,5),R(4,5)); c3,s3=(R(-5,13),R(-12,13))
sig2=(s2)**2+(s3)**2+(c2*s3-s2*c3)**2
e3v=c2*c3*(c2*c3+s2*s3)
print(f"det I_inplane = sum sigma^2 = {sig2} ;  2(1-e3) = {2*(1-e3v)} ;  equal: {sig2==2*(1-e3v)}")
# C at the four verified points for the report table
print(f"A: C = {float(sp.sqrt(R(28981090229,7233948160))):.6f}   B: C = {float(sp.sqrt(R(333432430282250,48648964964161))):.6f}")
