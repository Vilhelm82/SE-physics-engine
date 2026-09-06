import sympy as sp
from sympy import Rational as R
exec(open('cb_symbolic2.py').read().split('# EXACT verification')[0])

x,y,z=sp.symbols('g12 g13 g23')
def onshell(expr):
    """exact rewrite: cb->x, cg->y, sb^2->1-x^2, sg^2->1-y^2, sb*sg->z-x*y"""
    P=sp.Poly(sp.expand(expr),cb,sb,cg,sg)
    out=0
    for (a,p,b,q),coef in P.terms():
        if (p+q)%2: raise ValueError("odd parity term - reflection invariance violated")
        r=p%2
        out+=coef*x**a*y**b*(1-x**2)**((p-r)//2)*(1-y**2)**((q-r)//2)*(z-x*y)**r
    return sp.expand(out)

num,den=sp.fraction(C2)
D=onshell(den)
print("denominator on-shell equals 16(1-g12 g13 g23)^4:",
      sp.simplify(D-16*(1-x*y*z)**4)==0)
Pn=onshell(num)
Delta=1-x**2-y**2-z**2+2*x*y*z
q_,r_=sp.div(Pn,Delta,x)                     # reduce mod the boundary relation
Pred=sp.expand(r_)
print("degree of numerator after reduction mod Delta:",sp.total_degree(Pred))
Pf=sp.factor(Pred)
print("\nC^2 * 16(1-g12g13g23)^4  =  (mod Delta=0)")
sp.pprint(Pf,wrap_line=False)

# exact verification at configs A and B and one more random boundary point
for nm,sub,target in [("A",{cb:R(3,5),sb:R(4,5),cg:R(-5,13),sg:R(-12,13)},R(28981090229,7233948160)),
                      ("B",{cb:R(4,5),sb:R(3,5),cg:R(-12,13),sg:R(5,13)},R(333432430282250,48648964964161))]:
    g12=sub[cb]; g13=sub[cg]; g23=sub[cb]*sub[cg]+sub[sb]*sub[sg]
    got=Pred.subs({x:g12,y:g13,z:g23})/(16*(1-g12*g13*g23)**4)
    print(f"config {nm}: closed form matches exact 9x9 recipe: {sp.simplify(got-target)==0}")
t2,t3=R(2,7),R(-5,9)
c2v,s2v=(1-t2**2)/(1+t2**2),2*t2/(1+t2**2); c3v,s3v=(1-t3**2)/(1+t3**2),2*t3/(1+t3**2)
sub={cb:c2v,sb:s2v,cg:c3v,sg:s3v}
g12,g13,g23=c2v,c3v,c2v*c3v+s2v*s3v
lhs=sp.nsimplify(num.subs(sub)/den.subs(sub))
got=Pred.subs({x:g12,y:g13,z:g23})/(16*(1-g12*g13*g23)**4)
print(f"random boundary point: closed form matches component formula: {sp.simplify(got-lhs)==0}")
