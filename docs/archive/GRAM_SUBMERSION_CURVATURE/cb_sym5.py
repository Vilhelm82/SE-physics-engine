import sympy as sp, itertools
from sympy import Rational as R
x,y,z=sp.symbols('g12 g13 g23')
Pred=sp.sympify(open('Pred.srepr').read())
Psym=sp.expand(sum(Pred.subs(dict(zip((x,y,z),p)),simultaneous=True)
                   for p in itertools.permutations((x,y,z)))/6)
res=sp.polys.polyfuncs.symmetrize(Psym,x,y,z,formal=True)
sym,rem=res[0],res[1]
assert rem==0
defs={s:e for s,e in res[2]}
e1,e2,e3=sp.symbols('e1 e2 e3')
name={sp.Symbol('s1'):e1,sp.Symbol('s2'):e2,sp.Symbol('s3'):e3}
expr=sym.subs(name)
DeltaE=sp.expand(1-e1**2+2*e2+2*e3)
P3=sp.Poly(sp.expand(expr),e3)
Q3=sp.Poly(DeltaE,e3)
q_,r_=sp.div(P3,Q3)                       # over QQ(e1,e2)[e3]
fin=sp.simplify(sp.together(r_.as_expr()))
print("degree in e3:",sp.degree(fin,e3))
finf=sp.factor(sp.cancel(fin))
print("\n16 (1 - e3)^4 C^2  =")
sp.pprint(finf,wrap_line=False)
cb,sb,cg,sg=sp.symbols('c_b s_b c_g s_g')
ok=True
tests=[("A",{cb:R(3,5),sb:R(4,5),cg:R(-5,13),sg:R(-12,13)},R(28981090229,7233948160)),
       ("B",{cb:R(4,5),sb:R(3,5),cg:R(-12,13),sg:R(5,13)},R(333432430282250,48648964964161))]
for nm,sub,target in tests:
    g12=sub[cb];g13=sub[cg];g23=sub[cb]*sub[cg]+sub[sb]*sub[sg]
    E={e1:g12+g13+g23,e2:g12*g13+g12*g23+g13*g23,e3:g12*g13*g23}
    got=finf.subs(E)/(16*(1-E[e3])**4)
    m=sp.simplify(got-target)==0
    ok&=m; print(f"config {nm}: matches exact recipe: {m}")
if ok: open('Cb_symmetric.srepr','w').write(sp.srepr(finf))
