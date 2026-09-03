#!/usr/bin/env python3
# =============================================================================
# SECT-1 - THE TRICHOTOMY FROM ONE PAIRING  (2026-09-03)
#
# CLAIM: the three seat characters (free / floored / horizon-vanishing) are the
# three regimes of ONE bilinear form - RULE-1's bar pairing <M Mbar> on Cl(3) -
# and the hbar seat is DERIVED (the c seat read through its spinor square root),
# with the transport census (2,1,2) following from two-sided vs one-sided action.
# Will's criterion (2026-09-03): a seat whose behaviour must be forced in beyond
# its name is a failure signal. This suite checks that only a REGIME (a sign) is
# assigned per seat and every behaviour follows.
#
# INPUTS: BARE-1 (rotor); RULE-1 (bar pairing); E-1 (null ray); SEAT-hb (the
#   hbar state is a positive quadratic form - the regime label); KIN-2/e2 pinning
#   (the G regime label); thm_g G-3a/G-8 (squeeze, dictionary - RE-DERIVED here).
# COORDINATES rational: w = e^Lambda, w2 = e^l, t = tan phi.
# TWO PATHS: hand-rolled Cl(3) vs 2x2 quadratic-form matrices; nothing shared.
# =============================================================================
import sympy as sp, time
T0=time.time(); CH=[]
def check(t_,ok,n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t_}"+(f" - {n}" if n else ""), flush=True)
def blade_mul(a,b):
    seq=list(a)+list(b); sg=1
    for i in range(1,len(seq)):
        j=i
        while j>0 and seq[j-1]>seq[j]:
            seq[j-1],seq[j]=seq[j],seq[j-1]; sg=-sg; j-=1
    out,k=[],0
    while k<len(seq):
        if k+1<len(seq) and seq[k]==seq[k+1]: k+=2
        else: out.append(seq[k]); k+=1
    return tuple(out),sg
BL=[(),(1,),(2,),(3,),(1,2),(1,3),(2,3),(1,2,3)]
class MV:
    def __init__(s,d=None):
        s.d={}
        for k,v in (d or {}).items():
            v=sp.sympify(v)
            if v!=0: s.d[tuple(k)]=s.d.get(tuple(k),0)+v
    def __add__(s,o):
        d=dict(s.d)
        for k,v in o.d.items(): d[k]=d.get(k,0)+v
        return MV(d)
    def __sub__(s,o):
        d=dict(s.d)
        for k,v in o.d.items(): d[k]=d.get(k,0)-v
        return MV(d)
    def __mul__(s,o):
        if not isinstance(o,MV): return MV({k:v*o for k,v in s.d.items()})
        d={}
        for ka,va in s.d.items():
            for kb,vb in o.d.items():
                key,sg=blade_mul(ka,kb); d[key]=d.get(key,0)+sg*va*vb
        return MV(d)
    __rmul__=__mul__
    def bar(s): SG={0:1,1:-1,2:-1,3:1}; return MV({k:SG[len(k)]*v for k,v in s.d.items()})
    def rev(s): SG={0:1,1:1,2:-1,3:-1}; return MV({k:SG[len(k)]*v for k,v in s.d.items()})
    def c(s,key): return s.d.get(tuple(key),0)
    def grades(s): return sorted({len(k) for k,v in s.d.items() if sp.simplify(v)!=0})
def rz(e): return sp.expand(sp.numer(sp.cancel(sp.together(sp.expand(e)))))==0
def mvz(M): return all(rz(v) for v in M.d.values())
ONE=MV({():1}); E1,E2,E3=MV({(1,):1}),MV({(2,):1}),MV({(3,):1})
w,w2,t = sp.symbols('w w2 t', positive=True)
ch,sh = (w+1/w)/2,(w-1/w)/2                     # cosh, sinh of Lambda
A = MV({():(sp.sqrt(w)+1/sp.sqrt(w))/2,(1,):(sp.sqrt(w)-1/sp.sqrt(w))/2})  # half-rapidity rotor
def norm(M): return (M.bar()*M).c(())

print("="*78); print("PART 1 - one pairing, its signature, its invariances (Cl(3))"); print("="*78)
sig=[]
for k in BL:
    b=MV({k:1}); sig.append(int(sp.simplify(norm(b))))
para=sig[:4]; wick=sig[4:]
check("S-1  bar signature (4,4); paravector sector (1,3) = MINKOWSKI; Wick sector its negative",
      sig.count(1)==4 and para==[1,-1,-1,-1] and wick==[1,1,1,-1],
      f"per blade {sig}")
sigr=[int(sp.simplify((MV({k:1}).rev()*MV({k:1})).c(()))) for k in BL]
check("S-2  reversion pairing is (8,0) - the contrast: only the BAR pairing carries regimes", sigr==[1]*8)
P=[sp.Symbol(f'p{i}',real=True) for i in range(8)]
psi=MV({k:P[i] for i,k in enumerate(BL)})
check("S-3  Abar A = 1 and the one-sided action preserves psibar psi for GENERAL psi",
      mvz(A.bar()*A-ONE) and rz(norm(A*psi)-norm(psi)))
tt,v1,v2,v3=sp.symbols('tt v1 v2 v3',real=True)
X=MV({():tt,(1,):v1,(2,):v2,(3,):v3})
check("S-4  paravector norm is the scalar t^2-|v|^2 and is two-sided invariant",
      (X.bar()*X).grades()==[0] and rz(norm(X)-(tt**2-v1**2-v2**2-v3**2)) and rz(norm(A*X*A)-norm(X)))
check("S-5  psi psi~ -> A (psi psi~) A: the spinor square root intertwines half -> full",
      mvz((A*psi)*(A*psi).rev() - A*(psi*psi.rev())*A))

print(); print("="*78); print("PART 2 - three seats = three regimes of the ONE norm"); print("="*78)
ct,st=(1-t**2)/(1+t**2),2*t/(1+t**2)
Xc=ONE+ct*E1+st*E2
check("T-1  c seat, regime NULL: <X Xbar> = 0 on-state and after the pivot (E-1) -> FREE",
      rz(norm(Xc)) and rz(norm(A*Xc*A)))
Af=sp.Symbol('A0',positive=True)
Fvac=MV({():Af}); Fp=A*Fvac*A
check("T-2a hbar seat, regime TIMELIKE: the vacuum F = A0 boosts to A0(coshL + sinhL e1)",
      mvz(Fp-MV({():Af*ch,(1,):Af*sh})), "isotropic form acquires SHAPE under the pivot")
check("T-2b ... and its norm is pinned to A0^2 for every rapidity: FLOORED",
      rz(norm(Fp)-Af**2), "the floor is the invariant of a timelike paravector, not a declaration")
th=(w**2-1)/(w**2+1)                                # tanh Lambda
Pg=ONE+th*E1
check("T-3  G seat, regime TIMELIKE PINNED: <P Pbar> = 1-tanh^2 = sech^2 = (1+v)(1-v) -> 0 at the wall",
      rz(norm(Pg)-(1-th**2)) and rz(norm(Pg)-(1+th)*(1-th)) and sp.limit((1-th**2),w,sp.oo)==0,
      "section 6's protected product IS the bar norm of the pinned paravector")

print(); print("="*78); print("PART 3 - the hbar dictionary, DERIVED (2x2 forms; shares nothing with Cl(3))"); print("="*78)
c11,c12,c22=sp.symbols('c11 c12 c22',real=True)
Cm=sp.Matrix([[c11,c12],[c12,c22]])
Fmap=lambda C: ((C[0,0]+C[1,1])/2,(C[0,0]-C[1,1])/2,C[0,1])   # form -> paravector (f0,f1,f2)
f0,f1,f2=Fmap(Cm)
check("D-1  F(C) := ((C11+C22)/2, (C11-C22)/2, C12) has bar-norm f0^2-f1^2-f2^2 = det C EXACTLY",
      rz(f0**2-f1**2-f2**2-Cm.det()), "the Robertson-Schroedinger product is the Minkowski norm")
S=sp.diag(w2,1/w2)                                  # squeeze by rapidity l
g0,g1,g2=Fmap(S*Cm*S.T)
ch2,sh2=(w2**2+w2**-2)/2,(w2**2-w2**-2)/2           # cosh 2l, sinh 2l
check("D-2  squeeze by l on the form = BOOST by 2l on F: RAPIDITY DOUBLING is an identity",
      rz(g0-(ch2*f0+sh2*f1)) and rz(g1-(sh2*f0+ch2*f1)) and rz(g2-f2), "G-8's Lambda = 2 lambda, derived")
cf,sf=(1-t**2)/(1+t**2),2*t/(1+t**2)               # cos phi, sin phi
R=sp.Matrix([[cf,-sf],[sf,cf]])
h0,h1,h2=Fmap(R*Cm*R.T)
c2f,s2f=cf**2-sf**2,2*sf*cf
check("D-3  rotation by phi on the form = rotation by 2phi on (f1,f2): ANGLE DOUBLING is an identity",
      rz(h0-f0) and rz(h1-(c2f*f1-s2f*f2)) and rz(h2-(s2f*f1+c2f*f2)), "G-8's Theta = 2 phi, derived")
check("D-4  det C invariant under the squeeze: the floor survives every pivot",
      rz((S*Cm*S.T).det()-Cm.det()))

print(); print("="*78); print("PART 4 - the transport census from ONE rotor"); print("="*78)
u=sp.Symbol('u',real=True)
up=(u+th)/(1+u*th); s_vec=ch+sh*u
check("Q-1a vector angle (two-sided): du'/du = s^-2  -> Q_c = Q_G = 2",
      rz(sp.diff(up,u)-1/s_vec**2))
T=sp.Symbol('T',positive=True)                      # tan phi
Tp=T/w                                              # G-3a with 2l = Lambda
dphi=sp.diff(sp.atan(Tp),T)/sp.diff(sp.atan(T),T)
q=ch+sh*(1-T**2)/(1+T**2)                            # s at the doubled angle
check("Q-1b spinor angle (one-sided): dphi'/dphi = q^-1 with q = s(Lambda, cos 2phi) -> Q_hbar = 1",
      rz(dphi-1/q), "same rotor, half representation, half exponent: the census (2,1,2)")

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
THEOREM (SECT-1): under RULE-1's bar pairing Cl(3) splits (4,4); its paravector
  sector carries the Minkowski norm t^2-|v|^2, invariant under the two-sided
  rotor. The three seats are three REGIMES of that one norm:
    c     null           (norm 0 on-state)          -> free
    hbar  timelike       (norm pinned to A0^2 > 0)  -> floored
    G     timelike->null (norm sech^2 -> 0 at wall)  -> horizon-vanishing
  and the character of each seat is FORCED by invariance once its regime is
  named. The hbar seat is the paravector read through its spinor square root
  (S-5): the quadratic-form/squeeze picture of thm_g is the boost of a
  timelike paravector (D-1..4), which DERIVES the 2:1 dictionary in both halves
  and the census (2,1,2) (Q-1). Only a sign is assigned per seat.
  Will's criterion: the seat whose behaviour WAS forced beyond its name is G
  (KIN-2 by import, THM-K by four declarations); the strain law un-forces it.
COMPARISON-STAGE NAMES: paravector model of Minkowski space in Cl(3); the
  form -> vector map is the Stokes / Poincare-sphere construction (double
  angle); Robertson-Schroedinger.
HONEST SCOPE: 'timelike' for the hbar state = C positive definite = SEAT-hb's
  content, a regime label; A0^2 = hbar^2/4 remains the label tier. Single axis
  K = e1 throughout.
""")
