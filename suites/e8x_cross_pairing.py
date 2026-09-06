#!/usr/bin/env python3
# =============================================================================
# E-8-X - THE CROSS-SEAT PAIRING  (2026-08-30; queued same day, TESTING-SCHEDULE)
#
# QUESTION (Will's ruling, from the superposition objection): E-8's invariant
# rides a pair with BOTH slots on one seat. Does a pairing with slots on
# DIFFERENT seats exist; is its phase the PRED-1 class; does it carry through
# the wall?
#
# ALLOWED INPUTS: BARE-1 (rotor); RULE-1 (metric-office pairing, Clifford bar);
#   SEAT-hb (hbar state is a quadratic form; thm_g part 1); G-8 dictionary AS
#   STRUCTURE: the hbar pivot is the SAME rotor in the one-sided (half)
#   representation [thm_g.py, result]; KIN-2 pinning tanh(Lambda)=sqrt(rs/r)
#   [thm_g2, declared datum]; PRED-1 banked B=(S+iV)/4, J=V [bargmann.py].
# DESIGN INPUT (named): the cross-seat co-pivot is ONE rotor acting two-sided
#   on the c slot and one-sided on the hbar slot. Only link in the frozen
#   record; a different intertwiner would be a new rule (fork c).
# COORDINATES: rational throughout - w = e^(Lambda/2), t = tan(Theta/2);
#   zero tests are polynomial-numerator checks, exact by construction.
# BANNED until comparison: Pancharatnam, Berry, solid angle by name.
#
# FORK: (a) no invariant beyond tautology -> the net fails
#       (b) invariant exists, phase = the class -> the net closes
#       (c) needs a new rule -> named debt
# =============================================================================
import sympy as sp
import mpmath as mp
import random, time
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
                key,sg=blade_mul(ka,kb)
                d[key]=d.get(key,0)+sg*va*vb
        return MV(d)
    __rmul__=__mul__
    def bar(s):
        SG={0:1,1:-1,2:-1,3:1}
        return MV({k:SG[len(k)]*v for k,v in s.d.items()})
    def rev(s):
        SG={0:1,1:1,2:-1,3:-1}
        return MV({k:SG[len(k)]*v for k,v in s.d.items()})
    def c(s,key): return s.d.get(tuple(key),0)
    def grades(s): return sorted({len(k) for k,v in s.d.items() if sp.simplify(v)!=0})
def rzero(e):
    z=sp.cancel(sp.together(sp.expand(e)))
    return sp.expand(sp.numer(z))==0
def mv_rzero(M): return all(rzero(v) for v in M.d.values())

w,t=sp.symbols('w t', positive=True)     # w=e^(Lambda/2), t=tan(Theta/2)
ch2,sh2=(w+1/w)/2,(w-1/w)/2
ct,st=(1-t**2)/(1+t**2), 2*t/(1+t**2)
ONE=MV({():1}); E1,E2=MV({(1,):1}),MV({(2,):1})
A  = MV({():ch2,(1,):sh2})               # the ONE rotor
X  = ONE + ct*E1 + st*E2                 # c slot: null ray
Pa=[sp.Symbol(f'a{i}',real=True) for i in range(8)]
Pb=[sp.Symbol(f'b{i}',real=True) for i in range(8)]
psa=MV({k:Pa[i] for i,k in enumerate(BL)})    # hbar slots, fully general
psb=MV({k:Pb[i] for i,k in enumerate(BL)})
Xp=A*X*A; psap=A*psa; psbp=A*psb

print("="*78); print("PART 1 - existence, necessity, where the phase lives (exact)"); print("="*78)
check("X-0  reversion fixes the rotor: rev(A) = A (scalar+vector)",
      mv_rzero(A.rev()-A))
Fa=psa*psa.rev(); Fap=psap*psap.rev()
check("X-1  the quadratic hbar object intertwines half -> full: F' = A F A",
      mv_rzero(Fap-(A*Fa*A)),
      "one-sided on psi becomes two-sided on psi psi~: G-8's 2:1 as rep theory")
lin=((Xp*psap).c(())-(X*psa).c(()))
sub={Pa[0]:1,Pa[1]:sp.Rational(1,3),Pa[2]:0,Pa[3]:0,Pa[4]:0,Pa[5]:0,Pa[6]:0,Pa[7]:0,
     t:sp.Rational(1,2),w:sp.Rational(3,2)}
check("X-2  NO-GO: <X' psi'> != <X psi> generically - a LINEAR hbar slot fails",
      sp.simplify(lin.subs(sub))!=0,
      "the net cannot close on a one-sided slot; the quadratic is REQUIRED")
PS =(X*Fa.bar()).c(());  PI =(X*Fa.bar()).c((1,2,3))
PSp=(Xp*Fap.bar()).c(()); PIp=(Xp*Fap.bar()).c((1,2,3))
check("X-3  CROSS-SEAT INVARIANT exists: <X' F'bar>_S = <X F bar>_S, all w,t,psi",
      rzero(PSp-PS),
      "one rotor, two representations, one invariant - the net closes")
check("X-4  STATE LEMMA: psi psi~ has grades {0,1} only => <X Fbar>_I = 0",
      Fa.grades()==[0,1] and rzero(PI),
      "the hbar STATE pairs REAL: a single-seat reading carries no phase (PRED-1 echo)")
sub4={Pa[0]:1,Pa[1]:sp.Rational(1,2),Pa[2]:0,Pa[3]:0,Pa[4]:sp.Rational(1,2),
      Pa[5]:0,Pa[6]:0,Pa[7]:sp.Rational(1,3),w:1,t:sp.Rational(1,3)}
check("X-5  non-tautological: dPS/dt != 0 and dPS/dpsi != 0 at a generic sample",
      sp.simplify(sp.diff(PS,t).subs(sub4))!=0 and
      any(sp.simplify(sp.diff(PS,Pa[i]).subs(sub4))!=0 for i in range(8)),
      "content in BOTH slots; E-8 had both slots on one seat")
# the TRANSITION bilinear: a leg between two hbar states
Gab=psa*psb.rev(); Gabp=psap*psbp.rev()
QS =(X*Gab.bar()).c(());  QI =(X*Gab.bar()).c((1,2,3))
QSp=(Xp*Gabp.bar()).c(()); QIp=(Xp*Gabp.bar()).c((1,2,3))
check("X-6a the TRANSITION pairing is invariant too: S-part, all w,t,psi_a,psi_b",
      rzero(QSp-QS))
check("X-6b ... and its I-part: <X' (psi_a psi_b~)'bar>_I = <X (psi_a psi_b~)bar>_I",
      rzero(QIp-QI),
      "legs between hbar states pair COMPLEX: the phase lives on TRANSITIONS")
sub6={Pa[0]:1,Pa[1]:0,Pa[2]:sp.Rational(1,2),Pa[3]:0,Pa[4]:sp.Rational(1,3),Pa[5]:0,Pa[6]:0,Pa[7]:0,
      Pb[0]:1,Pb[1]:sp.Rational(1,5),Pb[2]:0,Pb[3]:sp.Rational(1,2),Pb[4]:0,Pb[5]:sp.Rational(1,7),Pb[6]:0,Pb[7]:0,
      w:sp.Rational(4,3),t:sp.Rational(2,5)}
check("X-6c the I-part is generically NONZERO on a transition",
      sp.simplify(QI.subs(sub6))!=0,
      "states read real, transitions carry phase: the split is structural")

print(); print("="*78); print("PART 2 - the wall (KIN-2 pinning; H-3b pattern cross-seat)"); print("="*78)
r,rs,Lm=sp.symbols('r r_s Lambda', positive=True)
vv=sp.sqrt(rs/r)
check("W-1  pinning translation: cosh(Lam(r)) = 1/sqrt(1-rs/r) -> oo at the wall",
      sp.simplify(sp.cosh(sp.atanh(vv)) - 1/sp.sqrt(1-rs/r))==0 and
      sp.limit(1/sp.sqrt(1-rs/r), r, rs, '+')==sp.oo,
      "r -> rs is Lambda -> oo is w -> oo: the c slot's reading diverges")
tprime=(w**2+w**-2)/2 + (w**2-w**-2)/2*ct
check("W-2  the c slot diverges: t'(w) -> oo as w -> oo; the pairing has dPS'/dw = 0",
      sp.limit(tprime.subs(t,sp.Rational(1,3)), w, sp.oo)==sp.oo and rzero(sp.diff(PSp,w)),
      "slots die/diverge at the wall; the cross-seat pair walks through - H-3b")

print(); print("="*78); print("PART 3 - numerics: the class, and the polarisation reading"); print("="*78)
mp.mp.dps=30; random.seed(11)
def mvec():
    v=[mp.mpf(random.gauss(0,1)) for _ in range(3)]
    n=mp.sqrt(sum(q*q for q in v)); return [q/n for q in v]
def mspin(a):
    th=mp.acos(a[2]); ph=mp.atan2(a[1],a[0])
    return [mp.cos(th/2), mp.e**(1j*ph)*mp.sin(th/2)]
def ovl(x,y): return mp.conj(x[0])*y[0]+mp.conj(x[1])*y[1]
def trip(a,b,c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])
            +a[2]*(b[0]*c[1]-b[1]*c[0]))
ok1=True
for _ in range(5):
    Avs=[mvec() for _ in range(3)]
    g=lambda i,j: sum(x*y for x,y in zip(Avs[i],Avs[j]))
    B=(1+g(0,1)+g(1,2)+g(0,2)+1j*trip(*Avs))/4
    a_,b_,c_=(mp.acos(g(1,2)),mp.acos(g(0,2)),mp.acos(g(0,1)))
    s_=(a_+b_+c_)/2
    E=4*mp.atan(mp.sqrt(mp.tan(s_/2)*mp.tan((s_-a_)/2)*mp.tan((s_-b_)/2)*mp.tan((s_-c_)/2)))
    ok1=ok1 and abs(2*abs(mp.arg(B))-E)<mp.mpf(10)**-24
check("Y-1  2|arg B| = spherical excess (l'Huilier, independent route), 5 triples",
      ok1, "the class IS half the area the seat cycle encloses - J's geometric name")
ok2=True
for eps in ('0.1','0.01','0.001'):
    a1=mvec(); e=mp.mpf(eps); td=mvec()
    a2=[(-x+e*y) for x,y in zip(a1,td)]
    n=mp.sqrt(sum(x*x for x in a2)); a2=[x/n for x in a2]
    gam=sum(x*y for x,y in zip(a1,a2))
    leg2=abs(ovl(mspin(a1),mspin(a2)))**2
    ok2=ok2 and abs(leg2-(1+gam)/2)<mp.mpf(10)**-24
check("Y-2  |<psi_a|psi_b>|^2 = (1+gamma)/2 exactly: crossed (antipodal) kills the",
      ok2, "LEG through the Kummer radicand 1+gamma -> 0; 90deg polarisation = 180deg sphere")
lamv=mp.mpf('0.83'); phv=mp.mpf('0.37')
q=lambda l,f: mp.cosh(2*l)+mp.sinh(2*l)*mp.cos(2*f)
ok3=(abs(q(lamv,phv+mp.pi/2)-q(-lamv,phv))<mp.mpf(10)**-25 and
     abs(mp.sqrt(q(lamv,0)*q(lamv,mp.pi/2))-1)<mp.mpf(10)**-25)
check("Y-3  turn the hbar reading 90deg: the pivot INVERTS (q -> q at -lambda);",
      ok3, "conjugate product = floor: the Malus structure was frozen as G-5/G-6")
Btr=(1+3*(-mp.mpf(1)/2))/4
check("Y-4  trine: B = -1/8 real, arg B = pi - PRED-1's 'reads pi on the far class'",
      abs(mp.arg(Btr)-mp.pi)<mp.mpf(10)**-25)

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
VERDICT (fork): (b). THE NET CLOSES - with a sharper split than designed.
  EXISTENCE: <X' F'bar>_S is invariant for one rotor acting two-sided on the
  c slot and one-sided on the hbar slot (X-3), non-tautological (X-5). It
  exists BECAUSE the hbar slot is quadratic (SEAT-hb): psi psi~ intertwines
  the half representation into the full one (X-1); a linear hbar slot
  provably fails (X-2). G-8's 2:1 dictionary is the vector/spinor split of
  ONE rotor.
  WHERE THE PHASE LIVES: the hbar STATE pairs real - psi psi~ has grades
  {0,1} and the I-part vanishes identically (X-4). The phase enters only
  through TRANSITIONS psi_a psi_b~ between two hbar states, whose pairing is
  invariant in BOTH parts and generically complex (X-6). States read real;
  transitions carry phase. This is tonight's PRED-1 verdict (modulus =
  gauge, phase = class) rederived from the algebra with no cocycle language.
  THE WALL: on the KIN-2-pinned family the c slot diverges at r = rs while
  the pairing is w-independent (W-1/2): the pair walks through the wall -
  H-3b's pattern with the slots on DIFFERENT seats.
  THE CLASS: cycle values are the banked B with 2|arg B| = the spherical
  excess (Y-1): J = V is the area-tangent identity; the medium's phase is
  half the area the seat cycle encloses.
  POLARISATION READING (Will's instinct, receipts): the hbar seat carries
  the cos-2phi Malus structure and the 90-degree turn inverts the pivot
  (Y-3 = frozen G-5/G-6); crossed = antipodal kills the LEG exactly through
  the Kummer radicand 1+gamma (Y-2) while the class survives off the wall;
  the wall itself reads pi (Y-4).
COMPARISON-STAGE NAMES (spoken here only): Pancharatnam 1956 (polarisation-
  cycle phase = half solid angle); Van Oosterom-Strackee 1983 (the area-
  tangent formula tan(Omega/2) = V/(1+Sum gamma) - PRED-1's J identity
  verbatim). MACHINERY INHERITED: polarimetric interferometry as the
  protocol family for PRED-1; geodesic-triangle area calculus for J.
HONEST SCOPE: the co-pivot intertwiner (one rotor, two reps) is the design
  input, licensed by G-8 but NAMED; single axis K = e1 (F-4's fence); the
  wall uses KIN-2 as declared datum; lineage above is comparison-stage
  naming from recall, not a database search tonight.
""")
