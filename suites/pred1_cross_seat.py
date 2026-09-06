#!/usr/bin/env python3
# =============================================================================
# PREDICTION-1, the cross-seat cycle test  (2026-08-30)
#
# KILL CONDITION (HANDOFF-2026-08-28, item 3): "if the cross-seat cycle is a
# coboundary, the whole framing is bookkeeping."
#
# ALLOWED INPUTS (all frozen, receipts named):
#   BARE-1   pivot = hyperbolic rotor                      [D1]
#   C-1      sigma_c(L;u) = cosh L + sinh L * u is a multiplicative cocycle
#            on the c boost groupoid                       [graph_cocycle.py]
#   G-8      cover dictionary Theta = 2 phi, Lambda = 2 lambda; the hbar
#            squeeze IS the c boost through the 2:1 cover  [thm_g.py, result]
#   Q-CENSUS transport exponents: Q_c = 2 (sphere, M-2), Q_hbar = 1 (circle,
#            G-3b), Q_G = 2 (shared sphere, thm_g2 uncounted note)
#   N-2      s = cosh(y+L)/cosh(y) in rapidity             [thm_e_addendum3]
#   BP       seat-cycle Bargmann holonomy B=(S+iV)/4       [bargmann.py, banked]
# BANNED: nothing physical is needed; this is structure only.
#
# PRE-REGISTERED FORK:
#   (a) every layer of the cross-seat cycle is a coboundary -> KILL fires
#   (b) some layer carries a non-trivial class -> PREDICTION-1 survives, class named
#   (c) the cycle is ill-defined without a new rule -> named debt
#
# TWO PATHS sharing only the CAS idea, not code: Path 1 exact sympy on the
# positive (modulus) layer; Path 2 mpmath numerics on maps + measures + the
# spinor phase layer, no sympy objects reused.
# =============================================================================
import sympy as sp
import mpmath as mp
import random, time
T0=time.time(); CH=[]
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}"+(f" - {n}" if n else ""), flush=True)

u  = sp.Symbol('u', real=True)
L, l, L1 = sp.symbols('Lambda ell Lambda1', positive=True)
boost = lambda R, x: (x + sp.tanh(R))/(1 + sp.tanh(R)*x)
sig   = lambda R, x: sp.cosh(R) + sp.sinh(R)*x

print("="*78); print("PART 1 - the modulus layer, exact"); print("="*78)

# P-0/1/2: work in rapidity, the record's own coordinate (N-1/N-2, addendum 3).
y0 = sp.Symbol('y0', real=True)
check("P-0  boost is translation in rapidity: boost(L, tanh y) = tanh(y+L)",
      sp.simplify((boost(L, sp.tanh(y0)) - sp.tanh(y0+L)).rewrite(sp.exp))==0,
      "N-1 re-verified; licenses the substitution used below")
# P-1: the presentation cocycle is ALREADY a coboundary on one seat.
lhs1 = sp.cosh(y0+L)/sp.cosh(y0)          # b(gu)/b(u), b = cosh(y) = (1-u^2)^(-1/2)
check("P-1  sigma_c(L;u) = b(g u)/b(u) with b = (1-u^2)^(-1/2), EXACT",
      sp.simplify((lhs1 - sig(L, sp.tanh(y0))).rewrite(sp.exp))==0,
      "N-2 restated: s = cosh(y+L)/cosh y. The c cocycle is a coboundary, always was")

# P-2: same certificate on the hbar seat through the cover, in T = tan(phi).
T = sp.Symbol('T', positive=True)
sin2f  = 2*T/(1+T**2)
Tp     = sp.exp(-2*l)*T                    # G-3a squeeze
sin2fp = 2*Tp/(1+Tp**2)
qT     = sp.cosh(2*l) + sp.sinh(2*l)*(1-T**2)/(1+T**2)
check("P-2  q(l;phi) = sin(2 phi)/sin(2 phi') : the hbar cocycle is the SAME",
      sp.simplify(sp.together(sin2f/sin2fp - qT).rewrite(sp.exp))==0,
      "coboundary pulled through the 2:1 cover - one certificate, two seats")

# P-3: the elementary cross-seat cycle closes on points.
u0 = sp.Symbol('u0', real=True)
u1 = boost(L, u0); u2 = boost(2*l, u1); u3 = boost(-(L+2*l), u2)
check("P-3  CLOSURE: boost(-(L+2l)) . boost(2l) . boost(L) = id on the point",
      sp.simplify(sp.together(u3 - u0).rewrite(sp.exp))==0,
      "c-pivot L, descend, squeeze l, ascend, c-pivot -(L+2l): the loop closes")

# P-4: presentation holonomy around the closed loop = 1 identically.
HP = sp.simplify(sig(-(L+2*l), u2)*sig(2*l, u1)*sig(L, u0))
check("P-4  presentation holonomy H_P = 1 IDENTICALLY (all L, l, u0)",
      sp.simplify(sp.expand_trig(sp.together(HP-1).rewrite(sp.exp)))==0,
      "forced by P-1: a coboundary telescopes around every closed loop")

# P-5: the ALL-IN measure holonomy (seat measures + honest cover densities),
#   in rapidity: u0 = tanh(ya), u1 = tanh(yb), u2 = tanh(yc), yb = ya+L, yc = yb+2l;
#   sqrt(1-u^2) = 1/cosh(y). Factors: c-legs sigma^2, hbar leg sigma(2l)^1,
#   descend pi/cosh(yb), ascend 2*cosh(yc)/pi.
ya = sp.Symbol('ya', real=True); yb = ya + L; yc = yb + 2*l
H_all = ( (sp.cosh(ya)/sp.cosh(yc))**2      # c-leg back: sigma(-(L+2l); u2)^2
        * (2*sp.cosh(yc)/sp.pi)             # ascend at phi2
        * (sp.cosh(yc)/sp.cosh(yb))         # hbar leg (Q=1)
        * (sp.pi/sp.cosh(yb))               # descend at phi1
        * (sp.cosh(yb)/sp.cosh(ya))**2 )    # c-leg out: sigma(L; u0)^2
check("P-5  ALL-IN modulus holonomy = 2 (the cover degree), CONSTANT in (L,l,u0)",
      sp.simplify((H_all-2).rewrite(sp.exp))==0,
      "no rapidity or base-point dependence survives the closed loop")
# the constant 2 is the degree of the cover, a topological constant, itself a
# coboundary on the two-object graph (absorb into b at the hbar objects).

# P-6: the within-seat truncation reproduces the naive 'nontrivial' answer -
# recorded to show WHERE a wrong verdict would have come from.
H_within = sig(-(L+2*l),u2)**2 * sig(2*l,u1)**1 * sig(L,u0)**2
gap = sp.simplify(sp.expand_trig(sp.together(H_within - sig(2*l,u1)**(-1)).rewrite(sp.exp)))
check("P-6  transitions dropped => H_within = sigma(2l;u1)^(-1): the DeltaQ ghost",
      gap==0,
      "exponent deficit Q_c-Q_hbar=1 on the foreign leg; PAID by the cover density")

print(); print("="*78); print("PART 2 - independent numerics + the phase layer"); print("="*78)
mp.mp.dps=30; random.seed(7)
def nboost(R,x): return (x+mp.tanh(R))/(1+mp.tanh(R)*x)
def nsig(R,x):  return mp.cosh(R)+mp.sinh(R)*x
worst=mp.mpf(0)
for _ in range(4):
    Lv=mp.mpf(random.uniform(0.2,1.7)); lv=mp.mpf(random.uniform(0.2,1.4))
    x0=mp.mpf(random.uniform(-0.9,0.9)); br=random.choice([0,1])
    # transport an interval, all-in bookkeeping, finite differences
    eps=mp.mpf(10)**-12
    def cycle(x):
        y1=nboost(Lv,x)
        ph1=mp.acos(y1)/2 + br*mp.pi/2          # both branches of the cover
        ph2=mp.atan(mp.e**(-2*lv)*mp.tan(ph1)) if abs(mp.cos(ph1))>1e-9 else ph1
        y2=mp.cos(2*ph2)
        return nboost(-(Lv+2*lv),y2), ph1, ph2
    e0,ph1,ph2=cycle(x0); e1,_,_=cycle(x0+eps)
    check_id = abs(e0-x0)
    # all-in factor product at the running points
    y1=nboost(Lv,x0); y2=mp.cos(2*ph2)
    Hnum=( nsig(-(Lv+2*lv),y2)**2 * (2/(mp.pi*mp.sqrt(1-y2**2)))
         * nsig(2*lv,y1) * (mp.pi*mp.sqrt(1-y1**2)) * nsig(Lv,x0)**2 )
    worst=max(worst, abs(Hnum-2), check_id)
check("N-1  numeric all-in holonomy = 2 on random loops, BOTH cover branches",
      worst<mp.mpf(10)**-24, f"worst deviation {mp.nstr(worst,3)}")

# N-2: the phase layer. Same geometric object (three seats/rays), spinor legs.
def mvec():
    v=[mp.mpf(random.gauss(0,1)) for _ in range(3)]
    n=mp.sqrt(sum(t*t for t in v)); return [t/n for t in v]
def mspin(a):
    th=mp.acos(a[2]); ph=mp.atan2(a[1],a[0])
    return [mp.cos(th/2), mp.e**(1j*ph)*mp.sin(th/2)]
def ovl(x,y): return mp.conj(x[0])*y[0]+mp.conj(x[1])*y[1]
def trip(a,b,c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])
            +a[2]*(b[0]*c[1]-b[1]*c[0]))
okph=okreph=okrev=okmod=True
for _ in range(5):
    Avs=[mvec() for _ in range(3)]
    Ss=[mspin(a) for a in Avs]
    B = ovl(Ss[0],Ss[1])*ovl(Ss[1],Ss[2])*ovl(Ss[2],Ss[0])
    S_=1+sum(sum(x*y for x,y in zip(Avs[i],Avs[j])) for i,j in ((0,1),(1,2),(0,2)))
    V_=trip(*Avs)
    okph  = okph  and abs(B-(S_+1j*V_)/4)<mp.mpf(10)**-25       # banked closed form
    ph=[mp.e**(1j*mp.mpf(random.uniform(0,6.28))) for _ in range(3)]
    Sg=[[p*c for c in s0] for p,s0 in zip(ph,Ss)]
    Bg=ovl(Sg[0],Sg[1])*ovl(Sg[1],Sg[2])*ovl(Sg[2],Sg[0])
    okreph= okreph and abs(B-Bg)<mp.mpf(10)**-25                 # rephasing = coboundary
    Sr=[mspin([-t for t in a]) for a in Avs]
    Br=ovl(Sr[0],Sr[1])*ovl(Sr[1],Sr[2])*ovl(Sr[2],Sr[0])
    okrev = okrev and abs(Br-mp.conj(B))<mp.mpf(10)**-25         # Z2 conjugates
    okmod = okmod and abs(mp.arg(B))>mp.mpf(10)**-6              # phase generically != 0
check("N-2a phase layer: B = (S+iV)/4 on the seat cycle (banked form re-verified)", okph)
check("N-2b vertex rephasings leave B invariant: positive-gauge CANNOT move it", okreph)
check("N-2c global reversal conjugates B: the one Z_2, click-monodromy shape", okrev)
check("N-2d arg B != 0 generically: the class is NON-TRIVIAL and equals arg(S+iV)", okmod)

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
VERDICT (pre-registered fork):
  MODULUS LAYER -> (a). The cross-seat cycle's positive layer is a coboundary
  - explicitly: b = (1-u^2)^(-1/2) per seat, cover density as the transition
  b, loop value = the cover degree, a constant. The sigma-framing of the
  cross-seat cycle IS bookkeeping. Kill condition FIRES on this layer.
  The DeltaQ 'holonomy' sigma^(Q_c-Q_hbar) is real as a DECOMPOSITION term
  (P-6) but it is paid exactly by the cover measure density (P-5): the
  exponent census (2,1,2) is a statement about how the books are kept, not
  about a class.
  PHASE LAYER -> (b). The same closed cycle carries arg B = arg(S+iV):
  invariant under every positive/vertex gauge (N-2b), conjugated by the one
  Z_2 (N-2c), generically nonzero (N-2d), and equal to the banked
  PREDICTION-1 observable (J = V at s = 1/2). The cross-seat holonomy CLASS
  exists and it is exactly the object PREDICTION-1 already names.
  NET: split verdict (a)+(b). What a seat can read alone (positive factors)
  is gauge; what survives the loop is orientation/phase - V, the sand.
NAMED SCOPE: transition map = the G-8 cover dictionary (the only seat
  transition in the frozen record); G-legs enter with Q=2 = c and add no
  modulus content; the spin-lift branch Z_2 (click paper) is consistent with
  N-2c but the fixed-locus comparison owed by G-8 remains owed.
QUEUED NEXT (Will, 2026-08-30): E-8-X, the cross-seat pairing - construct
  <X'_c, b'_hbar>-type pairings with slots on DIFFERENT seats; question:
  does a pairing invariant exist cross-seat, and is its phase this class?
""")
