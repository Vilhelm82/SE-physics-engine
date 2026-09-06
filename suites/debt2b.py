#!/usr/bin/env python3
# =============================================================================
# DEBT-2b - c(omega) FROM THE CHANNEL LAW  (2026-08-30)
#
# THE DEBT (RULING-1 / debt2.py): the bath coupling was DECLARED Ohmic-in-band,
# not derived. Discharge condition: derive c(omega) from F's channel law.
#
# ALLOWED INPUTS (frozen, receipts):
#   FAM     presentation family s(u') = 1/(gamma(1-beta u')), seat measure
#           du'/2                                          [E, thm_f]
#   M-2     du'/du = 1/s^2 - the office Jacobian           [thm_e_addendum2]
#   RULE-4  W = 1 on the SEAT measure: equal per-direction weight at every
#           gap                                            [thm_f, named rule]
#   LBL-2   per-direction Bose occupation n(E/T(u'))       [thm_f, declared]
#   GR      two-level golden rule, position coupling: per-mode weight
#           c(w0)^2/(2 w0), resonance at the PRESENTED frequency. Standard
#           import, debt2.py's FKM conventions - NAMED, same epistemic shape
#           as LBL-2.
#   F-16    <s^k> = sinh((k-1)L)/((k-1) sinh L)            [thm_f]
# BANNED until comparison: 'Ohmic', Caldeira-Leggett, spectral-density names.
#
# FORK: (a) c^2 g0 forced uniquely -> DEBT-2b discharged
#       (b) underdetermined -> debt stands, freedom named
#       (c) forced to something else -> debt2.py amended
# TWO PATHS: exact sympy (functional equation + closed-form integrals);
#            independent numeric mode-sum with the aberration map.
# =============================================================================
import sympy as sp
import mpmath as mp
import random, time
T0=time.time(); CH=[]
def check(t_,ok,n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t_}"+(f" - {n}" if n else ""), flush=True)

L  = sp.Symbol('Lambda', positive=True)
up = sp.Symbol('up', real=True)
E, y, a = sp.symbols('E y alpha', positive=True)
b, g = sp.tanh(L), sp.cosh(L)
s = 1/(g*(1-b*up))                      # the family, seat office (F's form)

print("="*78); print("PART 1 - the forced coupling (exact)"); print("="*78)

# B-1 the rate integrand, offices handled: modes counted on the source plane
# (du/2), channel weight demanded uniform on the seat plane (du'/2).
# du = s^2 du' (M-2 inverted); resonance delta eats dw0 giving 1/s; golden-rule
# amplitude^2 gives 1/(2 w0) = s/(2E). Per-(du'/2) weight at gap E:
#     W(u',E) = s^2 * F(E/s) * (s/E) * (1/s)  ->  (1/E) s^2 F(E/s)
# with F(w0) := c(w0)^2 g0(w0). RULE-4: W independent of u'.
F = sp.Function('F', positive=True)
Wgt = s**2 * F(E/s)
dW = sp.simplify(sp.diff(Wgt, up))
sprime = sp.diff(s, up)
check("B-1  d/du'[s^2 F(E/s)] = s'(u') * s * [2F(y) - y F'(y)],  y = E/s: EXACT",
      sp.simplify(sp.expand(dW - sprime*s*(2*F(E/s) - (E/s)*sp.diff(F(y),y).subs(y,E/s)))) == 0,
      "s' != 0 for Lambda > 0, so RULE-4 uniformity forces 2F = yF' pointwise")
sol = sp.dsolve(sp.Eq(2*F(y) - y*sp.diff(F(y),y), 0), F(y))
check("B-2  the ODE 2F = yF' has general solution F = C*y^2 - UNIQUE up to scale",
      sp.simplify(sp.diff(sp.simplify(sol.rhs/y**2), y)) == 0,
      "c(w)^2 g0(w) = C w^2: the coupling's radial law is FORCED, eta the only knob")

# B-3 sufficiency: with F = C w^2 the per-direction weight is u'-INDEPENDENT
C0 = sp.Symbol('C', positive=True)
Wfor = sp.simplify((s**2 * C0*(E/s)**2))
check("B-3  sufficiency: F = C w^2  =>  weight = C E^2, uniform in u', every gap",
      sp.simplify(Wfor - C0*E**2) == 0,
      "and the rate envelope is phi(E) = C E^2 * (1/E)-amplitude = C E: fixed too")

# B-4 the failed alternatives, recorded (where the wrong answers live):
for aa,nm in ((0,"flat c^2 g0"),(1,"F ~ w (one-office ghost)")):
    Wbad = sp.simplify(s**2*(E/s)**aa)
    check(f"B-4{'ab'[aa]} {nm}: weight = E^{aa} s^{2-aa} - u'-DEPENDENT, violates RULE-4",
          sp.simplify(sp.diff(Wbad,up)) != 0,
          "office mixing shifts the exponent by exactly the Jacobian power Q=2" if aa==1 else "")

# B-5 pivot eigenstructure on power-law weights: eigenvalue(alpha) = <s^alpha>
mom = lambda k: sp.simplify(sp.integrate(s**k, (up,-1,1))/2)
closed = sp.sinh((a-1)*L)/((a-1)*sp.sinh(L))
check("B-5a eigenvalue(alpha) = <s^(2-alpha)> = <s^alpha> (the k <-> 2-k mirror, C-3)",
      sp.simplify(sp.expand((sp.sinh((1-a)*L)/((1-a)*sp.sinh(L)) - closed).rewrite(sp.exp))) == 0,
      "band route and direction route agree BY the mirror: Renyi skew, third sighting")
tab = {1: L/sp.sinh(L), 2: sp.Integer(1), 3: sp.cosh(L), 4: sp.cosh(L)**2*(1+sp.tanh(L)**2/3)}
ok5 = all(sp.simplify(sp.expand((mom(k)-v).rewrite(sp.exp)))==0 for k,v in tab.items())
check("B-5b eigenvalue table = the F-12..15 moment table verbatim (alpha = 1,2,3,4)",
      ok5, "L/sinhL, 1, gamma, gamma^2(1+beta^2/3): the moments ARE the pivot eigenvalues")
# monotonicity lemma: h(k) = sinh(kL)/k is even and strictly increasing for k>0,
# since d/dk[k cosh k - sinh k] = k sinh k > 0 and it vanishes at 0. Hence
# <s^alpha> = 1 has EXACTLY the solutions |alpha-1| = 1.
xk = sp.Symbol('x', positive=True)
gmono = xk*sp.cosh(xk) - sp.sinh(xk)
check("B-6  fixed-point lemma: d/dx[x cosh x - sinh x] = x sinh x > 0, g(0) = 0",
      sp.simplify(sp.diff(gmono, xk) - xk*sp.sinh(xk)) == 0 and gmono.subs(xk,0) == 0,
      "=> sinh(kL)/k strictly increasing in |k| => <s^a>=1 iff |a-1|=1 iff a in {0,2}")
ok6 = (sp.simplify(mom(0)-1)==0 and sp.simplify(sp.expand((mom(2)-1).rewrite(sp.exp)))==0
       and sp.simplify(sp.expand((mom(1)-1).rewrite(sp.exp)))!=0
       and sp.simplify(sp.expand((mom(3)-1).rewrite(sp.exp)))!=0)
check("B-6b <s^0>=<s^2>=1 exactly; <s^1>,<s^3> != 1: exactly two blind exponents",
      ok6, "the invariant's two access routes (F part 1), now as bath fixed points")

# B-7 debt2.py's declared discretisation IS the derived law:
D, eta, w = sp.symbols('D eta w', positive=True)
ck2_g0 = (2/sp.pi*eta*D)*w**2 * (1/D)          # ck = sqrt(2/pi eta D) w, g0 = 1/D
check("B-7  debt2's ck = sqrt(2 eta D/pi) w with flat g0: c^2 g0 = (2 eta/pi) w^2",
      sp.simplify(ck2_g0 - 2*eta*w**2/sp.pi) == 0,
      "DECLARED there, DERIVED here: debt2.py needs zero amendment")

print(); print("="*78); print("PART 2 - independent numeric mode-sum (aberration map)"); print("="*78)
mp.mp.dps=25
worstU=mp.mpf(0); worstA=mp.mpf(0)
import math
for (xv,lv) in ((0.7,0.9),(2.5,0.4),(5.0,1.3)):
    # DETERMINISTIC grid on the SOURCE plane (du/2), binned by SEAT angle u'.
    # Sampling du and binning in u' supplies du/du' = s^2 by itself; the
    # per-mode weight must therefore be carried WITHOUT it (the B-4b lesson,
    # committed numerically on the first attempt and corrected here):
    #   per-mode weight at gap E:  F(E/s) * (1/(2 w0)) * (1/s) = C E/(2 s^2)
    Nu=200000; nb=20; bins=[0.0]*nb
    th=math.tanh(lv)
    for i in range(Nu):
        uu=-1.0+(2.0*(i+0.5))/Nu
        ss=math.cosh(lv)+math.sinh(lv)*uu
        upv=(uu+th)/(1.0+th*uu)
        wgt=xv/(2.0*ss*ss)
        k=min(nb-1,int((upv+1.0)*nb/2)); bins[k]+=wgt
    dens=[bins[k]*nb/(Nu) for k in range(nb)]
    m=sum(dens)/nb
    worstU=max(worstU, max(abs(d-m)/m for d in dens))
    Aq=mp.quad(lambda uv: 1/(mp.e**(mp.mpf(xv)*mp.cosh(mp.mpf(lv))*(1-mp.tanh(mp.mpf(lv))*uv))-1),[-1,1])/2
    Acl=mp.log((1-mp.e**(-mp.mpf(xv)*mp.e**mp.mpf(lv)))/(1-mp.e**(-mp.mpf(xv)*mp.e**(-mp.mpf(lv)))))/(2*mp.mpf(xv)*mp.sinh(mp.mpf(lv)))
    worstA=max(worstA,abs(Aq-Acl))
check("N-1  seat-office uniformity of the derived weight: flat in u' across 20 bins",
      worstU < 0.01, f"worst relative deviation {float(worstU):.2e} (grid, 2e5 rays x3)")
check("N-2  channel content intact: quadrature A = F-4 closed form at 3 (x,L) points",
      worstA < mp.mpf(10)**-20, f"worst {mp.nstr(worstA,3)}")

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
VERDICT (fork): (a) with a bonus. DEBT-2b DISCHARGED:
  RULE-4 (equal weight per SEAT direction, every gap) + M-2 (the office
  Jacobian du = s^2 du') + the golden rule force the functional equation
  s^2 F(E/s) = psi(E), whose only solutions are F(w) = c(w)^2 g0(w) = C w^2
  (B-1/2/3). The coupling's radial law is NOT a knob; eta (overall scale) is
  the only freedom. debt2.py's 'Ohmic-in-band, declared' becomes DERIVED with
  zero amendment (B-7).
  THE FIXED-POINT READING: the pivot acts on power-law bath weights with
  eigenvalue <s^alpha> - the F-12..15 moment table IS the eigenvalue table
  (B-5) - and the forced alpha = 2 is the unique nontrivial fixed point,
  i.e. the Jacobian-blind route of F-13 (B-6). The channel law selects the
  pivot-invariant bath: the sky every seat derives is the same sky.
  THE GHOST, recorded (B-4): demanding uniformity on the SOURCE measure
  instead yields F ~ w - one office error shifts the exponent by exactly the
  Jacobian power Q = 2. Same genus as addendum 1 and PRED-1's DeltaQ ghost;
  third specimen, now with the shift quantified.
PAYMENTS: DEBT-2b closed. Pressure-test 1 half-paid: the coupling exponent
  is forced (alpha = 2 = the Jacobian power) - a distinguished Mellin
  exponent exists and the pivot acts as a Mellin multiplier with symbol
  <s^alpha>; Q-RN adjacency strengthened, not closed.
HONEST SCOPE: Lambda > 0 required (no pivot, no constraint); eta free by
  scale; GR amplitude convention imported (named, LBL-2's epistemic shape);
  isotropic c (omega-only) is RULE-4's own isotropy; two-level probe is F's
  probe - the oscillator GLE inherits the same weight by D2-3's shared
  coefficients, not re-verified here.
COMPARISON STAGE (names spoken here only): F = C w^2 with flat band is
  J(w) = eta w - OHMIC, Caldeira-Leggett's standard case; the E-envelope
  phi(E) = C E is the Ohmic spontaneous+stimulated scaling; 'the moments are
  the eigenvalues' is Mellin diagonalisation of a multiplicative convolution.
""")
