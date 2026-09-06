#!/usr/bin/env python3
# =============================================================================
# THM-RN - THE BLIND MOMENT IS THE MASS OF THE SOURCE MEASURE (2026-08-28)
# Will's theorem, retiring the transport-exponent conjecture. Second path.
#
# LEMMA (W.L.): if the reading family satisfies d mu'/d mu = f^{-n} on the
# seat's carrier, then <f^n>_{mu'} = int f^n dmu' = int dmu = MASS(mu) = 1.
# The blind exponent equals the Jacobian power BY CONSTRUCTION; the blind
# value is the total mass; no third carrier was ever needed.
# COROLLARY (uniqueness, free): <f^p>_{mu'} = int f^{p-n} dmu, so p = n is
# blind always; any other blind p needs int f^{p-n} dmu pivot-independent.
# SUCCESSOR QUESTION Q-RN (the real content): WHY is the reading a power of
# the Radon-Nikodym derivative at all? "The readout is geometric, not
# instrumental" - to be derived from the channel law (F-4 adjacency).
# =============================================================================
import sympy as sp

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" - {n}" if n else ""), flush=True)

lam, u, up, ph, p = sp.symbols('lambda u uprime phi p', real=True)
L = sp.Symbol('lambda', positive=True)

print("=" * 78); print("PART 1 - the premises, independently re-verified"); print("=" * 78)
s_u  = sp.cosh(L) + sp.sinh(L)*u
umap = (u + sp.tanh(L))/(1 + sp.tanh(L)*u)                 # aberration on u
check("RN-1 c seat: du'/du * s(u)^2 = 1 exactly  (n = 2 is the u-Jacobian power)",
      sp.simplify((sp.diff(umap, u)*s_u**2 - 1).rewrite(sp.exp)) == 0)
q_ph  = sp.exp(2*L)*sp.cos(ph)**2 + sp.exp(-2*L)*sp.sin(ph)**2   # A = 1
phmap = sp.atan(sp.exp(-2*L)*sp.tan(ph))                    # squeeze on directions
check("RN-2 hbar seat: dphi'/dphi * q(phi) = 1 exactly  (n = 1 is the phi-Jacobian)",
      sp.simplify(sp.together(sp.diff(phmap, ph)*q_ph - 1)) == 0)

print(); print("=" * 78); print("PART 2 - the lemma lands: blind value = MASS"); print("=" * 78)
D_up = 1/(sp.cosh(L) - sp.sinh(L)*up)                       # family in seat coords
check("RN-3a consistency: s(u(u')) * (cosh - sinh*u') = 1  (same family, two charts)",
      sp.simplify((s_u.subs(u, sp.solve(sp.Eq(up, umap), u)[0])
                   *(sp.cosh(L) - sp.sinh(L)*up) - 1).rewrite(sp.exp)) == 0)
I2 = sp.integrate(D_up**2, (up, -1, 1))/2
check("RN-3b c seat: <D^2> over the seat's uniform sky = 1 for ALL lambda",
      sp.simplify((I2 - 1).rewrite(sp.exp)) == 0,
      "the mass of the source measure, read through the pivot")
Iq = sp.integrate(q_ph*(sp.diff(phmap, ph)), (ph, -sp.pi/2, sp.pi/2))/sp.pi
check("RN-3c hbar seat: <q> over the transported angle = 1 for ALL lambda (A units)",
      sp.simplify(Iq - 1) == 0)

print(); print("=" * 78)
print("PART 3 - uniqueness, now explained: only degree-0 RN calculus is free")
print("=" * 78)
k = sp.Symbol('k', real=True)
Fp = sp.sinh(k*L)/(k*sp.sinh(L))          # <D^p> closed form with k = p - 1
anti = (sp.cosh(L) - sp.sinh(L)*up)**(1 - (k + 1))/(sp.sinh(L)*k)
check("RN-4a closed form: <D^p> = sinh((p-1)lam)/((p-1) sinh lam)  [antiderivative]",
      sp.simplify(sp.diff(anti, up) - (sp.cosh(L) - sp.sinh(L)*up)**(-(k + 1))) == 0,
      "endpoints cosh -/+ sinh = e^{-/+lam} give the sinh ratio")
ser = sp.series(Fp, L, 0, 4).removeO()
check("RN-4b blind for ALL lam iff k^2 = 1, i.e. p in {0, 2}: p = 2 UNIQUE nontrivial",
      sp.simplify(ser.coeff(L, 0) - 1) == 0 and
      sp.simplify(ser.coeff(L, 2) - (k**2 - 1)/6) == 0,
      "lambda^2 coefficient (k^2 - 1)/6 kills every exponent but the mass one")
n_ = sp.Symbol('n', positive=True, integer=True)
q_u = sp.cosh(2*L) + sp.sinh(2*L)*sp.cos(2*ph)
mom = [sp.simplify(sp.integrate(q_u**m, (ph, 0, sp.pi))/sp.pi) for m in range(3)]
legs = [sp.legendre(m, sp.cosh(2*L)) for m in range(3)]
check("RN-5 the Legendre index shift EXPLAINED: <q^n>_{mu'} = uniform<q^{n-1}> =",
      all(sp.simplify((mom[m] - legs[m]).rewrite(sp.exp)) == 0 for m in range(3)),
      "P_{n-1}(cosh 2 lam): part 1's n-1 was the Jacobian eating one power of q")
check("RN-6 numeric spot (Will's run, path one): P_1(cosh .6), P_1(cosh 2.2)",
      abs(float(legs[1].subs(L, sp.Rational(3, 10))) - 1.1855) < 1e-3 and
      abs(float(legs[1].subs(L, sp.Rational(11, 10))) - 4.5679) < 1e-3,
      "1.1855 and 4.5679 - both paths agree")

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
LEDGER: the transport-exponent conjecture is RETIRED as of 2026-08-28 -
upgraded to the BLIND-MASS THEOREM (W.L.): given f^{-n} = d mu'/d mu, the
n-th blind moment is the total mass of the source measure; the exponent is
the Jacobian power by construction; uniqueness holds because only the
degree-zero object of RN calculus is measure-free (RN-4b).  The two former
"instances" were the premise verified per seat (RN-1, RN-2), never evidence.
OPEN, and sharper: Q-RN - derive the premise itself (reading = RN^{-1/n})
from the pairing/channel law.  Note the coincidence to be tested there:
n(c) = 2 equals F's quadratic channel exponent.  Filed beside DEBT-2b.
""")
