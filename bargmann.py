#!/usr/bin/env python3
# =============================================================================
# THE BARGMANN PHASE OF THE SEAT CYCLE  (2026-08-28, the encore)
# Lift each axis to its minimal spinor (half-angle, the 2:1 cover); the cycle
# B = <a1|a2><a2|a3><a3|a1> has |B| Gram-determined and arg(B) = -Omega/2,
# the SIGNED solid angle of the frame triangle [Bargmann; Pancharatnam;
# Mukunda-Simon].  Van Oosterom-Strackee: tan(Omega/2) = V / (1 + Sum gamma).
# The phase generator is V: THE MEDIUM.  Lineage loud; model-native content:
# the numerator = J, the sheet-detector, the branch null, the trine pi-event.
# =============================================================================
import sympy as sp
import mpmath as mp
import random

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" - {n}" if n else ""), flush=True)

print("=" * 78); print("PART 1 - the lift (symbolic)"); print("=" * 78)
t1, p1, t2, p2 = sp.symbols('theta1 phi1 theta2 phi2', real=True)
def vec(th, ph): return sp.Matrix([sp.sin(th)*sp.cos(ph), sp.sin(th)*sp.sin(ph), sp.cos(th)])
def spinor(th, ph): return sp.Matrix([sp.cos(th/2), sp.exp(sp.I*ph)*sp.sin(th/2)])
ov = (spinor(t1, p1).conjugate().T * spinor(t2, p2))[0, 0]
lhs = sp.expand_trig(sp.expand(sp.re(ov)**2 + sp.im(ov)**2))
rhs = sp.expand_trig(sp.expand((1 + (vec(t1, p1).T*vec(t2, p2))[0, 0])/2))
check("BP-1 minimal lift: |<a|b>|^2 = (1 + a.b)/2, exact  (half-angle = the cover)",
      sp.simplify(lhs - rhs) == 0)

print(); print("=" * 78)
print("PART 2 - the phase = -Omega/2, two paths (30-digit)"); print("=" * 78)
mp.mp.dps = 30
def mvec():
    v = [mp.mpf(random.gauss(0, 1)) for _ in range(3)]
    n = mp.sqrt(sum(x*x for x in v)); return [x/n for x in v]
def mspin(a):
    th = mp.acos(a[2]); ph = mp.atan2(a[1], a[0])
    return [mp.cos(th/2), mp.e**(1j*ph)*mp.sin(th/2)]
def bargmann(A):
    s = [mspin(a) for a in A]
    ovl = lambda x, y: mp.conj(x[0])*y[0] + mp.conj(x[1])*y[1]
    return ovl(s[0], s[1]) * ovl(s[1], s[2]) * ovl(s[2], s[0])
def gdata(A):
    g = [sum(A[0][k]*A[1][k] for k in range(3)), sum(A[0][k]*A[2][k] for k in range(3)),
         sum(A[1][k]*A[2][k] for k in range(3))]
    V = (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]) - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
         + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))
    return g, V
def Omega(A):
    g, V = gdata(A); return 2*mp.atan2(V, 1 + g[0] + g[1] + g[2])

oct_ = [[0.0,0,1.0],[1.0,0,0],[0.0,1.0,0]]
anchor = abs(mp.arg(bargmann(oct_)) - mp.pi/4) < mp.mpf(10)**-25 and \
         abs(Omega(oct_) - mp.pi/2) < mp.mpf(10)**-25
check("BP-2a ANCHOR (octant): Omega = pi/2, arg(B) = +pi/4: this chart gives +Omega/2",
      anchor, "sign is chart convention; structure unchanged  [Bargmann, lineage]")
random.seed(11)
ok = True; worst = mp.mpf(0)
for _ in range(8):
    A = [mvec(), mvec(), mvec()]
    d = abs(mp.e**(1j*(mp.arg(bargmann(A)) - Omega(A)/2)) - 1)
    worst = max(worst, d); ok = ok and d < mp.mpf(10)**-25
check("BP-2b arg(B) = +Omega/2, tan(Omega/2) = V/(1+Sum gamma): spinor vs vector path",
      ok, f"8 random frames, worst {mp.nstr(worst,3)}  [VOS + Bargmann, lineage]")
A = [mvec(), mvec(), mvec()]; B0 = bargmann(A)
Bc = bargmann([A[1], A[2], A[0]]); Bs = bargmann([A[1], A[0], A[2]])
check("BP-3 circulation preserves B; a SWAP conjugates it: untwisted/twisted as",
      abs(Bc - B0) < mp.mpf(10)**-25 and abs(Bs - mp.conj(B0)) < mp.mpf(10)**-25,
      "complex conjugation of the seat-cycle holonomy - the click split, realised")
Ar = [[-x for x in a] for a in A]
gA, _ = gdata(A); gR, _ = gdata(Ar)
check("BP-4 global reversal: G fixed, Phi -> -Phi: a G-BLIND sign - the phase is",
      max(abs(gA[k]-gR[k]) for k in range(3)) < mp.mpf(10)**-28
      and abs(bargmann(Ar) - mp.conj(B0)) < mp.mpf(10)**-25,
      "the interferometric detector for the orientation sheet [B'-monodromy's meter]")
cop = [[1.0,0,0], [mp.mpf('-0.5'), mp.sqrt(3)/2, 0], [mp.mpf('0.2'), mp.mpf('0.9'), 0]]
n3 = mp.sqrt(cop[2][0]**2 + cop[2][1]**2); cop[2] = [x/n3 for x in cop[2]]
check("BP-5 detector NULLS on the branch locus: coplanar frame => Phi = 0",
      abs(mp.arg(bargmann(cop))) < mp.mpf(10)**-25,
      "the sheets meet where the meter reads zero")
print()
print("  BP-6 the pi-event, located: tilt the symmetric trine out of plane:")
prev = 0
okp = True
for eps in ('0.3', '0.1', '0.03', '0.01'):
    e = mp.mpf(eps)
    a3 = [mp.mpf('-0.5'), -mp.sqrt(3)/2, e]
    n3 = mp.sqrt(sum(x*x for x in a3)); a3 = [x/n3 for x in a3]
    A = [[1.0,0,0], [mp.mpf('-0.5'), mp.sqrt(3)/2, 0], a3]
    Ph = mp.arg(bargmann(A))
    print(f"      eps = {eps:>5}:  Phi = {mp.nstr(Ph, 10)}")
    okp = okp and abs(Ph) > abs(prev); prev = Ph
check("BP-6 |Phi| -> pi approaching the trine: the pi-phase event has an address",
      okp and abs(prev) > 3.0, f"|Phi(0.01)| = {mp.nstr(abs(prev),8)} -> pi")
tau = mp.mpf('0.3')
def frameT(tv):
    a3 = [mp.mpf('-0.5'), mp.mpf('-0.43'), tv]
    n3 = mp.sqrt(sum(x*x for x in a3)); return [[1.0,0,0],[mp.mpf('-0.5'), mp.sqrt(3)/2, 0],[x/n3 for x in a3]]
Pp, Pm = mp.arg(bargmann(frameT(tau))), mp.arg(bargmann(frameT(-tau)))
P0 = mp.arg(bargmann(frameT(mp.mpf(0))))
check("BP-7 crossing: Phi(-tau) = -Phi(tau) mod 2pi; Phi(0) in {0, pi} - the branch",
      abs(mp.e**(1j*(Pp + Pm)) - 1) < mp.mpf(10)**-25 and
      (abs(P0) < mp.mpf(10)**-25 or abs(abs(P0) - mp.pi) < mp.mpf(10)**-25),
      f"locus is Z2-GRADED by the meter: 0-class / pi-class, split at the trine; Phi(0) = {mp.nstr(P0,8)}")
print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
BANKED: PREDICTION-1.  J = (1 + g12 + g13 + g23) * tan(Omega/2),  Omega = +2 Phi (this chart)
(minimal lift s = 1/2; general lift multiplies Phi by 2s - charge is a labelled-
tier datum, candidate tie to the hbar quantum: CONJECTURE).  Mathematics exact
(BP-2); physics status: the medium is interferometrically accessible via the
seat-cycle holonomy - the model's first clean prediction-shaped statement.
The phase generator is V: the sand IS the phase.  It nulls on the branch,
saturates to pi at the trine, conjugates under reversals, and reads the sheet.
""")
