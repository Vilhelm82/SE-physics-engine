#!/usr/bin/env python3
# =============================================================================
# THE SEAT-CHANGE COCYCLE AND THE GRAPH DESCENT THEOREM  (2026-08-28, encore 2)
# Will's synthesis: presentation factor = multiplicative cocycle on a seat
# groupoid; complexified legs tau_ij; cycle products = holonomies; descent.
# DAG checked first: 'cocycle' 0 hits; 'holonomy' = the discrete-transport
# analogy only. This layer is NOT in the corpus. Deriving.
# =============================================================================
import sympy as sp, mpmath as mp, random

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" - {n}" if n else ""), flush=True)

u = sp.Symbol('u', real=True)
L1, L2 = sp.symbols('lambda1 lambda2', positive=True)
s  = lambda L, x: sp.cosh(L) + sp.sinh(L)*x
ab = lambda L, x: (x + sp.tanh(L))/(1 + sp.tanh(L)*x)     # aberration action

print("=" * 78); print("PART 1 - the groupoid cocycle, verified clause by clause"); print("=" * 78)
lhs = s(L2, ab(L1, u)) * s(L1, u)
check("C-1  COCYCLE: sigma_{g2}(g1 x) sigma_{g1}(x) = sigma_{g2 g1}(x), exact",
      sp.simplify((lhs - s(L1 + L2, u)).rewrite(sp.exp)) == 0,
      "the presentation factor is a multiplicative cocycle on the boost groupoid")
check("C-2  INVERSE: sigma_{g^-1}(g x) = sigma_g(x)^{-1}  (seat reversal law)",
      sp.simplify((s(-L1, ab(L1, u))*s(L1, u) - 1).rewrite(sp.exp)) == 0,
      "reversal sends sigma -> sigma^{-1}, as claimed")
k = sp.Symbol('k', real=True)
Jk = sp.sinh((k + 1)*L1)/((k + 1)*sp.sinh(L1))            # int s^k du/2, closed
check("C-3  MIRROR = the k <-> -k-2 symmetry of the spectrum: J(k) = J(-k-2)",
      sp.simplify(Jk - Jk.subs(k, -k - 2)) == 0,
      "with k = p - Qm, (p,m)->(-p,1-m) is EXACTLY this reflection")
print("       HONEST DOWNGRADE: this is the alpha <-> 1-alpha skew symmetry of")
print("       the Renyi/Chernoff integral between the two seat measures; the")
print("       fixed point m = 1/2 is the Bhattacharyya/Hellinger midpoint.")
print("       The half-seat reflection loses 'unlocated' status: it is the")
print("       Renyi skew symmetry, instantiated on the Doppler pair. [lineage]")
check("C-4  det-1 multipliers: the pole pair e^{+-lambda} multiplies to 1",
      sp.simplify(s(L1, 1)*s(L1, -1) - 1) == 0,
      "loxodromic multiplier product = 1: E-8c named in Moebius language")

print(); print("=" * 78)
print("PART 2 - THE GRAPH DESCENT THEOREM (constructive, n = 4 and n = 5)")
print("=" * 78)
# Binet-Cauchy: [u1u2u3][w1w2w3] = det(u_i . w_j) - the relation engine
U = [sp.Matrix(sp.symbols(f'u{i}_1:4')) for i in range(1, 4)]
W = [sp.Matrix(sp.symbols(f'w{i}_1:4')) for i in range(1, 4)]
lhsBC = sp.Matrix.hstack(*U).det() * sp.Matrix.hstack(*W).det()
rhsBC = sp.Matrix([[(U[i].T*W[j])[0,0] for j in range(3)] for i in range(3)]).det()
check("D-1  BINET-CAUCHY: [u1 u2 u3][w1 w2 w3] = det(u_i . w_j), exact polynomial",
      sp.expand(lhsBC - rhsBC) == 0,
      "=> ALL triple volumes on shared rays are rationally related through ONE")
print("       radical: V(1jk) = M_jk / V(123) with M_jk Gram-rational. The rank-3")
print("       ambient is the mechanism; this is what a general CP^N graph lacks.")

mp.mp.dps = 30
def mvec():
    v = [mp.mpf(random.gauss(0, 1)) for _ in range(3)]
    n = mp.sqrt(sum(x*x for x in v)); return [x/n for x in v]
def mspin(a):
    th = mp.acos(a[2]); ph = mp.atan2(a[1], a[0])
    return [mp.cos(th/2), mp.e**(1j*ph)*mp.sin(th/2)]
def ovl(x, y): return mp.conj(x[0])*y[0] + mp.conj(x[1])*y[1]
def gd(a, b): return sum(x*y for x, y in zip(a, b))
def trip(a, b, c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))
def Sfun(g12, g23, g13): return 1 + g12 + g23 + g13
random.seed(41); worst4 = worst5 = mp.mpf(0)
for _ in range(6):
    A = [mvec() for _ in range(5)]
    S1 = mspin(A[0]); sp2, sp3, sp4, sp5 = map(mspin, A[1:])
    g = lambda i, j: gd(A[i], A[j])
    V1 = trip(A[0], A[1], A[2])                       # the ONE radical V0
    # 4-cycle direct
    B4 = ovl(S1, sp2)*ovl(sp2, sp3)*ovl(sp3, sp4)*ovl(sp4, S1)
    # descent formula: only Gram data + V1 appear
    Sa = Sfun(g(0,1), g(1,2), g(0,2)); Sb = Sfun(g(0,2), g(2,3), g(0,3))
    M  = (mp.matrix([[1, g(0,2), g(0,3)],
                     [g(0,1), g(1,2), g(1,3)],
                     [g(0,2), 1, g(2,3)]]))
    Mdet = mp.det(M)                                  # = V(123) V(134), rational
    D1 = 1 - g(0,1)**2 - g(0,2)**2 - g(1,2)**2 + 2*g(0,1)*g(0,2)*g(1,2)
    B4f = ((Sa*Sb - Mdet) + 1j*V1*(Sb + Sa*Mdet/D1))/(8*(1 + g(0,2)))
    worst4 = max(worst4, abs(B4 - B4f))
check("D-2  4-CYCLE DESCENT: B(1234) = [Gram-rational] + [Gram-rational]*V0, to",
      worst4 < mp.mpf(10)**-25, f"1e-25 (worst {mp.nstr(worst4,3)}) - degree <= 2 over F, constructive")

random.seed(43)
for _ in range(6):
    A = [mvec() for _ in range(5)]
    S = [mspin(a) for a in A]
    g = lambda i, j: gd(A[i], A[j])
    V1 = trip(A[0], A[1], A[2])
    B5 = (ovl(S[0],S[1])*ovl(S[1],S[2])*ovl(S[2],S[3])
          *ovl(S[3],S[4])*ovl(S[4],S[0]))
    # anchored decomposition: B(12345) = B(123) B(134) B(145)
    #                                     / (|<a1|a3>|^2 |<a1|a4>|^2)
    def Btri(i, j, k):
        Sv = Sfun(g(i,j), g(j,k), g(i,k)); Vv = trip(A[i], A[j], A[k])
        return (Sv + 1j*Vv)/4
    # every V(1jk) reduced to V1 by Binet-Cauchy:
    def Vred(j, k):
        Mm = mp.matrix([[1, g(0,j), g(0,k)],
                        [g(0,1), g(1,j), g(1,k)],
                        [g(0,2), g(2,j), g(2,k)]])
        return mp.det(Mm)/V1
    def BtriRed(j, k):
        return (Sfun(g(0,j), g(j,k), g(0,k)) + 1j*Vred(j, k))/4
    B5f = (BtriRed(1,2)*BtriRed(2,3)*BtriRed(3,4)
           / (((1+g(0,2))/2) * ((1+g(0,3))/2)))
    worst5 = max(worst5, abs(B5 - B5f))
check("D-3  5-CYCLE DESCENT: built from Gram data + the SAME single V0, to 1e-25",
      worst5 < mp.mpf(10)**-25, f"worst {mp.nstr(worst5,3)} - beta_1 grew; the field did not")

print(); print("=" * 78)
print("PART 3 - coboundaries, the Z_2, and where beta_1 actually lives")
print("=" * 78)
random.seed(47); A = [mvec() for _ in range(5)]; S = [mspin(a) for a in A]
ph = [mp.e**(1j*mp.mpf(random.uniform(0, 6.28))) for _ in range(5)]
Sg = [[p*c for c in s_] for p, s_ in zip(ph, S)]
B0 = ovl(S[0],S[1])*ovl(S[1],S[2])*ovl(S[2],S[3])*ovl(S[3],S[4])*ovl(S[4],S[0])
Bg = ovl(Sg[0],Sg[1])*ovl(Sg[1],Sg[2])*ovl(Sg[2],Sg[3])*ovl(Sg[3],Sg[4])*ovl(Sg[4],Sg[0])
check("E-1  vertex rephasings are COBOUNDARIES: 5-cycle holonomy gauge-invariant",
      abs(B0 - Bg) < mp.mpf(10)**-25, "open legs are cocycle data; closed cycles are cohomology")
Ar = [[-x for x in a] for a in A]; Sr = [mspin(a) for a in Ar]
Br = ovl(Sr[0],Sr[1])*ovl(Sr[1],Sr[2])*ovl(Sr[2],Sr[3])*ovl(Sr[3],Sr[4])*ovl(Sr[4],Sr[0])
check("E-2  ONE Z_2 rules every cycle: global reversal conjugates the 5-cycle too",
      abs(Br - mp.conj(B0)) < mp.mpf(10)**-25,
      "Gal(F(V0)/F) acts by simultaneous conjugation of ALL holonomies")

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
GRAPH DESCENT THEOREM (candidate, banked tonight):
  For a connected graph on n rays in CP^1 (unit vectors in R^3), generic Gram:
  (i)   individual legs need the conjugate Kummer radicals sqrt(1+g_ij) - the
        leg field grows with the graph;
  (ii)  vertex rephasings are coboundaries; the U(1) holonomy TORUS has
        dimension beta_1 = |E| - |V| + 1  [Will's beta_1, placed: it counts
        holonomy COORDINATES, i.e. cohomology];
  (iii) yet EVERY closed-cycle Bargmann product lies in F(V0), a SINGLE
        quadratic extension, with explicit [rational] + [rational]*V0 forms
        via anchored triangle decomposition + Binet-Cauchy (D-1..D-3):
        THE FIELD DOES NOT GROW WITH beta_1. Degree <= 2, always.
  (iv)  Gal(F(V0)/F) = Z_2 = global reversal, conjugating all cycles at once.
  MECHANISM AND BOUNDARY: the collapse is forced by the rank-3 ambient
  (Binet-Cauchy relations among triple volumes). For rays in CP^N, N >= 2,
  those relations weaken - the general-N cycle field is OPEN and likely
  grows. That is the 'considerably subtler' part, located precisely.
  STATUS: constructive for cycles through n = 5 + the general mechanism;
  full arbitrary-graph write-up, genericity conditions (V0 != 0, overlaps
  != 0), and a novelty search are next-session work.
""")
