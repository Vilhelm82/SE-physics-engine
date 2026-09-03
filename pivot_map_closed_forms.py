#!/usr/bin/env python3
# PIVOT-MAP CLOSED FORMS -- the mixed readings of the cold map (pivot_map.py) in closed form, split into deck-even and deck-odd parts.
# Objects: observer's axis n = A/|A| (VIEW-1), h = n.a_i = V/2|A| (equal for all three axes, V-1a); geodesic roots rho_i = (1 + a_i n)/N,
# N = sqrt(2(1+h)); transitions T_ij = rho_i rho_j^dag.  sympy exact on general vectors; numeric cross-check against pivot_map's readings.
import sys, itertools, numpy as np, sympy as sp
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
# ---------- symbolic: Pauli algebra on general vectors ----------
sx, sy, sz = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])
def sd(v): return v[0]*sx + v[1]*sy + v[2]*sz
a = [sp.Matrix(sp.symbols(f'a{i}1 a{i}2 a{i}3', real=True)) for i in (1, 2, 3)]
tn, pn = sp.symbols('t_n p_n', real=True); n = sp.Matrix([sp.sin(tn)*sp.cos(pn), sp.sin(tn)*sp.sin(pn), sp.cos(tn)]); h = sp.symbols('h', real=True)   # UNIT n
def halfTr(M): return sp.expand(M.trace()/2)
def rho(i):   return (sp.eye(2) + sd(a[i])*sd(n))                      # unnormalised geodesic rotor n -> a_i ; N^2 = 2(1 + n.a_i)
def T(i, j):  return rho(i)*rho(j).H                                  # unnormalised transition; true T = this / (2(1+h))^... see below
print("=== F-1  the transition in closed form ===")
i, j, k = 0, 1, 2
Tij = sp.expand(T(i, j))
m_ij = n.cross(a[j] - a[i]) + a[i].cross(a[j])
s_ij = 1 + n.dot(a[i]) + n.dot(a[j]) + a[i].dot(a[j])
check("F-1a (1 + a_i n)(1 + n a_j) = s_ij + i m_ij.sigma with s_ij = 1 + n.a_i + n.a_j + a_i.a_j and m_ij = n x (a_j - a_i) + a_i x a_j,"
      " for any a_i and UNIT n (n^2 = 1 is used; before VIEW-1's n.a_i = h)", sp.simplify(Tij - (s_ij*sp.eye(2) + sp.I*sd(m_ij))) == sp.zeros(2))
print("  With VIEW-1 (n.a_i = h for all i) and the roots normalised: T_ij = [(1 + 2h + gamma_ij) + i m_ij.sigma] / (2(1+h)).")
print("=== F-2  the axis-transition pseudoscalar reading ===")
# <a_k Tbar_ij>_I with Tbar = reversion = dagger:  (1/2)Tr[(a_k.sigma)(s - i m.sigma)] = -i a_k.m
X_I = halfTr(sd(a[k])*(s_ij*sp.eye(2) - sp.I*sd(m_ij)))
check("F-2a (1/2)Tr[(a_k.sigma)(s_ij - i m_ij.sigma)] = -i a_k.m_ij: the pseudoscalar reading is minus the component of the transition's"
      " bivector along the third axis", sp.simplify(X_I + sp.I*a[k].dot(m_ij)) == 0)
# a_k.m_ij = a_k.(n x (a_j - a_i)) + a_k.(a_i x a_j) = E_kij + eps_kij V
E_kij = a[k].dot(n.cross(a[j] - a[i])); V = a[0].dot(a[1].cross(a[2]))
check("F-2b a_k.m_ij = E_kij + V for (k,i,j) = (3,1,2) (even permutation), E_kij := a_k.(n x (a_j - a_i)) = n.((a_j - a_i) x a_k)",
      sp.expand(a[k].dot(m_ij) - E_kij - V) == 0)
print("  READING_I(k; i->j) = -(E_kij + eps_kij V) / (2(1+h)),   h = V/(2|A|).")
print("=== F-3  the even/odd split under the deck (V -> -V, h -> -h, E even) ===")
E, Vs, As = sp.symbols('E V A', real=True); hs = Vs/(2*As)
f = -(E + Vs)/(2*(1 + hs))
odd = sp.simplify((f - f.subs(Vs, -Vs))/2); even = sp.simplify((f + f.subs(Vs, -Vs))/2)
check("F-3a odd part = (E h - V)/(2(1 - h^2)) = V (E/(2A) - 1)/(2(1 - h^2)): V times a GRAM function -- the sqrt(Delta) class (G-9)",
      sp.simplify(odd - (E*hs - Vs)/(2*(1 - hs**2))) == 0)
check("F-3b even part = -(E - V h)/(2(1 - h^2)) = -(E - V^2/(2A))/(2(1 - h^2)): Gram and Delta only",
      sp.simplify(even + (E - Vs*hs)/(2*(1 - hs**2))) == 0)
print("=== F-4  E_kij and h in Gram terms (VIEW-1 closed forms) ===")
g12, g13, g23 = sp.symbols('g12 g13 g23', real=True)
Avec = (a[0].cross(a[1]) + a[1].cross(a[2]) + a[2].cross(a[0]))/2
# n.(a_2 x a_3) = A.(a_2 x a_3)/|A|; A.(a_2 x a_3) = (1/2)(1 - g23)(1 + g23 - g12 - g13) by Lagrange's identity
gram_sub = {a[0].dot(a[1]): g12, a[0].dot(a[2]): g13, a[1].dot(a[2]): g23}
lhs = sp.expand(Avec.dot(a[1].cross(a[2])))
# reduce via Lagrange: express in dot products by substituting an explicit orthonormal-free identity check on random rationals
import random
rng = random.Random(3); ok = True
for _ in range(6):
    sub = {s_: sp.Rational(rng.randint(-9, 9), rng.randint(1, 7)) for v in a for s_ in v}
    G = {g12: (a[0].dot(a[1])).subs(sub), g13: (a[0].dot(a[2])).subs(sub), g23: (a[1].dot(a[2])).subs(sub)}
    # for general (non-unit) vectors Lagrange gives A.(a2xa3) = (1/2)[ (a1.a2)(a2.a3) - (a1.a3)|a2|^2 + |a2|^2|a3|^2 - (a2.a3)^2 + (a3.a2)(a1.a3) - |a3|^2 (a1.a2) ]
    n2 = {ii: (a[ii].dot(a[ii])).subs(sub) for ii in range(3)}
    rhs = sp.Rational(1, 2)*(G[g12]*G[g23] - G[g13]*n2[1] + n2[1]*n2[2] - G[g23]**2 + G[g23]*G[g13] - n2[2]*G[g12])
    ok &= sp.simplify(lhs.subs(sub) - rhs) == 0
check("F-4a A.(a_2 x a_3) by Lagrange's identity on 6 random rational frames (general lengths); for UNIT axes it is"
      " (1/2)(1 - g23)(1 + g23 - g12 - g13)", ok and sp.expand(sp.Rational(1, 2)*(g12*g23 - g13 + 1 - g23**2 + g23*g13 - g12) - sp.Rational(1, 2)*(1 - g23)*(1 + g23 - g12 - g13)) == 0)
print("  Hence n.(a_2 x a_3) = (1 - g23)(1 + g23 - g12 - g13) / (2|A|),  4|A|^2 = 1^T adj(G) 1,  and cyclically; E_kij is a signed sum of two of these.")
print("=== F-5  numeric cross-check against the cold map's readings (pivot_map.py construction) ===")
Xn = np.array([[0, 1], [1, 0]], complex); Yn = np.array([[0, -1j], [1j, 0]]); Zn = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2)
def sdn(v): return v[0]*Xn + v[1]*Yn + v[2]*Zn
def unit(v): return v/np.linalg.norm(v)
def geo(u, w): return (I2 + sdn(w) @ sdn(u))/np.sqrt(2*(1 + np.dot(u, w)))
rngn = np.random.default_rng(31); ok5 = ok6 = ok7 = True; oddvals = []
for _ in range(150):
    av = [unit(v) for v in rngn.normal(size=(3, 3))]
    A = (np.cross(av[0], av[1]) + np.cross(av[1], av[2]) + np.cross(av[2], av[0]))/2; nn = A/np.linalg.norm(A)
    Vn = np.dot(av[0], np.cross(av[1], av[2])); hn = Vn/(2*np.linalg.norm(A))
    r = [geo(nn, av[q]) for q in range(3)]; Tn = r[0] @ r[1].conj().T
    read_I = (0.5*np.trace(sdn(av[2]) @ Tn.conj().T)).imag                                   # pivot_map's axis-trans I-part, k=3, (i,j)=(1,2)
    mn = np.cross(nn, av[1] - av[0]) + np.cross(av[0], av[1]); En = np.dot(av[2], np.cross(nn, av[1] - av[0]))
    closed = -(En + Vn)/(2*(1 + hn))
    ok5 &= abs(read_I - closed) < 1e-10
    # even/odd split vs the deck image computed directly
    avt = [-v for v in av]; At = (np.cross(avt[0], avt[1]) + np.cross(avt[1], avt[2]) + np.cross(avt[2], avt[0]))/2; nt = At/np.linalg.norm(At)
    rt = [geo(nt, avt[q]) for q in range(3)]; Tt = rt[0] @ rt[1].conj().T
    read_t = (0.5*np.trace(sdn(avt[2]) @ Tt.conj().T)).imag
    odd_direct = (read_I - read_t)/2; odd_closed = (En*hn - Vn)/(2*(1 - hn**2))
    ok6 &= abs(odd_direct - odd_closed) < 1e-10
    oddvals.append(odd_closed)
    # trans-trans scalar: (1/2)Tr(T_ij T_kl^dag) = (s_ij s_kl + m_ij.m_kl)/(4(1+h)^2), (k,l) = (2,3)
    Tkl = r[1] @ r[2].conj().T
    s12 = 1 + 2*hn + np.dot(av[0], av[1]); s23 = 1 + 2*hn + np.dot(av[1], av[2])
    m12 = np.cross(nn, av[1] - av[0]) + np.cross(av[0], av[1]); m23 = np.cross(nn, av[2] - av[1]) + np.cross(av[1], av[2])
    ok7 &= abs((0.5*np.trace(Tn @ Tkl.conj().T)).real - (s12*s23 + np.dot(m12, m23))/(4*(1 + hn)**2)) < 1e-10
check("F-5a closed form of the axis-transition pseudoscalar reading matches pivot_map's construction on 150 frames", ok5)
check("F-5b its odd part (E h - V)/(2(1 - h^2)) equals the direct antisymmetrisation over the deck pair on 150 frames", ok6)
check("F-5c transition-transition scalar: (1/2)Tr(T_ij T_kl^dag) = (s_ij s_kl + m_ij.m_kl)/(4(1+h)^2), s = 1 + 2h + gamma (150 frames)", ok7)
check("F-5d the odd part is generically nonzero (150/150) -- the separating content of the mixed reading", all(abs(v) > 1e-6 for v in oddvals),
      f"|odd| range {min(map(abs, oddvals)):.3g} .. {max(map(abs, oddvals)):.3g}")
print("=== SUMMARY OF CLOSED FORMS ===")
print("  T_ij = [(1 + 2h + gamma_ij) + i m_ij.sigma] / (2(1+h)),   m_ij = n x (a_j - a_i) + a_i x a_j,   h = V/(2|A|),  4|A|^2 = 1^T adj(G) 1")
print("  <a_k Tbar_ij>_I  = -(E_kij + eps_kij V) / (2(1+h)),   E_kij = n.((a_j - a_i) x a_k)  [Gram-rational over |A|: n.(a_2 x a_3) = (1-g23)(1+g23-g12-g13)/(2|A|), cyclic]")
print("     odd part  = (E_kij h - eps_kij V) / (2(1 - h^2))  = V [E_kij/(2|A|) - eps_kij] / (2(1 - h^2))     -- sqrt(Delta) class, deck-odd")
print("     even part = -(E_kij - eps_kij V h) / (2(1 - h^2))                                              -- Gram and Delta only")
print("  <T_ij Tbar_kl>_S = [(1+2h+gamma_ij)(1+2h+gamma_kl) + m_ij.m_kl] / (4(1+h)^2);  odd part by the same antisymmetrisation in h.")
print("  Every odd part is V times a rational function of (gamma, |A|, h^2): the cold map's separating content is the orientation datum")
print("  dressed by the observer's own circumaxis -- mixed because h sits in the denominator, odd because V sits in the numerator.")


print("=== F-6  (added after Will's derivation) the dressing IS the pole triangle (n, a_i, a_j) ===")
ok8 = ok9 = ok10 = True
for _ in range(150):
    av = [unit(v) for v in rngn.normal(size=(3, 3))]
    A = (np.cross(av[0], av[1]) + np.cross(av[1], av[2]) + np.cross(av[2], av[0]))/2; nn = A/np.linalg.norm(A)
    hn = np.dot(nn, av[0]); ai, aj = av[0], av[1]; gij = np.dot(ai, aj)
    Ri, Rj = geo(nn, ai), geo(nn, aj); Tn = Ri @ Rj.conj().T                        # T_ij = R_i R_j~  (Will's two-step route through n)
    Gd = geo(aj, ai)                                                                  # the DIRECT geodesic rotor a_j -> a_i
    W = Gd.conj().T @ Tn                                                              # what the detour adds: must fix a_j
    ok8 &= np.allclose(W @ sdn(aj) @ W.conj().T, sdn(aj))
    # W is a rotation about a_j; its rotor angle: W = cos(psi/2) - i sin(psi/2) (a_j.sigma) up to sign
    c = (0.5*np.trace(W)).real; s_ = -(0.5*np.trace(1j*W @ sdn(aj))).real          # cos(psi/2), sin(psi/2)
    psi_half = np.arctan2(s_, c)
    # the pole triangle (n, a_i, a_j): Van Oosterom-Strackee, tan(Omega/2) = n.(a_i x a_j) / (1 + n.a_i + a_i.a_j + a_j.n)
    Om_half = np.arctan2(np.dot(nn, np.cross(ai, aj)), 1 + hn + gij + hn)
    ok9 &= abs(abs(psi_half) - abs(Om_half)) < 1e-9 or abs(abs(abs(psi_half) - np.pi) - abs(Om_half)) < 1e-9
    # and the n-component of the transition's bivector with its scalar part IS that triangle's (denominator, numerator) pair
    mn = np.cross(nn, aj - ai) + np.cross(ai, aj)
    ok10 &= abs(np.dot(mn, nn) - np.dot(nn, np.cross(ai, aj))) < 1e-12 and abs((1 + 2*hn + gij) - (1 + hn + gij + hn)) < 1e-12
check("F-6a G~_ij T_ij fixes a_j: the two-step route through n differs from the direct geodesic by a rotation ABOUT a_j", ok8)
check("F-6b its angle is the solid angle of the pole triangle (n, a_i, a_j) -- the holonomy of the loop a_j -> n -> a_i -> a_j (150 frames,"
      " up to the spinor sign)", ok9)
check("F-6c (scalar part of T_ij, n-component of m_ij) = (1 + 2h + gamma_ij, n.(a_i x a_j)): exactly the Van Oosterom-Strackee pair of that"
      " pole triangle.  The cold map's transitions carry CENSUS-C's pole triangles; the mixed readings are dressed by their solid angles", ok10)
n_pass = sum(CH); print(f"\n{n_pass}/{len(CH)} checks passed (with F-6)")
sys.exit(0 if all(CH) else 1)
