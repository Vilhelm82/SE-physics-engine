#!/usr/bin/env python3
# VIEW-1 incidence: the forced pivot circle through the three alignments (Will's incidence theorem, 2026-09-04).
# sympy exact; numpy for random frames.  Exit 0 iff every check passes.
# VIEW-1 [declared]: one pivot orbit must present all three constants.  Then C = S^2 ∩ Aff(a_1,a_2,a_3) [theorem | VIEW-1].
import sys, itertools, sympy as sp, numpy as np
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
a = [sp.Matrix(sp.symbols(f'a{i}1 a{i}2 a{i}3', real=True)) for i in (1, 2, 3)]      # GENERAL vectors, not unit
V = a[0].dot(a[1].cross(a[2]))
Avec = (a[0].cross(a[1]) + a[1].cross(a[2]) + a[2].cross(a[0]))/2
G = sp.Matrix(3, 3, lambda i, j: a[i].dot(a[j]))
one = sp.Matrix([1, 1, 1])
print("=== V-1  closed forms (any three vectors) ===")
check("V-1a 2 A.a_i = V identically for i = 1,2,3 (equidistance automatic)", all(sp.expand(2*Avec.dot(a[i]) - V) == 0 for i in range(3)))
check("V-1b 4|A|^2 = 1^T adj(G) 1 identically", sp.expand(4*Avec.dot(Avec) - (one.T*G.adjugate()*one)[0, 0]) == 0)
check("V-1c hence h^2 = V^2/(4|A|^2) = Delta/(1^T adj(G) 1) with Delta = det G = V^2", sp.expand(G.det() - V**2) == 0)
# the circumaxis in the DUAL basis: a_i x a_j = V eps_ijk a^k, so A_eps has dual coordinates (V/2)(e2 e3, e3 e1, e1 e2)
dual = [a[1].cross(a[2])/V, a[2].cross(a[0])/V, a[0].cross(a[1])/V]
ok_dual = True
for eps in [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]:
    ae = [eps[i]*a[i] for i in range(3)]
    Ae = (ae[0].cross(ae[1]) + ae[1].cross(ae[2]) + ae[2].cross(ae[0]))/2
    coords = (eps[1]*eps[2], eps[2]*eps[0], eps[0]*eps[1])
    Ae_pred = (V/2)*(coords[0]*dual[0] + coords[1]*dual[1] + coords[2]*dual[2])
    ok_dual &= all(sp.simplify(Ae[k] - Ae_pred[k]) == 0 for k in range(3))
check("V-1d POLES (theorem, whole cell): the circumaxis of the sign class eps is A_eps = (V/2) (e2e3 a^1 + e3e1 a^2 + e1e2 a^3) in the"
      " dual basis -- its dual coordinates ARE the Cayley node gamma_ij = eps_i eps_j.  'Axes ~ sum eps_i a_i' holds only at the trine", ok_dual)
print("=== V-2  the trine and the four body diagonals ===")
E = [sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])]
sub = {a[i][k]: E[i][k] for i in range(3) for k in range(3)}
hT = (V/(2*sp.sqrt(Avec.dot(Avec)))).subs(sub)
check("V-2a at the orthonormal frame h = 1/sqrt(3) and n = (1,1,1)/sqrt(3)", sp.simplify(hT - 1/sp.sqrt(3)) == 0
      and all(sp.simplify((Avec/sp.sqrt(Avec.dot(Avec))).subs(sub)[k] - 1/sp.sqrt(3)) == 0 for k in range(3)))
diags = set()
for eps in [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]:
    ae = [eps[i]*E[i] for i in range(3)]
    Ae = (ae[0].cross(ae[1]) + ae[1].cross(ae[2]) + ae[2].cross(ae[0]))/2
    diags.add(tuple(sp.sign(x) for x in Ae))
check("V-2b the four pole-resolved triples give the four body diagonals (up to sign)", len(diags) == 4 and all(abs(x) == 1 for d in diags for x in d), f"{sorted(diags)}")
print("=== V-3  degeneration and spin holonomy ===")
# unit-vector cell: parametrise a_1 = e1, a_2 = (c12, s12, 0), a_3 general unit; Delta -> 0 as the third leaves coplanarity
t2, t3, p3 = sp.symbols('t2 t3 p3', real=True)
u = [sp.Matrix([1, 0, 0]), sp.Matrix([sp.cos(t2), sp.sin(t2), 0]), sp.Matrix([sp.sin(t3)*sp.cos(p3), sp.sin(t3)*sp.sin(p3), sp.cos(t3)])]
subu = {a[i][k]: u[i][k] for i in range(3) for k in range(3)}
Vu = sp.simplify(V.subs(subu)); A2u = sp.simplify((Avec.dot(Avec)).subs(subu)); hu = Vu/(2*sp.sqrt(A2u))
check("V-3a h -> 0 exactly on the branch locus V = 0 (great circle iff coplanar), with h^2 = Delta/(1^T adj G 1) LINEAR in Delta:"
      " h ~ sqrt(Delta), the square class of V itself", sp.simplify(Vu - sp.sin(t2)*sp.cos(t3)) == 0
      and sp.simplify(hu**2*A2u*4 - Vu**2) == 0 and sp.simplify(hu.subs(t3, sp.pi/2)) == 0)
# spin lift: transport around a small circle of colatitude cos^-1(h) subtends cap 2 pi (1 - h); spin-1/2 lift exp(-i (cap/2) sigma.n)
h, al, be = sp.symbols('h alpha beta', real=True)
n = sp.Matrix([sp.sin(al)*sp.cos(be), sp.sin(al)*sp.sin(be), sp.cos(al)])
sig = [sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])]
sn = n[0]*sig[0] + n[1]*sig[1] + n[2]*sig[2]
I2 = sp.eye(2)
def rot(theta): return sp.cos(theta)*I2 + sp.I*sp.sin(theta)*sn          # exp(i theta sigma.n), (sigma.n)^2 = 1
check("V-3b (sigma.n)^2 = 1 for unit n", sp.simplify(sn*sn - I2) == sp.zeros(2))
R_cap = rot(-sp.pi*(1 - h)); R_will = -rot(sp.pi*h)
check("V-3c the cap lift exp(-i pi (1-h) sigma.n) equals Will's R = -exp(i pi h sigma.n) identically", sp.simplify(R_cap - R_will) == sp.zeros(2))
check("V-3d h -> 1 gives R -> +I; h -> 0 (the branch locus) gives R -> -I: the deck is reached exactly at V = 0",
      sp.simplify(R_will.subs(h, 1) - I2) == sp.zeros(2) and sp.simplify(R_will.subs(h, 0) + I2) == sp.zeros(2))
print("=== V-4  cyclic order and the deck (random frames, both lifts) ===")
rng = np.random.default_rng(4)
ok_out, ok_plus, ok_B = True, True, True
for _ in range(200):
    fr = [v/np.linalg.norm(v) for v in rng.normal(size=(3, 3))]
    for eps in [(1, 1, 1), (-1, -1, -1)]:                         # a state and its deck image tau: V -> -V
        f = [eps[i]*fr[i] for i in range(3)]
        Vn = np.dot(f[0], np.cross(f[1], f[2])); An = (np.cross(f[0], f[1]) + np.cross(f[1], f[2]) + np.cross(f[2], f[0]))/2
        nn = An/np.linalg.norm(An); hn = np.dot(nn, f[0])
        tri = np.cross(f[1] - f[0], f[2] - f[0])                     # oriented normal of the triangle a1 a2 a3
        ok_plus &= np.sign(np.dot(tri, nn)) > 0                      # seen from +n: ALWAYS counterclockwise (n is defined by A)
        outside = np.sign(hn)*nn                                     # the cap side: where an outside viewer sits
        ok_out &= np.sign(np.dot(tri, outside)) == np.sign(Vn)      # seen from OUTSIDE the sphere: the deck
        gam = [np.dot(f[0], f[1]), np.dot(f[0], f[2]), np.dot(f[1], f[2])]
        B = (1 + sum(gam) + 1j*Vn)/4
        ok_B &= np.sign(np.angle(B)) == np.sign(Vn) == np.sign(hn)
check("V-4a read from the +n side the cyclic order of (a1,a2,a3) is ALWAYS counterclockwise -- n is the triangle's own oriented normal,"
      " so that reading carries no deck information", ok_plus)
check("V-4b read from OUTSIDE the sphere (the cap side, sgn(h) n) the cyclic order equals sgn V, flips under tau, and equals sgn arg B"
      " with B = (1 + sum gamma + iV)/4: the geometric twin is the OUTSIDE reading, i.e. sgn h", ok_out and ok_B)
n_pass = sum(CH); print(f"\n{n_pass}/{len(CH)} checks passed")
sys.exit(0 if all(CH) else 1)
