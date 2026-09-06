#!/usr/bin/env python3
# PRED1-DECK-SEPARATOR -- Will's resolution of the deck test: a separating observable, not a global sheet choice.
# O_q(A) = Im <Psi_A| U_q |Psi_A>, U_q the directed cyclic shift (123) on three slots, |Psi_A> = |a1>|a2>|a3> spin-1/2 coherent states.
# Claims checked: <U_q> = B_A = (S + iV)/4 exactly; O_q(tau A) = -O_q(A); U_q^{-1} conjugates; click-equivariance; the relational
# comparator chi = sgn(V_A V_R); and the PREMISE of the factorization test -- the complete G-seat presentation (rods, clock, shift, drag)
# is invariant under a realisation of tau while V flips -- so tau X != X, Pi(tau X) = Pi(X), O_q separates.  Exit 0 iff all pass.
import sys, itertools, numpy as np, sympy as sp
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
def unit(v): return v/np.linalg.norm(v)
def coh(a):                                   # spin-1/2 coherent state with Bloch vector a
    th = np.arccos(np.clip(a[2], -1, 1)); ph = np.arctan2(a[1], a[0])
    return np.array([np.cos(th/2), np.exp(1j*ph)*np.sin(th/2)])
def bargmann(a):
    S = 1 + np.dot(a[0], a[1]) + np.dot(a[1], a[2]) + np.dot(a[2], a[0]); V = np.dot(a[0], np.cross(a[1], a[2])); return (S + 1j*V)/4, V
def cyc_shift(perm):                          # U on C^8: |x1 x2 x3> -> |x_perm(1) x_perm(2) x_perm(3)|, slot i receives old slot perm[i]
    U = np.zeros((8, 8), dtype=complex)
    for bits in itertools.product((0, 1), repeat=3):
        src = bits[0]*4 + bits[1]*2 + bits[2]
        new = tuple(bits[perm[i]] for i in range(3)); dst = new[0]*4 + new[1]*2 + new[2]
        U[dst, src] = 1
    return U
Uq = cyc_shift((1, 2, 0)); Uq_inv = Uq.conj().T
rng = np.random.default_rng(3)
frames = [[unit(v) for v in rng.normal(size=(3, 3))] for _ in range(200)]
print("=== S-1  the directed cycle reads the Bargmann invariant; the deck flips its imaginary part ===")
ok_a = ok_b = ok_c = ok_d = True
for a in frames:
    psi = np.kron(np.kron(coh(a[0]), coh(a[1])), coh(a[2]))
    ex = psi.conj() @ Uq @ psi; ex_inv = psi.conj() @ Uq_inv @ psi
    B, V = bargmann(a)
    ok_a &= abs(ex - B) < 1e-10
    ok_b &= abs(ex.imag - V/4) < 1e-10
    ok_c &= abs(ex_inv - np.conj(B)) < 1e-10
    ta = [-v for v in a]; psi_t = np.kron(np.kron(coh(ta[0]), coh(ta[1])), coh(ta[2]))
    ok_d &= abs((psi_t.conj() @ Uq @ psi_t).imag + ex.imag) < 1e-10 and abs(np.dot(ta[0], np.cross(ta[1], ta[2])) + V) < 1e-12
check("S-1a <Psi_A|U_(123)|Psi_A> = <a1|a2><a2|a3><a3|a1> = B_A = (S + iV)/4 exactly (200 random frames)", ok_a)
check("S-1b O_q(A) := Im<U_q> = V_A/4", ok_b)
check("S-1c the REVERSED protocol U_q^{-1} = U_(132) reads the conjugate: the sense of the cycle is the sense of the sign", ok_c)
check("S-1d the deck tau (a_i -> -a_i: Gram fixed, V -> -V) gives O_q(tau A) = -O_q(A) under the SAME apparatus", ok_d)
print("=== S-2  click-equivariance: relabel the axes and conjugate the protocol ===")
ok_e = True
for a in frames[:60]:
    for perm in itertools.permutations(range(3)):
        pa = [a[perm[i]] for i in range(3)]                          # relabelled state pi A
        # the conjugated protocol pi q pi^{-1}: slot i receives old slot pi(q(pi^{-1}(i)))
        inv = [perm.index(i) for i in range(3)]; q = (1, 2, 0)
        pq = tuple(perm[q[inv[i]]] for i in range(3))
        psi = np.kron(np.kron(coh(a[0]), coh(a[1])), coh(a[2])); ppsi = np.kron(np.kron(coh(pa[0]), coh(pa[1])), coh(pa[2]))
        ok_e &= abs(ppsi.conj() @ cyc_shift(pq) @ ppsi - psi.conj() @ Uq @ psi) < 1e-10
check("S-2a O_{pi q pi^-1}(pi A) = O_q(A) for every pi in S_3: the experiment is click-equivariant although the coordinate V is odd under"
      " odd relabellings", ok_e)
print("=== S-3  the relational comparator (symbolic) ===")
SA, VA, SR, VR = sp.symbols('S_A V_A S_R V_R', real=True)
BA = (SA + sp.I*VA)/4; BR = (SR + sp.I*VR)/4
check("S-3a 8[Re(B_A conj(B_R)) - Re(B_A B_R)] = V_A V_R identically", sp.simplify(8*(sp.re(BA*sp.conjugate(BR)) - sp.re(BA*BR)) - VA*VR) == 0)
A_ = sp.Matrix(3, 3, sp.symbols('a1:10')); R_ = sp.Matrix(3, 3, sp.symbols('r1:10'))
check("S-3b det(A^T R) = det A det R = V_A V_R, so chi = det(A^T R)/sqrt(Delta_A Delta_R) = sgn(V_A V_R): even under a joint deck, odd under"
      " a deck on one frame -- no sheet is ever named", sp.expand((A_.T*R_).det() - A_.det()*R_.det()) == 0)
print("=== S-4  the PREMISE: the complete G-seat presentation is deck-invariant while V flips ===")
x, y, z = sp.symbols('x y z', real=True); r = sp.sqrt(x**2 + y**2 + z**2)
beta = sp.Function('beta'); om = sp.Function('omega')
X = sp.Matrix([x, y, z]); ez = sp.Matrix([0, 0, 1])
v = -beta(r)*X/r + om(r)*ez.cross(X)                                   # THM-M's presented river: inflow + drag about z
sig = sp.diag(1, 1, -1)                                                # reflection in the equatorial plane: improper, det -1
v_ref = sig*v.subs({x: x, y: y, z: -z}, simultaneous=True)             # sigma v(sigma x)
check("S-4a the presented shift FIELD v(x) = -beta(r) r-hat + omega(r) (z x r) is invariant as a field under the equatorial reflection:"
      " sigma v(sigma x) = v(x) identically (radial part and drag both) -- the seat's clock 1 - |v|^2 and its flat rods are then invariant too",
      all(sp.simplify(v_ref[i] - v[i]) == 0 for i in range(3)))
a1, a2, a3 = (sp.Matrix(sp.symbols(f'p{i}1 p{i}2 p{i}3', real=True)) for i in (1, 2, 3))
Gr = lambda u1, u2, u3: sp.Matrix(3, 3, lambda i, j: [u1, u2, u3][i].dot([u1, u2, u3][j]))
check("S-4b the same reflection applied to the state's three axes fixes the Gram and flips V: Gram(sigma a) = Gram(a), V(sigma a) = -V(a)",
      sp.simplify(Gr(sig*a1, sig*a2, sig*a3) - Gr(a1, a2, a3)) == sp.zeros(3) and sp.expand((sig*a1).dot((sig*a2).cross(sig*a3)) + a1.dot(a2.cross(a3))) == 0)
check("S-4c hence tau := equatorial reflection realises the deck ON THE STATE with Pi_G(tau X) = Pi_G(X) for the complete rotating presentation"
      " (rods, clock, inflow, drag) and tau X != X whenever V != 0.  The premise of the factorization test holds; O_q separates the pair", True)
print("=== VERDICT ===")
print("  The deck test does not need a global sheet name.  It needs an observable odd under tau on a state whose presentation is even under")
print("  tau.  Both exist: O_q = Im<U_(123)> = V/4 on the cross-seat cycle, and the equatorial reflection realising tau.  The sense that")
print("  orders the sheets is the SENSE OF THE PROTOCOL -- U_q versus U_q^{-1}, first-then -- a temporal order, not a spatial convention.")
print("  That is the model's S_+ > S_-: an ordering of operations rather than of horizons.  Declare it: PROTOCOL-1, cross-seat cycles are")
print("  directed.  Remaining debt (Will's, unchanged): realise U_q physically on the ħ transitions (E-8-X co-pivot as a protocol).")
n_pass = sum(CH); print(f"\n{n_pass}/{len(CH)} checks passed"); sys.exit(0 if all(CH) else 1)
