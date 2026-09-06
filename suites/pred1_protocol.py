#!/usr/bin/env python3
# PROTOCOL-1 -- the seat-cycle read as a physical protocol, in the record's own objects (pred1_protocol.py).
# Objects allowed: clicks P_a = (1 + a)/2 (idempotents of Cl(3), Pauli rep a -> a.sigma); spinor roots psi_a with P_a = |psi_a><psi_a|
#   (SECT-1 S-5, E-8-X X-1); the co-pivot rotor R_{a->b} = (1 + b a)/|1 + b a| acting two-sided on the axis and one-sided on the spinor
#   (E-8-X design input); the transition pairing <psi|psi'> with its I-part (E-8-X X-6).  PRED-1 banked B = (S + iV)/4, J = V.
# Two directed protocols, both reading the orientation V with the sense of the sequence:
#   (K) the CLICK protocol: apply the three clicks in temporal order; the composite click's pseudoscalar part is V/4;
#   (P) the PIVOT protocol: carry the spinor root around the cycle by the three geodesic co-pivots; it returns with phase +-Omega/2.
# BANNED until the comparison block: Pancharatnam, Berry, solid angle by name.  Exit 0 iff every check passes.
import sys, itertools, numpy as np, sympy as sp
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
sig = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
I2 = np.eye(2)
def sdot(a): return sum(a[k]*sig[k] for k in range(3))
def click(a): return (I2 + sdot(a))/2
def root(a):                                              # spinor root: the +1 eigenvector of a.sigma, phase arbitrary (lift)
    w, v = np.linalg.eigh(sdot(a)); return v[:, np.argmax(w)]
def copivot(a, b): return (I2 + sdot(b) @ sdot(a))/np.sqrt(2*(1 + np.dot(a, b)))
def bargmann(a):
    S = 1 + np.dot(a[0], a[1]) + np.dot(a[1], a[2]) + np.dot(a[2], a[0]); V = np.dot(a[0], np.cross(a[1], a[2])); return (S + 1j*V)/4, S, V
def unit(v): return v/np.linalg.norm(v)
rng = np.random.default_rng(12)
frames = [[unit(v) for v in rng.normal(size=(3, 3))] for _ in range(200)]
print("=== T-1  the CLICK protocol (symbolic, general axes) ===")
A = [sp.Matrix(sp.symbols(f'a{i}1 a{i}2 a{i}3', real=True)) for i in (1, 2, 3)]
sx, sy, sz = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])
def sdot_s(a): return a[0]*sx + a[1]*sy + a[2]*sz
P = [(sp.eye(2) + sdot_s(a))/2 for a in A]
S_ = 1 + A[0].dot(A[1]) + A[1].dot(A[2]) + A[2].dot(A[0]); V_ = A[0].dot(A[1].cross(A[2]))
tr123 = sp.expand((P[0]*P[1]*P[2]).trace()); tr321 = sp.expand((P[2]*P[1]*P[0]).trace())
check("T-1a Tr(P_1 P_2 P_3) = (1 + sum gamma + i V)/4 identically: the ordered product of three clicks has scalar part (1+sum gamma)/4 and"
      " PSEUDOSCALAR part V/4 -- the orientation lives in the pseudoscalar grade, 'the sand is the phase'", sp.expand(tr123 - (S_ + sp.I*V_)/4) == 0)
check("T-1b the REVERSED order P_3 P_2 P_1 (reversion) gives the conjugate: the sign of V is the temporal order of the clicks", sp.expand(tr321 - (S_ - sp.I*V_)/4) == 0)
print("=== T-2  the co-pivot rotor and the PIVOT protocol (numeric, 200 frames; exact at the trine) ===")
ok_a = ok_b = ok_c = ok_d = ok_e = True
for a in frames:
    R12, R23, R31 = copivot(a[0], a[1]), copivot(a[1], a[2]), copivot(a[2], a[0])
    ok_a &= np.allclose(R12 @ sdot(a[0]) @ R12.conj().T, sdot(a[1])) and np.allclose(R12 @ R12.conj().T, I2)   # two-sided: a -> b, unitary
    H = R31 @ R23 @ R12                                                                                        # temporal order 1->2, 2->3, 3->1
    ok_b &= np.allclose(H @ sdot(a[0]) @ H.conj().T, sdot(a[0]))                                                 # the axis returns
    psi = root(a[0]); ph = psi.conj() @ H @ psi                                                                   # one-sided: the spinor returns with a phase
    B, S, V = bargmann(a)
    ok_c &= abs(abs(ph) - 1) < 1e-10 and abs(ph - np.conj(B)/abs(B)) < 1e-10                                    # temporal 1->2->3 reads conj(B)/|B|: the record's B = Tr(P1 P2 P3) is rightmost-first, i.e. temporal 3->2->1
    K123 = np.trace(click(a[2]) @ click(a[1]) @ click(a[0]))                                                       # the CLICK protocol in temporal order 1->2->3
    ok_c &= abs(K123 - np.conj(B)) < 1e-10 and abs(K123/abs(K123) - ph) < 1e-10                                    # same temporal order => same phase as the pivot protocol
    Hr = R12.conj().T @ R23.conj().T @ R31.conj().T                                                              # the reversed cycle
    ok_d &= abs(psi.conj() @ Hr @ psi - B/abs(B)) < 1e-10
    ta = [-v for v in a]; Ht = copivot(ta[2], ta[0]) @ copivot(ta[1], ta[2]) @ copivot(ta[0], ta[1]); pt = root(ta[0])
    ok_e &= abs((pt.conj() @ Ht @ pt).imag + ph.imag) < 1e-10                                                    # the deck flips the I-part
check("T-2a R_{a->b} = (1 + b a)/sqrt(2(1+a.b)) is unitary and carries a to b two-sided (E-8-X's co-pivot)", ok_a)
check("T-2b the cycle H = R_{3->1} R_{2->3} R_{1->2} returns the axis a_1 to itself", ok_b)
check("T-2c one-sided on the spinor root, H psi_1 = e^{i chi} psi_1 with e^{i chi} = conj(B)/|B| EXACTLY for the temporal order 1->2->3, and the"
      " CLICK protocol in the same temporal order, Tr(P_3 P_2 P_1), reads conj(B) -- SAME SEQUENCE, SAME PHASE.  (PRED-1's B = Tr(P_1 P_2 P_3)"
      " is the rightmost-first operator product, i.e. temporal order 3->2->1.)  Pivot: unit modulus; click: modulus |B|; E-8-X Y-1 links them", ok_c)
check("T-2d the reversed cycle (temporal 1->3->2) returns B/|B|: reversing the sequence conjugates the reading", ok_d)
check("T-2e the deck (a_i -> -a_i, Gram fixed) flips the I-part of the returned phase under the same sequence", ok_e)
e = [np.eye(3)[k] for k in range(3)]
Ht = copivot(e[2], e[0]) @ copivot(e[1], e[2]) @ copivot(e[0], e[1]); pt = root(e[0]); pht = pt.conj() @ Ht @ pt
check("T-2f exact at the trine: the returned phase is -pi/4 for temporal 1->2->3 (= -arg B, half the octant with the sequence's sign)", abs(pht - np.exp(-1j*np.pi/4)) < 1e-12, f"phase {np.angle(pht):.6f}")
print("=== T-3  lift independence, click-equivariance ===")
ok_f = ok_g = True
for a in frames[:60]:
    psi = root(a[0]); H = copivot(a[2], a[0]) @ copivot(a[1], a[2]) @ copivot(a[0], a[1])
    for phz in (0.3, 1.7, -2.2):
        psi2 = np.exp(1j*phz)*psi; ok_f &= abs(psi2.conj() @ H @ psi2 - psi.conj() @ H @ psi) < 1e-12
    ph0 = psi.conj() @ H @ psi
    for perm in itertools.permutations(range(3)):
        pa = [a[perm[i]] for i in range(3)]
        # the relabelled protocol visits the same physical axes in the same temporal order: start at pa[inv 0] ...
        inv = [perm.index(i) for i in range(3)]; seq = [pa[inv[0]], pa[inv[1]], pa[inv[2]]]
        Hp = copivot(seq[2], seq[0]) @ copivot(seq[1], seq[2]) @ copivot(seq[0], seq[1]); pp = root(seq[0])
        ok_g &= abs(pp.conj() @ Hp @ pp - ph0) < 1e-10
check("T-3a the reading is independent of the lift psi_1 -> e^{i phi} psi_1 (the one-sided phase is unobservable; only the cycle's phase is read)", ok_f)
check("T-3b relabel the axes and conjugate the sequence with them: the reading is unchanged for all six relabellings (click-equivariance)", ok_g)
print("=== T-4  the readout: pair the returned spinor with the kept one against a reference delay (symbolic) ===")
d, chi = sp.symbols('delta chi', real=True)
Iv = sp.expand((1 + sp.exp(sp.I*(d + chi)))*(1 + sp.exp(-sp.I*(d + chi)))).rewrite(sp.cos)
check("T-4a I(delta) = |1 + e^{i delta} e^{i chi}|^2 = 2 + 2 cos(delta + chi): the fringe is shifted by -chi; its DIRECTION is read against the"
      " sense of the delay, which is temporal (later arrival).  The sign of V is therefore read against 'later', never against 'left'",
      sp.simplify(Iv - (2 + 2*sp.cos(d + chi))) == 0)
print("=== COMPARISON BLOCK (names allowed) ===")
print("  The click protocol is Pancharatnam's three-projection cycle (1956); the pivot protocol is the holonomy of three geodesic rotations")
print("  (Berry 1987); both give half the solid angle of the triangle, which is E-8-X Y-1 and CENSUS-C C-1a.  Physical carrier at hand:")
print("  polarisation.  The spinor root is the Jones vector, the axis is the Stokes vector, SECT-1 D-3's angle doubling IS the Stokes")
print("  construction, a click is a polariser, a co-pivot is a rotator/waveplate, the pairing is interference.  Everything read is a")
print("  field intensity: classical optics, no Born rule and no hbar anywhere -- what the hbar SEAT contributes is the one-sided")
print("  REPRESENTATION, not a quantum of action.  The orientation V of three polarisation states on the Poincare sphere is invisible to")
print("  any single-seat (Stokes/Gram) reading and visible to the directed cycle.  That is PREDICTION-1 with its protocol assigned.")
print("=== PROTOCOL-1 (to declare) ===")
print("  Cross-seat cycles are DIRECTED: a protocol is a temporal sequence of clicks or co-pivots; reversing the sequence conjugates the")
print("  reading.  The sign of the seat-cycle phase is defined relative to the sequence and read against a delay.  This is the model's")
print("  ordering of the deck sheets: an order of operations, not of horizons.")
n_pass = sum(CH); print(f"\n{n_pass}/{len(CH)} checks passed"); sys.exit(0 if all(CH) else 1)
