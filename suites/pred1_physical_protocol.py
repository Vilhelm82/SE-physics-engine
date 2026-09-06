#!/usr/bin/env python3
# PRED1-PHYSICAL-PROTOCOL -- Will's three physical readings of the directed permutation (controlled-U_B Hadamard test; spin-chirality
# operator; three-photon tritter), verified numerically.  Status: physical realisations of the WITNESS (ordered three-body operation),
# admissible in the record only through the preparation-interface declaration (see verdict).  Exit 0 iff all checks pass.
import sys, itertools, numpy as np
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2)
def unit(v): return v/np.linalg.norm(v)
def sdot(a): return a[0]*X + a[1]*Y + a[2]*Z
def coh(a):
    th = np.arccos(np.clip(a[2], -1, 1)); ph = np.arctan2(a[1], a[0]); return np.array([np.cos(th/2), np.exp(1j*ph)*np.sin(th/2)])
def kron(*ms):
    out = np.array([[1.0+0j]])
    for m in ms: out = np.kron(out, m)
    return out
def perm_op(perm):                                            # U|x1 x2 x3> = |x_perm(1) x_perm(2) x_perm(3)>: slot i receives old slot perm[i]
    U = np.zeros((8, 8), complex)
    for bits in itertools.product((0, 1), repeat=3):
        src = bits[0]*4 + bits[1]*2 + bits[2]; new = tuple(bits[perm[i]] for i in range(3)); U[new[0]*4 + new[1]*2 + new[2], src] = 1
    return U
UB = perm_op((1, 2, 0))                                       # U_B|x1,x2,x3> = |x2,x3,x1>
def bargmann(a):
    S = 1 + np.dot(a[0], a[1]) + np.dot(a[1], a[2]) + np.dot(a[2], a[0]); V = np.dot(a[0], np.cross(a[1], a[2])); return (S + 1j*V)/4, S, V
rng = np.random.default_rng(5)
frames = [[unit(v) for v in rng.normal(size=(3, 3))] for _ in range(200)]
print("=== H  controlled-U_B Hadamard test ===")
ok_tr = ok_y = ok_x = ok_inv = True
for a in frames:
    P = [(I2 + sdot(v))/2 for v in a]; B, S, V = bargmann(a)
    ok_tr &= abs(np.trace(UB @ kron(*P)) - np.trace(P[0] @ P[1] @ P[2])) < 1e-12 and abs(np.trace(P[0] @ P[1] @ P[2]) - B) < 1e-12
    psi = kron(coh(a[0]).reshape(2, 1), coh(a[1]).reshape(2, 1), coh(a[2]).reshape(2, 1)).flatten()
    anc = np.array([1, 1])/np.sqrt(2)
    state = np.kron(np.array([1, 0]), psi) + np.kron(np.array([0, 1]), UB @ psi); state = state/np.sqrt(2)       # controlled-U on |+>|psi>
    rho_anc = np.einsum('ij,kj->ik', state.reshape(2, 8), state.reshape(2, 8).conj())
    ok_y &= abs(np.trace(rho_anc @ Y).real - V/4) < 1e-12
    ok_x &= abs(np.trace(rho_anc @ X).real - S/4) < 1e-12
    state2 = (np.kron(np.array([1, 0]), psi) + np.kron(np.array([0, 1]), UB.conj().T @ psi))/np.sqrt(2)
    rho2 = np.einsum('ij,kj->ik', state2.reshape(2, 8), state2.reshape(2, 8).conj())
    ok_inv &= abs(np.trace(rho2 @ Y).real + V/4) < 1e-12
check("H-1 Tr[U_B (P_1 x P_2 x P_3)] = Tr(P_1 P_2 P_3) = B exactly (200 frames)", ok_tr)
check("H-2 ancilla |+>, controlled-U_B, measure Y: <Y> = Im B = V/4;  X: <X> = Re B = (1 + sum gamma)/4", ok_y and ok_x)
check("H-3 controlled-U_B^dagger reads <Y> = -V/4: U_B is the transition word 3->2->1 (reads B), U_B^dagger is 1->2->3 (reads conj B) --"
      " matches pred1_protocol T-2c", ok_inv)
print("=== C  the spin-chirality operator ===")
def three(A, Bm, C): return kron(A, Bm, C)
chi = (three(X, Y, Z) + three(Y, Z, X) + three(Z, X, Y) - three(X, Z, Y) - three(Y, X, Z) - three(Z, Y, X))/(2*np.sqrt(3))
ev = np.sort(np.linalg.eigvalsh(chi).round(9))
check("C-1 chi-hat has spectrum {-1 (x2), 0 (x4), +1 (x2)}: +-1 on the two mixed-symmetry doublets, 0 on the symmetric quadruplet",
      np.allclose(ev, [-1, -1, 0, 0, 0, 0, 1, 1]), f"{ev}")
ok_c = True
for a in frames[:60]:
    psi = kron(coh(a[0]).reshape(2, 1), coh(a[1]).reshape(2, 1), coh(a[2]).reshape(2, 1)).flatten(); B, S, V = bargmann(a)
    ok_c &= abs((psi.conj() @ chi @ psi).real - (2/np.sqrt(3))*V/4) < 1e-12
check("C-2 <chi-hat> = (2/sqrt 3) Im B on the product state: chi-hat = (1/2 sqrt 3) sigma_1.(sigma_2 x sigma_3) reads V/(2 sqrt 3)", ok_c)
w, Vc = np.linalg.eigh(chi); expo = Vc @ np.diag(np.exp(1j*(2*np.pi/3)*w)) @ Vc.conj().T
which = 'U_B' if np.allclose(expo, UB) else ('U_B^dagger' if np.allclose(expo, UB.conj().T) else 'neither')
check("C-3 exp(+i (2pi/3) chi-hat) equals U_B or U_B^dagger (the 3-cycle acts as e^{+-2pi i/3} on the doublets, 1 on the quadruplet)",
      which != 'neither', f"exp(+i 2pi/3 chi-hat) = {which}")
print("=== T  three-photon tritter, internal states = the spinor roots ===")
om = np.exp(2j*np.pi/3); Ut = np.array([[1, 1, 1], [1, om, om**2], [1, om**2, om]])/np.sqrt(3)
def tritter_probs(a):
    # first quantisation: photon = external mode (3) x internal (2); bosonic symmetrised 3-photon state; apply U on externals
    phis = [np.kron(np.eye(3)[j], coh(a[j])) for j in range(3)]                       # input: photon j in mode j, internal |a_j>
    psi = np.zeros(6**3, complex)
    for perm in itertools.permutations(range(3)):
        psi += kron(*[phis[perm[k]].reshape(6, 1) for k in range(3)]).flatten()
    psi /= np.sqrt(6)                                                                  # orthogonal input modes -> normalised
    Uext = np.kron(Ut, I2); psi_out = kron(Uext, Uext, Uext) @ psi
    P = {}
    for occ in [(3, 0, 0), (0, 3, 0), (0, 0, 3), (2, 1, 0), (1, 2, 0), (2, 0, 1), (1, 0, 2), (0, 2, 1), (0, 1, 2), (1, 1, 1)]:
        tot = 0.0
        for modes in itertools.product(range(3), repeat=3):
            if tuple(modes.count(k) for k in range(3)) != occ: continue
            proj = kron(*[np.kron(np.outer(np.eye(3)[m], np.eye(3)[m]), I2) for m in modes])
            tot += (psi_out.conj() @ proj @ psi_out).real
        P[occ] = tot
    return P
ok_t = ok_norm = ok_deck = True; other_pairs = {}; signs_210 = set()
for a in frames[:40]:
    P = tritter_probs(a); B, S, V = bargmann(a)
    ok_norm &= abs(sum(P.values()) - 1) < 1e-10
    c210 = (3*np.sqrt(3)/2)*(P[(2, 1, 0)] - P[(1, 2, 0)])*4/V if abs(V) > 1e-3 else np.nan
    c201 = (3*np.sqrt(3)/2)*(P[(2, 0, 1)] - P[(1, 0, 2)])*4/V if abs(V) > 1e-3 else np.nan
    ok_t &= (np.isnan(c210) or abs(abs(c210) - 1) < 1e-9) and (np.isnan(c201) or abs(c201 - 1) < 1e-9)
    signs_210 = signs_210 | {round(c210)} if not np.isnan(c210) else signs_210
    Pd = tritter_probs([-v for v in a])                                                  # the deck-reflected triple: same Gram, V -> -V
    ok_deck &= abs((P[(2, 1, 0)] - P[(1, 2, 0)]) + (Pd[(2, 1, 0)] - Pd[(1, 2, 0)])) < 1e-9 and abs(P[(1, 1, 1)] - Pd[(1, 1, 1)]) < 1e-9
    for pair in [((2, 0, 1), (1, 0, 2)), ((0, 2, 1), (0, 1, 2))]:
        other_pairs.setdefault(pair, []).append((3*np.sqrt(3)/2)*(P[pair[0]] - P[pair[1]])*4/V if abs(V) > 1e-3 else np.nan)
check("T-1 the ten output patterns sum to 1 (bosonic symmetrisation and the projectors are consistent)", ok_norm)
check("T-2 Im B = +-(3 sqrt 3 / 2)(P_(210) - P_(120)) EXACTLY with a frame-independent sign fixed by the tritter's chirality (omega vs conj omega);"
      " with U_jk = omega^{(j-1)(k-1)} the pair (2,0,1)-(1,0,2) reads +V/4 with coefficient exactly 1 -- the device's handedness is the directed"
      " reference here, as the permutation order is in the Hadamard test", ok_t and len(signs_210) == 1, f"sign of the (210)-(120) coefficient: {signs_210}")
check("T-3 the deck-reflected triple (same Gram, V -> -V) reverses the imbalance and leaves P_(111) unchanged: pairwise HOM content is even,"
      " the imbalance is the odd datum", ok_deck)
print("    other imbalance pairs, coefficient x4/V:", {k: (round(np.nanmin(v), 4), round(np.nanmax(v), 4)) for k, v in other_pairs.items()})
print("=== G  what pairwise tests give ===")
ok_g = True
for a in frames[:60]:
    P = [(I2 + sdot(v))/2 for v in a]; B, S, V = bargmann(a)
    Gam = np.array([[2*np.trace(P[i] @ P[j]).real - 1 for j in range(3)] for i in range(3)])
    ok_g &= np.allclose(Gam, [[np.dot(a[i], a[j]) for j in range(3)] for i in range(3)]) and abs(np.linalg.det(Gam) - V**2) < 1e-10 \
            and abs(np.sqrt(np.linalg.det(Gam))/4 - abs(B.imag)) < 1e-10
check("G-1 Gamma_ij = 2 Tr(P_i P_j) - 1 = gamma_ij (pairwise SWAP/HOM), det Gamma = V^2, so pairwise tests give |Im B| = sqrt(det Gamma)/4;"
      " only the directed three-body operation supplies the sign", ok_g)
print("=== VERDICT ===")
print("  All three readings are exact.  Each is a physical realisation of the ordered three-body operation Tr(U_B P_1 x P_2 x P_3), which")
print("  pivot_map.py showed the record's bilinear map does NOT contain.  Their admissibility therefore rests on ONE declaration: that the")
print("  seat's spinor roots |a_i> are physically preparable hbar-qubit inputs and a directed three-body operation on them is a model")
print("  operation (Will's boundary).  With it, the protocol debt is closed physically; without it, the debt is the preparation interface")
print("  and nothing else -- not the map C -> K, not the readout.  Pairwise HOM/SWAP is the Gram; the sign is three-body and directed.")
n_pass = sum(CH); print(f"\n{n_pass}/{len(CH)} checks passed"); sys.exit(0 if all(CH) else 1)
