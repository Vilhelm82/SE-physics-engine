#!/usr/bin/env python3
# PRED1-OPERATIONAL-CHAIN (Will's closure of the U_q debt):  X --pivoted observation--> (P_1, P_2, P_3) --directed cycle test--> Im B.
# The first arrow is seated-root physics: at three ORDERED pivots the observer records the presentation of the aligned axis,
# rho_X(i) = P(Pi_{s,theta_i}(X)) = P_i, into an ordered register; X -> X absolutely (nothing moves).  The second arrow is readout.
# Two requirements are made explicit and checked: (R1) the registers hold the presentations in ONE common frame -- the observer's own,
# which the pivot does not rotate (one root, one frame); (R2) the temporal order of observation is the direction of the test.
import sys, itertools, numpy as np
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.array([[1, 0], [0, -1]], complex); I2 = np.eye(2)
def unit(v): return v/np.linalg.norm(v)
def sdot(a): return a[0]*X + a[1]*Y + a[2]*Z
def P_of(a): return (I2 + sdot(a))/2
def kron(*ms):
    out = np.array([[1.0+0j]])
    for m in ms: out = np.kron(out, m)
    return out
def perm_op(perm):
    U = np.zeros((8, 8), complex)
    for bits in itertools.product((0, 1), repeat=3):
        src = bits[0]*4 + bits[1]*2 + bits[2]; new = tuple(bits[perm[i]] for i in range(3)); U[new[0]*4 + new[1]*2 + new[2], src] = 1
    return U
UB = perm_op((1, 2, 0))
def rot(n, th): return np.cos(th/2)*I2 - 1j*np.sin(th/2)*sdot(n)
def circumaxis(a): return unit((np.cross(a[0], a[1]) + np.cross(a[1], a[2]) + np.cross(a[2], a[0]))/2)
def read(regs): return np.trace(UB @ kron(*regs)).imag
rng = np.random.default_rng(8)
frames = [[unit(v) for v in rng.normal(size=(3, 3))] for _ in range(200)]
print("=== O-1  the chain: three ordered presentations of one source, common frame, directed test ===")
ok_a = ok_b = True
for a in frames:
    V = np.dot(a[0], np.cross(a[1], a[2])); regs = [P_of(a[i]) for i in range(3)]
    ok_a &= abs(read(regs) - V/4) < 1e-12
    Rg = rot(unit(rng.normal(size=3)), rng.uniform(0, 2*np.pi))                        # the observer's frame is a gauge if COMMON
    ok_b &= abs(read([Rg @ r @ Rg.conj().T for r in regs]) - V/4) < 1e-12
check("O-1a registers rho_X(i) = P(a_i) in the observer's frame: Im Tr[U_B x_i rho_X(i)] = V(X)/4 (200 frames)", ok_a)
check("O-1b a common rotation of all three registers leaves the reading unchanged: the observer's frame is a gauge, only the RELATIVE"
      " frame across the registers is physical (U_B commutes with R x R x R)", ok_b)
print("=== O-2  requirement R1: one root, one frame ===")
ok_c = ok_d = True
for a in frames[:80]:
    V = np.dot(a[0], np.cross(a[1], a[2])); n = circumaxis(a)
    # WRONG storage: each presentation recorded in its own pivoted frame (the pivot rotor at theta_i applied to register i)
    ths = sorted(rng.uniform(0, 2*np.pi, 3)); Rs = [rot(n, t) for t in ths]
    regs_piv = [Rs[i] @ P_of(a[i]) @ Rs[i].conj().T for i in range(3)]
    ok_c &= abs(read(regs_piv) - V/4) > 1e-6
    # WRONG storage, extreme case: each pivot brings the aligned axis onto the line of sight; the registers are then identical
    ok_d &= abs(read([P_of(n)]*3)) < 1e-12
check("O-2a if each presentation is stored in ITS OWN pivoted frame the reading is not V/4 (80/80 frames differ): the pivot must select which"
      " axis is presented without rotating the frame the registers share", ok_c)
check("O-2b extreme case -- each pivot puts the aligned axis on the line of sight, so the three registers are identical -- reads exactly 0:"
      " three copies of one presentation carry no relational obstruction", ok_d)
print("=== O-3  the deck through the chain ===")
ok_e = True
for a in frames[:80]:
    regs = [P_of(a[i]) for i in range(3)]; regs_t = [P_of(-a[i]) for i in range(3)]
    Gam = np.array([[2*np.trace(regs[i] @ regs[j]).real - 1 for j in range(3)] for i in range(3)])
    Gam_t = np.array([[2*np.trace(regs_t[i] @ regs_t[j]).real - 1 for j in range(3)] for i in range(3)])
    ok_e &= np.allclose(Gam, Gam_t) and abs(read(regs) + read(regs_t)) < 1e-12
check("O-3a Gamma(tau X) = Gamma(X) (every pairwise reading of the registers agrees) and O(tau X) = -O(X): the source is fixed, the ordinary"
      " presentation is indistinguishable, the operation exposes the transition-level difference", ok_e)
print("=== O-4  requirement R2: the temporal order of observation is the direction ===")
ok_f = ok_g = True
for a in frames[:80]:
    V = np.dot(a[0], np.cross(a[1], a[2])); regs = [P_of(a[i]) for i in range(3)]
    ok_f &= abs(read(regs[::-1]) + V/4) < 1e-12                                        # observe in the reverse order -> registers reversed
    for perm in [(1, 2, 0), (2, 0, 1)]:                                                 # cyclic re-orderings of the observation: same reading
        ok_g &= abs(read([regs[p] for p in perm]) - V/4) < 1e-12
check("O-4a observing the three alignments in the REVERSE temporal order feeds the registers reversed and reads -V/4: the direction of the"
      " test is inherited from the order of observation (PROTOCOL-1 as a fact about repeated observation, not an added axiom)", ok_f)
check("O-4b a cyclic re-ordering of the observations (starting at a different alignment) reads the same V/4: only the cyclic sense matters", ok_g)
print("=== VERDICT ===")
print("  The chain closes the U_q operational debt at the theoretical level, on two model axioms that are already in the record: one root")
print("  with one frame (the pivot selects the presented axis, it does not rotate the registers' frame -- O-2), and X -> X absolutely (the")
print("  source is unchanged by observation, so three ordered presentations of it can be taken -- O-1/O-3).  The direction of the test is")
print("  the temporal order of the observations (O-4).  What remains is apparatus wiring: repeated runs or stored temporal modes feeding")
print("  three ordered registers.  Not a preparation principle.")
n_pass = sum(CH); print(f"\n{n_pass}/{len(CH)} checks passed"); sys.exit(0 if all(CH) else 1)
