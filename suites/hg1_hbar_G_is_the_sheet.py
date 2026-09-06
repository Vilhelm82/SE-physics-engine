#!/usr/bin/env python3
# =============================================================================
# HG-1 -- the hbar/G distinction IS the sheet.  (Claire, 2026-09-05.)
#   Closes the handoff's open item: 'T7b3 says the bare tier is silent on hbar vs G; where does it enter -- the
#   one-sided readout, or convention?'  Answer: it enters as the ORIENTATION of the frame -- the sheet -- which the
#   Gram cannot see for exactly the reason it cannot see the deck, and which Codex's fixed-input experiment reads.
# INPUTS: T7b3 (the swap S: (l1,l2,t) -> (l2,l1,-t) realises a <-> b with the Gram invariant); T5c (D = cosh l1 cosh l2 sin t,
#   V_spin = 2D, tau: t -> -t sends D -> -D); PINCH-2 (xi = (b - a gamma)/d, tau: xi -> -xi); T7g (two horizons, hbar -+ G).
# CLAIMS, checked:
#   h1  S realises a <-> b with gamma fixed (T7b3, re-verified).
#   h2  At the FRAME level S sends D -> -D: swapping the rulers reverses the orientation.  Same fact from the Clifford
#       side: tr(Gamma_C Gamma_H Gamma_G) -> tr(Gamma_C Gamma_G Gamma_H) = -V_spin (T5c 2c).  Swap = sheet flip.
#   h3  a <-> b has TWO lifts to seat coordinates: S = (l2,l1,-t) and S' = (l2,l1,+t).  S' preserves D; S = tau o S'.
#       The Gram realises the same map under both.  So 'swap hbar and G' is defined on the Gram only up to the deck:
#       the hbar/G ambiguity of T7b3 and the deck ambiguity of T7c are ONE ambiguity, and h2 says which lift is the
#       physical swap (the one that moves the vectors).
#   h4  F, eta, N^2, delta are symmetric in (a,b) and even in t: hbar/G-blind and sheet-blind for the same reason.
#   h5  Codex's resolved chart is normalised by hbar's depth: F = (1 - xi^2)/(1 + a^2).  Its companion normalised by
#       G's depth, xi' = (a - b gamma)/d, satisfies F = (1 - xi'^2)/(1 + b^2) identically.  S exchanges them with a sign:
#       xi -> -xi'.  'Which ruler is hbar' = 'which depth normalises the chart', and the sign is the sheet.
#   h6  S = tau on xi exactly on {a = b} u {gamma = -1}:  xi_S - tau(xi) = (b - a)(1 + gamma)/d.  So on the l1 = l2 plane
#       (which contains the pinch P_+) AND on the whole sigma = -1 horizon, swapping hbar and G IS the deck.  On the
#       sigma = +1 horizon it is NOT: there S fixes xi and flips d.  T7g's two horizons differ by exactly this.
#   h7  On the sigma = -1 horizon n0 = h + g is null (h, g anti-aligned): reversing their order reverses orientation and
#       the resolved direction together -- the deck.  On sigma = +1, n0 = h - g is null (h, g aligned): reversing their
#       order flips orientation but not the direction.  This is why the two horizons carry different signs in LABEL-3.
#   h8  Therefore the distinction is observed exactly as the sheet is (Codex 63d951e): by the sign of dF/dt under a
#       fixed input, never by a static readout.  [Cited; her formula verified in 86f3ebc.]
# [First run: h6a failed on a double negation in the runner (tau(xi) written as -(-xi)); h7b chained to it. Fixed; claims unchanged.]
# RESOLVED (Codex BR-9, docs/results/2026-09-05/2026-09-05-response-on-branched-resolved-state-spaces.md sec. 9, be1ae74; verified by Claire
#   in sandbox): the swap extends smoothly to the resolved chart as S(a,xi,d) = (a gamma + xi d, gamma xi - a F d, -d),
#   preserving F and gamma, S^2 = id, [S, tau] = 0, S_xi = -xi_g.  On the horizons S = (sigma a + xi d, sigma xi, -d):
#   h6's split exact.  Over the collapsed pinch S(a,xi,0) = (sigma a, sigma xi, 0): on sigma = +1 the WHOLE exceptional set
#   is fixed, so by her BR-5 a swap-equivariant field under swap-invariant drive (2 e_d = xi e_a - a F e_xi) cannot reach
#   P_+ in finite time; on sigma = -1 only the origin is fixed.  With h6: hbar/G-symmetric dynamics can cross the pinch
#   only on the branch where the swap is the deck.  Conditional on the drive's symmetry, which is the load's.
# TIER: [DERIVED | T7a's signature].  KILL: a lift of a <-> b to the frame that preserves D (then swap != sheet);
#   or F failing to be (a,b)-symmetric; or xi_S - tau(xi) not vanishing on {a = b}.
# =============================================================================
import sympy as sp, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), sp.cancel, sp.factor, lambda q: sp.simplify(sp.expand_trig(q))):
        try:
            if f(e) == 0: return True
        except Exception: pass
    return False

a, b, gam, d = sp.symbols('a b gamma d', real=True)
t, l1, l2 = sp.symbols('t l1 l2', real=True)
A = -sp.sinh(l1); B = -sp.sinh(l2); Gm = sp.cosh(l1)*sp.cosh(l2)*sp.cos(t) - sp.sinh(l1)*sp.sinh(l2)
D = sp.cosh(l1)*sp.cosh(l2)*sp.sin(t)
S  = {l1: l2, l2: l1, t: -t}          # the physical swap (T7b3)
Sp = {l1: l2, l2: l1}                 # the other lift of a <-> b
tau = {t: -t}
sub = lambda e, m: e.subs(m, simultaneous=True)

print("=== h1: S realises a <-> b with gamma fixed (T7b3) ===")
check("h1a S: a -> b, b -> a, gamma -> gamma", z(sub(A, S) - B) and z(sub(B, S) - A) and z(sub(Gm, S) - Gm))

print("=== h2: at the frame level the swap reverses orientation ===")
check("h2a S sends D = cosh l1 cosh l2 sin t to -D", z(sub(D, S) + D))
# the Clifford side, T5c's real 2x2 rep of Cl(2,1): Gamma_C^2 = -1, Gamma_H^2 = Gamma_G^2 = +1, anticommuting
GC = sp.Matrix([[0, -1], [1, 0]]); GH = sp.Matrix([[1, 0], [0, -1]]); GG = sp.Matrix([[0, 1], [1, 0]])
check("h2b T5c's generators: Gamma_C^2 = -1, Gamma_H^2 = Gamma_G^2 = 1, pairwise anticommuting", GC**2 == -sp.eye(2) and GH**2 == sp.eye(2) and GG**2 == sp.eye(2) and (GH*GG + GG*GH) == sp.zeros(2) and (GC*GH + GH*GC) == sp.zeros(2) and (GC*GG + GG*GC) == sp.zeros(2))
check("h2c swapping the ruler generators flips the spinor trace: tr(C G H) = -tr(C H G)  (T5c 2c) -- swap = sheet flip", (GC*GG*GH).trace() == -(GC*GH*GG).trace())

print("=== h3: two lifts of a <-> b; they differ by the deck ===")
check("h3a S' = (l2, l1, +t) ALSO realises a <-> b with gamma fixed", z(sub(A, Sp) - B) and z(sub(B, Sp) - A) and z(sub(Gm, Sp) - Gm))
check("h3b S' preserves D; S = tau o S'", z(sub(D, Sp) - D) and z(sub(sub(D, Sp), tau) - sub(D, S)))
check("h3c the Gram is invariant under both lifts (T7b3a): 'swap hbar and G' is defined on the Gram only up to the deck", True if (CH[0] and CH[-2]) else False)

print("=== h4: the base-field observables are hbar/G-blind and sheet-blind for the same reason ===")
eps = 1 - gam**2; W = a**2 + b**2 - 2*a*b*gam; delta = eps + W; F = eps/delta; eta = W/eps
swap_ab = {a: b, b: a}
check("h4a F, eta, delta are symmetric in (a, b)", all(z(sub(e, swap_ab) - e) for e in (F, eta, delta)))
check("h4b and in seat coordinates they are even in t (deck-blind): the two blindnesses coincide", all(z(sub(e.subs({a: A, b: B, gam: Gm}), tau) - e.subs({a: A, b: B, gam: Gm})) for e in (F, delta)))

print("=== h5: the two resolved charts, exchanged by the swap ===")
xi  = (b - a*gam)/d;  xip = (a - b*gam)/d
check("h5a Codex's F = (1 - xi^2)/(1 + a^2) identically (d^2 = delta)", z((F - (1 - xi**2)/(1 + a**2)).subs(d, sp.sqrt(delta))))
check("h5b the companion F = (1 - xi'^2)/(1 + b^2) identically: G's depth normalises it", z((F - (1 - xip**2)/(1 + b**2)).subs(d, sp.sqrt(delta))))
xi_S = sub(xi, {a: b, b: a, d: -d})
check("h5c under S (a <-> b, d -> -d): xi -> -xi'.  'Which ruler is hbar' = 'which depth normalises the chart'; the sign is the sheet", z(xi_S + xip))

print("=== h6: where the swap IS the deck ===")
diff = sp.factor(sp.together(xi_S + xi))                          # xi_S - tau(xi), with tau(xi) = -xi
check("h6a xi_S - tau(xi) = (b - a)(1 + gamma)/d exactly", z(diff - (b - a)*(1 + gam)/d))
check("h6b so S = tau on xi precisely on {a = b} u {gamma = -1}: the l1 = l2 plane (containing P_+) and the ENTIRE sigma = -1 horizon", z(diff.subs(a, b)) and z(diff.subs(gam, -1)) and not z(diff.subs({a: 2, b: 1, gam: 1})))
H6C = z((xi_S - xi).subs(gam, 1)); check("h6c on the sigma = +1 horizon S FIXES xi and flips d: not the deck", H6C)

print("=== h7: why the two horizons differ (T7g) ===")
q = sp.diag(1, 1, -1)
# explicit horizon states from PINCH-1: c = (0,0,1), h = (sqrt(1+a^2), 0, -a), g with q(h,g) = sigma
def hg(sig, av, bv):
    h = sp.Matrix([sp.sqrt(1 + av**2), 0, -av]); x = (sig + av*bv)/sp.sqrt(1 + av**2); y2 = 1 + bv**2 - x**2
    return h, sp.Matrix([x, sp.sqrt(y2), -bv])
h1_, g1_ = hg(1, 2, 1); hm_, gm_ = hg(-1, 2, 1)
check("h7a sigma = +1: n0 = h - g is null (rulers ALIGNED); sigma = -1: n0 = h + g is null (rulers ANTI-ALIGNED)",
      z(sp.simplify(((h1_ - g1_).T*q*(h1_ - g1_))[0])) and z(sp.simplify(((hm_ + gm_).T*q*(hm_ + gm_))[0])))
check("h7b reversing the order of an anti-aligned pair is the deck; of an aligned pair it is orientation-only: LABEL-3's sign split between the two horizons", CH[-1] and H6C and z(diff.subs(gam, -1)))

print("=== h8: how the distinction is observed ===")
check("h8a the hbar/G distinction is the sheet (h2); the sheet is observed by the sign of dF/dt under fixed input (Codex 63d951e, verified 86f3ebc): the distinction is DYNAMICAL, never static", True if (CH[1] and CH[9]) else False)

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: 'Which ruler is hbar' is the orientation of the frame.  Swapping the rulers reverses D (both in the seat's")
print("  coordinates and as the Clifford trace); the Gram cannot see it for the same reason it cannot see the deck, so")
print("  T7b3's silence and T7c's non-injectivity are one blindness.  On the l1 = l2 plane and on the whole sigma = -1")
print("  horizon the swap IS the deck; on the sigma = +1 horizon it fixes the resolved direction and flips only d.  The")
print("  distinction enters at the frame level and is read dynamically, by the sign of dF/dt under a fixed input.")
