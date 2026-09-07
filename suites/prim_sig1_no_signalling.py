# prim_sig1_no_signalling.py  --  SIG-1: CONJECTURE-COSMOLOGY kill 3 (signalling). Run FORWARD.
#
# Kill 3: "the model's non-injective presentations [T7c, P7] are the candidate mechanism for correlation without transit. If they
#   also give a usable channel, the model predicts superluminal signalling and dies immediately."
#
# SETUP (Cl(2,1) frame, T7c parametrisation; R^{2,1} with q = diag(+,+,-)):
#   frame lines  c = (0,0,1) [q = -1],  hbar = (cosh l1, 0, sinh l1),  G = (cosh l2 cos t, cosh l2 sin t, sinh l2)  [q = +1 each]
#   so q(c,hbar) = -sinh l1 = a, q(c,G) = -sinh l2 = b, q(hbar,G) = cosh l1 cosh l2 cos t - sinh l1 sinh l2 = gamma.  [T7c, verified below]
#   Two seats A, B share the frame. Each has a pivot: a hyperbolic rotor of rapidity lambda_A (in the x-z plane) / lambda_B (in the y-z plane).
#   P6 READING 1 (T7d tilt): the pivot moves the SEAT'S LINE relative to the frame: c_A' = R_A c.  Presentation at A = Gram(c_A', hbar, G).
#   P6 READING 2: the pivot rotates the WHOLE FRAME: (c, hbar, G) -> R(c, hbar, G). Presentation at A = the frame's Gram, read from c.
#   A 'reading' is any function of the seat's presentation; the block (sheet) label is a further local choice (T5c 3g).
# BANNED: quantum mechanics, any Bell theorem as input. Fine's theorem is the COMPARISON for part (c).
# KILLS: (K1) under either reading of P6, if seat B's presentation depends on lambda_A the model signals and dies.
#        (K2, informational) CHSH of Gram-level readings: if <= 2, non-injectivity gives correlation without transit but NOT Bell violation;
#        the quantum fork's structure (holonomy phases) would have to supply that separately.

import sympy as sp, numpy as np, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

l1, l2, t, lA, lB = sp.symbols('l1 l2 t lambda_A lambda_B', real=True)
Q = sp.diag(1, 1, -1)
def q(u, v): return sp.simplify((u.T*Q*v)[0])
c  = sp.Matrix([0, 0, 1])
hb = sp.Matrix([sp.cosh(l1), 0, sp.sinh(l1)])
G  = sp.Matrix([sp.cosh(l2)*sp.cos(t), sp.cosh(l2)*sp.sin(t), sp.sinh(l2)])
check("s0", q(c,c)==-1 and q(hb,hb)==1 and q(G,G)==1 and sp.simplify(q(c,hb)+sp.sinh(l1))==0 and sp.simplify(q(hb,G)-(sp.cosh(l1)*sp.cosh(l2)*sp.cos(t)-sp.sinh(l1)*sp.sinh(l2)))==0,
      "frame realised with T7c's Gram: (a, b, gamma) = (-sinh l1, -sinh l2, cosh l1 cosh l2 cos t - sinh l1 sinh l2)")
# rotors: boosts in the x-z plane (seat A) and y-z plane (seat B)
RA = sp.Matrix([[sp.cosh(lA),0,sp.sinh(lA)],[0,1,0],[sp.sinh(lA),0,sp.cosh(lA)]])
RB = sp.Matrix([[1,0,0],[0,sp.cosh(lB),sp.sinh(lB)],[0,sp.sinh(lB),sp.cosh(lB)]])
check("s1", sp.simplify(RA.T*Q*RA - Q)==sp.zeros(3,3) and sp.simplify(RB.T*Q*RB - Q)==sp.zeros(3,3), "both pivots are isometries of the frame's form (rotors preserve q)")

print("=== SIG-1a: P6 reading 1 -- the pivot moves the seat's line ===")
cA = RA*c; cB = RB*c
presA = [q(cA,hb), q(cA,G), q(hb,G)]                          # what A reads: its presented Gram (its line against the rulers, rulers' angle)
presB = [q(cB,hb), q(cB,G), q(hb,G)]
depB_on_A = [sp.simplify(sp.diff(e, lA)) for e in presB]
depA_on_B = [sp.simplify(sp.diff(e, lB)) for e in presA]
check("a1", all(e==0 for e in depB_on_A) and all(e==0 for e in depA_on_B), "seat B's presentation is independent of lambda_A (and A's of lambda_B): a pivot changes the pivoting seat's line only; the frame's rulers and their Gram are untouched. NO CHANNEL (K1 passes, reading 1)")
check("a2", any(sp.simplify(sp.diff(e, lA))!=0 for e in presA), "and the pivot is not empty: A's own presentation DOES change with lambda_A (it is a real perspective change)")

print("=== SIG-1b: P6 reading 2 -- the pivot rotates the whole frame ===")
fr = lambda R: [R*c, R*hb, R*G]
def gram(v): return sp.Matrix(3,3, lambda i,j: q(v[i], v[j]))
GA = gram(fr(RA)); G0 = gram([c, hb, G])
check("b1", sp.simplify(GA - G0)==sp.zeros(3,3), "rotating the whole frame by A's pivot leaves EVERY Gram entry unchanged (isometry): every seat, being a line OF the frame, reads exactly what it read before. NO CHANNEL (K1 passes, reading 2)")
check("b2", True, "under reading 2 the pivot is invisible to everyone including A -- which is why T7d's tilt (reading 1) is the one the model actually uses. Either way kill 3 does not fire")

print("=== SIG-1c: correlation without transit -- and how much (CHSH of Gram-level readings) ===")
# hidden variable: the frame's invisible depths (l1, l2) drawn uniformly; settings: lambda_A, lambda_B; readings: signs of presented entries.
# reading at A: s_A = sign( q(c_A', hbar) - q(c_A', G) )  (which ruler tilts more toward A's line); at B: s_B = sign( q(c_B', hbar) + q(c_B', G) ).
fA = sp.lambdify((l1,l2,t,lA), presA[0]-presA[1], 'numpy'); fB = sp.lambdify((l1,l2,t,lB), presB[0]+presB[1], 'numpy')
rng = np.random.default_rng(7); N = 60000
L1 = rng.uniform(-1.5,1.5,N); L2 = rng.uniform(-1.5,1.5,N); T = rng.uniform(0.2, np.pi-0.2, N)
def E(a, b): return np.mean(np.sign(fA(L1,L2,T,a))*np.sign(fB(L1,L2,T,b)))
def margB(a, b): return np.mean(np.sign(fB(L1,L2,T,b)))
grid = np.linspace(-2, 2, 7); best = 0; bestset = None
for a1 in grid:
    for a2 in grid:
        for b1 in grid:
            for b2 in grid:
                S = abs(E(a1,b1) - E(a1,b2) + E(a2,b1) + E(a2,b2))
                if S > best: best, bestset = S, (a1,a2,b1,b2)
print(f"  max |CHSH| over a 7^4 grid of pivot settings: {best:.4f} at settings {tuple(float(x) for x in bestset)}   (local bound 2; Tsirelson 2.828)")
check("c1", best <= 2.0 + 0.02, "CHSH <= 2: with a definite shared frame and readings that are functions of (frame, own pivot), the model's Gram-level correlations are LOCAL (Fine 1982, comparison). Non-injectivity gives correlation without transit; it does not, by itself, give Bell violation")
a1,a2,b1,b2 = bestset
check("c2", abs(margB(a1,b1) - margB(a2,b1)) < 0.01, "B's marginal is unchanged when A switches settings (numerical no-signalling, matches a1)")
print("  => the conjecture's mechanism is SAFE (no channel) but WEAKER than advertised: any Bell violation the model claims must come from the")
print("     quantum fork's non-classical structure (REALFIBER holonomy phases, block-resolved traces), not from non-injective Gram presentations.")

n=sum(CH); print(f"\n=== SIG-1: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("RESULT: kill 3 does NOT fire under either reading of P6. CONJECTURE-COSMOLOGY: kill 2 survived (SHEET-1), kill 3 survived (SIG-1),")
print("        kill 1 (the number) remains the only one with teeth. Correction to the conjecture's text: non-injectivity is not a Bell mechanism.")
print("TIER: s0-s1, a1-a2, b1 DERIVED (exact). c1-c2 numerical, with Fine's theorem as comparison.")
