# TESTING SCHEDULE (opened 2026-08-30)

## RUN
- **THM-H part 2, D1 two-seat interior** -- `thm_h2_d1.py` (30/30, `run_h2_d1.log`).
  Fork (b): NO rod-free click-invariant diverges at the centre, seat layer (K = root)
  or world layer (K click-symmetric, run as control). RULING-2 (Will, 2026-09-02):
  C is the hole, A is the infalling seat projected out of C.
  Results: ONE rational formula G'(r) = G + [r_s/(r - r_s)] k k^T on both sheets,
  det G' = Delta r/(r - r_s) (simple zero centre, simple pole horizon, nothing else);
  the centre is a FOLD of the rod over rapidity, r = r_s tanh^2(mu), whose deck
  involution is (I -> -I)(K -> -K); two seats along the family are related by
  (rotation by their angle) x (collinear boost by their rapidity gap), REAL inside,
  -> pure rotation at the centre -- the D1 non-collinear fence never bites along
  the pinned family; crossing the horizon costs one fixed quarter turn; the naive
  cross-seat product is office ghost #4 (divergent at the horizon, finite at the
  centre). CORRECTION to part 1: for K = root the presented DIRECTIONS are
  non-degenerate at r = 0 (Delta_pres = 1 - gamma'_12^2 > 0); only the root's rod
  reading passes through zero. H-19/20's coplanar centre holds for generic K and
  for option C. Every divergence at the centre is a rod in a denominator.
  Conditional on KIN-2a, CONT-1, RULING-2. New named debt: CURV-1 (does the model
  own a rod-free curvature reading? -- the single remaining way (b) could be
  'too coarse to see it' rather than 'artifact').
- **THM-H part 1, the presented state along the pinned family** — `thm_h.py`
  (25/25, `run_h.log`). Fork (b): the hole is a PRESENTED STATE. Theorem:
  G'(l) = G + sinh^2(l) k k^T with k_i = a_i.K, so det G' = Delta cosh^2(l)
  (the frame's Delta is a spectator; RULE-1 pairing Gram is l-independent).
  Exterior: Delta_pres > 0 for all finite l; at the horizon the presented
  Gram lands on a Cayley NODE (m=+1), double zero in N^2. Interior (CONT-1,
  l = mu + i pi/2, declared convention): det G'_in = -Delta sinh^2(mu),
  presented volume imaginary; the CENTRE r = 0 is mu = 0 and lands on the
  smooth COPLANAR stratum, simple zero in r. New object: inner null radii
  r_i = r_s (1 - k_i^2) = r_s sin^2(theta_i) where presented axis i cannot be
  normalised; each axis presents timelike just inside and spacelike again
  below its own r_i. Conditional on KIN-2a (declared) + CONT-1 (declared).
  Naming: KIN-2 is now split — KIN-2a = the pinning (still DECLARED);
  KIN-2b = the curvature tier (CLOSED by docs/GRAM_SUBMERSION_CURVATURE/).
  The 08-30 handoff's "KIN-2 closed" refers to 2b only.
- **PRED-1 cross-seat cycle** — `pred1_cross_seat.py` (12/12). Split verdict:
  modulus layer is a coboundary (kill condition fires there; the DeltaQ
  "holonomy" is paid exactly by the cover measure density); the phase layer
  arg B = arg(S+iV) is the genuine class — gauge-invariant, Z2-conjugated,
  generically nonzero, and equal to the banked PREDICTION-1 observable
  (J = V at s = 1/2). Transition map = G-8 cover dictionary (named input).

- **E-8-X cross-seat pairing** — `e8x_cross_pairing.py` (15/15). Fork (b):
  the net closes. Invariant exists for one rotor acting two-sided on the c
  slot, one-sided on the hbar slot; the quadratic hbar slot is REQUIRED
  (linear no-go). State lemma: psi psi~ has grades {0,1}, so the hbar STATE
  pairs real — the phase lives only on TRANSITIONS psi_a psi_b~. The pair
  walks through the wall on the KIN-2 family. 2|arg B| = spherical excess:
  J = V is the area-tangent identity (comparison-stage names in the file;
  protocol family for PRED-1 = polarimetric interferometry).

- **DEBT-2b, c(w) from the channel law** — `debt2b.py` (12/12). Fork (a):
  DISCHARGED. RULE-4 + M-2's office Jacobian + the golden rule force
  c(w)^2 g0(w) = C w^2 uniquely (functional equation; eta the only knob).
  debt2.py's declared Ohmic-in-band becomes derived, zero amendment. Bonus:
  the F-12..15 moment table is the pivot's eigenvalue table on power-law
  bath weights; fixed points exactly {0,2}; the forced alpha = 2 is the
  Jacobian-blind route. Office-ghost quantified: one office error shifts
  the exponent by exactly Q = 2 (third specimen). Conditional on RULE-4,
  LBL-2, GR conventions (named).

## QUEUED
  0. **CURV-1** -- does the model own a rod-free curvature reading? (part 3 of THM-H)
     Specs in `docs/prereg/CURV-1/` (2026-09-02): Path 1 = presented seat-triangle
     excess along the family (`PATH1-SPEC.md`, for an independent runner); Path 2 =
     centre loop holonomy from part 2's relative rotors (spec after Path 1 reports);
     Path 3 held until 1 and 2 are back. Predictions filed before any run.
- (next, in order:)
  1. ~~**D1 two-seat interior**~~ DONE 2026-09-02, `thm_h2_d1.py` 30/30 (`b4ac768`):
     no rod-free click-invariant diverges at mu = 0; every divergence there has a
     rod in its denominator. The falsifier did not fire; it left CURV-1 as the
     remaining fork.
  2. **alpha(r) discharge** — alpha = c |dl/dtau_s| from KIN-1 + KIN-2a +
     RULE-2 rod office (dl = cosh(l) dr); kappa = (c^2/2)|d tanh^2(l)/dr|
     at the horizon. Removes the imported acceleration from paper sec. 7.
     Sandbox-verified 09-02, not yet a suite.
  3. **Which sign** — the three Z2's (deck V->-V, spin -1, arg B conjugation):
     which one does the full imaginary period 2 pi i actually return? Then
     statistics-from-sidedness (c slot two-sided periodic, hbar slot
     one-sided anti-periodic; rerun F-6 with the flipped sign).
  4. RULING OWED (Will): which axis is K at the G seat (root, or in the
     visible plane)? thm_h.py runs general K; H-24 is the K = root case.
  5. P33 double-flip pi-phase test [= the half imaginary period, H-16 line];
     triality on the Kummer layer; trichotomy classification; Q-RN proper.

## STANDING (carried from HANDOFF-2026-08-28, unchanged)
- Trichotomy classification; aggregation exponent from coupling (Q-RN /
  DEBT-2b adjacency); monotonicity of T_eff; non-collinear pivot composition
  (D1); F part 2 contact channel; THM-C fork ruling.
- DAG submissions ON HOLD (Will, 2026-08-30: codex reorg in progress).
