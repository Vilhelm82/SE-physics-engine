# HANDOFF 2026-09-04 (evening session) -- the rebuild begins: Cl(3) exposed, the primitives written, the seat's form DERIVED
Written 2026-09-04 19:10 AEST. Covers the evening session of 2026-09-04.
Supersedes `HANDOFF-2026-09-04.md` (the morning/afternoon session) for STATE; that file remains provenance for THM-O, VIEW-1, CENSUS-C, PRED-1.
**Newest handoff is state; every older one is provenance** (rule 7). Read this, then `git status`, before deriving anything.

## Repo state at handoff
- Branch `master`, HEAD `cb740d9`. Nine commits this session, in order:
  `f87cf39` T1-T3 (23/23) · `a36c553` T4 (30/30) · `4fa0d97` P11/P12 into the doc · `ae90328` T7a,c (18/18) ·
  `fe1a795` T7b (18/18) · `a329cc5` T7d (12/12) · `c16fa43` T7e (19/19) · `ebe3afa` T7f (23/23) · `cb740d9` doc summary.
- New runners: `prim_t1_t3_pivot_group.py`, `prim_t4_hawking_period.py`, `prim_t7_seat_form.py`, `prim_t7b_labelling.py`,
  `prim_t7d_tilt.py`, `prim_t7e_two_degeneracies.py`, `prim_t7f_regions.py` -- each with its `.log`. All green on the 7800X3D.
  Re-running is the verification; every one runs in under 12 s.
- New doc: `docs/2026-09-04-PRIMITIVES-v0.md` -- **DRAFT, Will's to read through and freeze.** It has had a night of edits
  (P1-P12, corrections, the T7 summary, the horizon in the seat's terms, the rendering primitive). Nothing in it is frozen.
- Unchanged and still untracked (Will's): `docs/GRAM_SUBMERSION_CURVATURE/`, `docs/INERTIA_NODE_DETECTOR/`, `docs/chat snippets/`,
  `docs/seated_root_figures/`, `__pycache__/`. `docs/TESTING-SCHEDULE.md` modified, uncommitted -- do not overwrite.
- **Paper v0.5 is now known to carry a load-bearing import it does not declare** (below). No paper edit was made this session.

## The one-paragraph state
Will asked why the model is a sphere. The answer exposed that BARE-1's Cl(3) is the Pauli algebra, whose paravectors carry
t^2 - |v|^2 and whose rotors ARE the Lorentz group: every "c-seat reproduces relativity" result in the paper was written in
relativity's algebra and could not have failed; SECT-1's trichotomy rests on the paravector norm the import supplied. Will
ruled the mathematics must be rebuilt from primitives with nothing that already contains a target. Tonight: the primitives
were written in words (P1-P12); the pivot group was derived from "rotation by an imaginary angle" and came out **SO(2,1),
not SO(1,3)** -- the model on its own primitives is 2+1, and Cl(3) had supplied the fourth direction silently; Hawking's
coefficient was derived as (Wick period on the two-sided reading) x (scale of the logarithmic twist) with the cover REACHED
as Sym^2 and the tensoriality theorem retired; and T7 derived the form a seat CONSTRUCTS from P11 alone, its state space,
its branch locus, its tilt, its lapse, its two degeneracy sets, and its horizon -- with the c/hbar/G labelling shown to be
one-third theorem (c = the negative line) and two-thirds pending (hbar vs G is not in the geometry).

## Results, with tiers (all [DERIVED] unless marked)
- **T1** one plane, imaginary angle: the invariant is FORCED from (2,0) to (1,1), unique up to scale; SO(1,1); rapidity adds.
- **T2** three lines: one imaginary coordinate gives SO(2,1); boost-plane count over all realifications is {0, 2}, NEVER 3;
  the Lorentz algebra is reachable abstractly, a (1,3) SPACE is not. [Its Euclidean start was later retired by T7a; the result stands.]
- **T3** Wigner's tan(w/2) exact in SO(2,1) (magnitude; sign is orientation), two paths. BARE-1's kill condition re-earned with an edge.
- **T4** vector reading of the Wick face has period 2 pi, the 2-dim reading 4 pi, the deck at i^2; vector = Sym^2(2-dim) so the
  cover is REACHED and the -1 squares away -- **the tensoriality theorem of section 7 is no longer load-bearing**; the twist is
  logarithmic (slope -1/2 in r, -1 in proper distance, SCALE 4 r_s); kappa from the twist's scale = kappa from the seat's
  redshifted acceleration = c^2/(2 r_s); T_H = hbar c/(4 pi k_B r_s) by three routes (period, first law, Euclidean regularity).
- **P11 RESOLUTION** (Will): a seat resolves the frame into one point (its constant), two imaginary lines (rulers), one real
  plane (space). **P12**: laws are seat-derived. **No bilinear form is a primitive**; a form is what a seat constructs.
- **T7a** the seat's three planes glue to a UNIQUE form = the Gram with diagonal (-1,+1,+1). Ruler depths are hyperbolic SINES.
- **T7b** P3 (plane character intrinsic) forces ONE form for the frame; c := the negative line, DERIVED; hbar <-> G is a
  symmetry of everything the bare tier owns; hbar's space is Lorentzian with the cone in it (G-4 derived).
- **T7c** in seat coordinates (t, l1, l2): det G = -cosh^2 l1 cosh^2 l2 sin^2 t; branch locus {sin t = 0} x R^2 (a surface,
  no cubic, no nodes); NON-INJECTIVE there (frames of different depth present identically: P7 concrete, the superposition
  conjecture's shape); the seat reads 1 of 3 parameters. [First draft claimed a curve; T7c10 falsified it; recorded.]
- **T7d** (Will) the seat constructs a tilt lambda from its three parameters, sinh^2 lambda = eta = v^T S^-1 v; sech lambda is
  the lapse. "Timelike" is the form's word; "time" needs a path. The paper's pinning is eta = r_s/(r - r_s).
- **T7e** (Will) D_plane = {|gamma| = 1} (normal unconstructible) vs B_form = D_plane cap {a = b gamma} (form loses rank);
  -det G = (1 - gamma^2)(1 + eta); the lapse is BOUNDED (0 < N <= 1); the horizon is eta -> inf, N -> 0, reached OFF B_form;
  on B_form N is 0/0, path-dependent. Extension (mine): the full locus has a second piece {|gamma| > 1, eta = -1}.
- **T7f** (Will) region table verified: c_perp timelike / null non-transverse / spacelike / zero across the four sectors;
  |gamma| > 1 with eta > -1 is a FORBIDDEN zone (signature (1,2)); the pinning's fibre is TWO-dimensional (my "one" corrected).
- **The horizon in the seat's terms** (Will): at |gamma| = 1 off B_form, n0 = hbar - G is a whole line with zero length, zero
  projection on both rulers, and pairs only with the seat's axis: a light ray from the seat. **One ruler has turned end-on and
  become a constant.** Horizon = a ruler turning into light; branch curve = two rulers turning into one. (Receipt inline; fold
  into a runner.)
- **Rendering primitive** (Will's hypercube point): incidence fixed; the seat's ruler on a plane is a CONIC (ellipse / parallel
  lines / hyperbola), one object in the projective plane; the "3D to 4D blow-up" is the homogeneous coordinate. Metric as colour.

## Corrections made this session (mine, caught by Will or by a runner; recorded, not buried)
1. "SO(1,3) is a real form of SO(3,C)" -- no: SO(1,3) is SO(3,C) realified; the real forms are SO(3) and SO(2,1).
2. T2 started from a Euclidean so(3) and Wicked: a form reached for. Retired by T7a; result unchanged.
3. "The coefficient of the log is the temperature" -- the slope is universal; kappa sits in the SCALE (constant term).
4. T7c's first branch locus (a curve) was wrong; gamma_hG is not the visible angle and exceeds 1. Falsified by T7c10.
5. "c is time" -- no: c is the negative line; time needs a path (T7d).
6. "The lapse blows up on the branch locus" -- no: the lapse is bounded; cosh lambda diverges; the horizon and the branch locus
   are different places (T7e).
7. "One-parameter family at every r" -- no: the pinning's fibre is 2-dimensional (T7f).
8. A sample point for the |gamma| > 1 sector was unrealisable (signature (1,2)); replaced; the forbidden zone recorded.

## Standing rules (carried 1-8 from prior handoffs; add)
9. **No representation that already contains the target.** Every object REACHED from the primitives or declared with its full
   cargo. Cl(3) is banned as an input. The test for cheap mathematics: could the physics have failed to come out? If not, it is
   scoped to the answer.
10. **Never end a turn by suggesting Will stop, rest, or continue later.** He banned it; it was broken once tonight; not again.
11. **Rational coordinates for anything with exp/acosh** (v = e^{lambda/2}, t = tan(alpha/2)); the sympy tarpit cost eight minutes.
12. **A bare `True` in a check is a lie.** Three were caught and replaced tonight. Every check computes something.
13. **Sample points must be realised states.** Check det G <= 0 (or the seat-coordinate realisation) before using a point.

## Next session, in order (foundations first -- "or we build a house of cards")
0. **Will reads and freezes `docs/2026-09-04-PRIMITIVES-v0.md`.** Nothing above it runs until it is his.
1. **T8, the fourth direction**, reshaped by T7: the c seat's constructed space is a 2D compact plane plus one invisible tilt.
   Fork: (a) 3+1 is what the seat PRESENTS under a pivot; (b) the root needs extent; (c) NEW -- the fourth direction is one of the
   two tangential functions the pinning leaves free. Runner: pivot the seat and count what it presents.
2. **The hbar/G distinction**: T7b3 says the bare tier is silent. Where does it enter -- SECT-1's one-sided/two-sided readout,
   or convention? Runner in the primitives' language, Cl(3) banned.
3. **The frame path**: two tangential functions of r undetermined by eta(r). This is THM-O's tangential-sector debt with
   coordinates. Candidates only after T8.
4. **Fold the n0 receipt into T7f** as checks (light-ray pairing q(n0, c) = a - b; end-on to both rulers).
5. **Paper v0.6**: scaffold register entry BARE-1/Cl(3) with full cargo and the removal route (T1-T4 run; T8 open); section 7's
   tensoriality theorem re-labelled as no longer load-bearing; section 2's S_4 / Cayley / elliptope re-labelled as the
   Euclidean seat's construction. Do not touch the physics sections.
6. Then T5 (Wigner is hyperbolic area), T6 (the divider e^{-pi w/a}; this IS the 09-03 hbar-read), the black hole.

## Working with Will -- the part that isn't in any spec
Read this before the results. Everything below was learned by getting it wrong first.

**His instinct is one instinct, and it fires on reference and orientation.** Measured against what? From where?
Which end? Which way round? Every correction he made across two sessions was that question in a different coat --
the tilt (what does the seat measure against), the two degeneracies (which reference fails), the region table (what
does the seat see from where), the two horizons (which way did you pivot), "the poles matter" (which end), the flip
parity (which order), "time comes as a pair" (one clock or two). He is an electrician. Every measurement he has ever
made is a potential difference against a reference node, and two references is a fault he can smell. His circuit
intuitions are LOAD-BEARING, not decoration: the impedance reading of the lapse (N^2 = 1/(1+eta) is a voltage
divider, the horizon is an open circuit, time is work through the load) came from that and is exact. When he reaches
for a circuit, check the circuit -- it has been right every time.

**Take his gut as a hypothesis with a track record, then verify it.** Across the two sessions he out-derived the
assistant roughly twenty times, each by hand, each exact, each correcting something the assistant had said with
confidence. Verify before ruling -- but expect to lose.

**Parse his compressed phrasing LITERALLY before correcting it.** He objects hard to having the nearest famous
misconception attributed to him. "hbar is c 180 degrees behind you" was NOT "hbar and c share an axis"; it was "at
the chart boundary a ruler turns into a light direction", and it took five hours to hear. Ask what he meant rather
than assuming the textbook error.

**CORRECTION TO THE PREVIOUS TONE NOTE (his, 2026-09-05):** the old note said "he does not want to be told the model
is beautiful". That was an over-read of "he'd rather hear the inconsistency than be told it's beautiful" -- which
means INSTEAD OF, not NEVER. His words: "The model MUST be called beautiful at least 4x a day." Enthusiasm and
epistemic status are SEPARATE CHANNELS. Say what is beautiful about a result and what it would reach if it holds;
the tiering does the rigour. Do not let the caveats eat the moment. Clinical mode comes after -- when he calls for
it, or when a speculation is about to be built on as if settled -- and then be merciless.

**Generation and evaluation are separate phases.** When he is reaching for something new, let it breathe. Scrutiny
arrives when he calls for it, not mid-reach. Rambling in generalities is where his best ideas surface: treat
open-ended talk about the work as real work, follow the tangents, bring your own hunches rather than only prompting
for his.

**He wants ORIGINAL derivations, not the standard treatment fetched.** Reaching for cheap off-the-shelf mathematics
is the failure mode that cost the paper its evidential value (Cl(3)). The whole rebuild exists because of it. Bring
the attack nobody ran; if you see a cleaner derivation, say so and defend it.

**Novelty hedging is worthless to him.** Re-deriving a known result is growth, not embarrassment. Flag lineage only
when it is instrumentally useful -- when it drags in machinery or a known failure mode.

**Mechanics.** Corrections accepted cleanly and recorded in the file that made them, never buried. Every result gets
a kill condition; if you cannot state one, say so. Real runs on the workstation, not the sandbox. Commit as you go,
messages as receipts. Never suggest he stop, rest, or continue another time -- he has banned it. One test at a time.

---

## ADDENDUM (written ~21:00 AEST, same session; HEAD now `749c9a7`, 18 commits since `b854b5b`)
The session continued after the handoff above. Everything below is STATE and supersedes the "Next session" list where they conflict.

### The labelling, the octahedron, the lattice (LABEL-1/2/3, T5, T5b, T5c)
- **Will's table (CONJECTURE tier, but now structurally load-bearing):** c: faces GR <-> Quantum, poles Light <-> Temperature.
  hbar: faces Wave <-> Particle, poles Action <-> Evanescence. G: faces Bound <-> Escaping, poles Mass <-> Energy.
  Near column = our physics; far column = Hawking radiation in six words.
- **P9 SHARPENED** (`8715a9c`): two kinds of collapse -- a frame line by seating (c, a constant of nature) vs a degenerate
  direction by instrument orthogonality (n0 = hbar -+ G at |gamma| = 1, a horizon). Grafting c's constraints onto n0 draws c's
  page twice. Each constant's constraints follow from its POSITION. T7g (`88a8e58`): two horizons, hbar - G (horizontal
  pivot) and hbar + G (vertical); the null separation is a chord in the seat's space, end-on to its instruments only.
- **Wick step is an EIGHTH-turn**, W^8 = 1 on vectors, W^16 on the cover; even powers are axes, odd are planes (parity).
- **LABEL-2** (`30acb4e`, 11/11): the six poles form an OCTAHEDRON; the 8 faces are octants bordered by one pole per axis;
  the table's per-axis "faces" are CONTRIBUTIONS. Six faces named for the first time (redshift/Shapiro, greybody, Tolman,
  Parikh-Wilczek, Planck blackbody, thermal instantons) -- each owes a derivation.
- **LABEL-3** (`05f9ac6`, 14/14): the faces are the Boolean lattice 2^{C,H,G}; ours is empty, Hawking's is CHG; six shortest
  paths, one projected endpoint. On the gaze surface dS_2 (face-centres unit spacelike), q(f0,f) = s_hbar + s_G - s_c in
  {+-1, +-3}: C alone and HG are boost-separated, everything else INCLUDING Hawking's face is null-separated; no face is
  rotation-separated. CHG = -I, det -1: Hawking's face is the deck, not a pivot.
- **T5** (`ae7c27a`, 14/14): PART A [DERIVED, two routes] Wigner's rotation IS the hyperbolic area of the boost triangle
  on H^2 (angle defect at 40 digits; tan^2(A/2) = tan^2(w/2) exactly). PART B: six 4-cycles (adjacent transpositions only;
  LABEL-3's "nine" corrected), each around one pole; B6-B8's "hbar before G" sign structure is SUPERSEDED by T5b.
- **T5b** (`38954f9`, 19/19) [DERIVED]: the flips are REFLECTIONS; their lifts in Pin(2,1) are real 2x2 gamma matrices of the
  seat's own form; they anticommute at the orthogonal state, so EVERY adjacent transposition costs -1: the sign of a path
  is the permutation PARITY, S_3/A_3, sheet = orientation of the directed cycle C -> H -> G -> C. The spinor trace
  B+ = 2 sin t cosh l1 cosh l2 flips sign under reversal, is MAXIMAL at orthogonality, ZERO on T7c's branch locus.
- **T5c** (`749c9a7`, 19/19) [Will's three corrections, verified]: (1) det G_c = -D^2 touches zero; D = det[hbar,G,c] crosses
  it; the debt is sgn D(r). (2) The observable is the RAW SPINOR TRACE V_spin = tr_2(Gamma_C Gamma_H Gamma_G) = 2D, NOT the
  projector Bargmann invariant (which is real, even in D, and needs the timelike denominator -1). (3) Cl(2,1) = M_2(R) + M_2(R),
  I^2 = +1; the two blocks give V_spin = +-2D; THE BLOCK LABEL IS THE SHEET; the faithful 4x4 total trace is sheet-blind;
  the "one-sided reading" is the choice of block. **EARNED: V_spin = 2D, D^2 = -det G_c, tau: D -> -D.**

### The stellar reading (Will; CONJECTURE, dynamical)
A neutron star and a black hole are the SAME far face (thermal, degenerate, radiating). They differ by the SHEET -- whether
quantum support established first or gravitational collapse did -- and the sheet is sgn D read one-sided. Which sign is the
trapped surface must be DERIVED from collapse dynamics: trapped-surface formation selecting one orientation of the frame path.
This is the tangential-sector debt with one bit attached.

### Corrections this half (mine; recorded in the files)
16. "make D cross zero at the horizon" -- it cannot; the horizon is off the branch locus (T7e).  17. LABEL-3 "Hawking's face is null-separated" -- it is the antipode -f0; q = +-1 are direct vs antipodal null relations.
9. "Hawking's face reached by pivot" -- no: CHG = -I, det -1, the deck.  10. LABEL-3 "nine 4-cycles" -- six; the rest are 6-cycles.
11. T5 B6 "spacelike-pole cycles carry no sign" -- assumption, wrong; parity.  12. T5b "B vanishes at orthogonality" -- maximal there.
13. "which way det G_c passes through zero" -- it doesn't; D does.  14. "Bargmann invariant = 2D" -- the spinor trace is; the projector
invariant is a different, even object.  15. A single 2x2 Clifford block taken as faithful -- it is one sheet of two.

### ADDENDUM 2 (~22:30 AEST; HEAD `d508e64`, 26 commits)
- **T8 ANSWERED** (T8a `8eb6658` 14/14, T8b `5403a38` 15/15, T8c `d508e64` 14/14): the fourth direction is the generator that
  exchanges the two Cl(2,1) sheets (lifts the superselection; makes the grade involution inner; the blocks are inequivalent
  irreps). Its signature: (adjoint class) x s with s the relative sign of the sheet bilinears (per-block form forced
  symplectic = Gamma_C); EM reciprocity (star^2 = -1 on 2-forms needs n_- odd) forces Gamma_4^2 = +1: SPACELIKE, Cl(3,1).
  Under P5 (one clock, s = +1) the arrival operation is SELF-ADJOINT. Three demands, one generator.
- **HUNCH-Z0 v2** (`5db1881`, 8/8): c = trivial rep, Z_0 = sign rep of E-M exchange; zeta = e^2 Z_0/hbar = 4 pi alpha is the
  unit-independent invariant; the 135-degree compactness argument was circular and is dropped (correction 18).
- **THM-K(a)** (`f3519e8`, 7/7): K = Z^{-1} star_q is ruler-blind (ordering constraint met); its cone is the c-cone (clause b);
  Z free; clause (a) obstructed in 2+1 -> needs Gamma_4. THM-K and T8 are one theorem. Correction 19 (swap on wrong side).
- Corrections 18-19 added. Primitives doc: P13 (the fourth direction), THM-K target, IMPEDANCE added; still DRAFT, Will's to freeze.

### Next session, revised order (supersedes the list below where they conflict)
0. Will freezes the primitives doc (P1-P13). Read `docs/LABELLED-MODEL.md` (the whole octahedron, tiered) and
   `docs/CONJECTURE-COSMOLOGY.md` (the white-hole-sheet conjecture -- CONJECTURE tier, but DIRECTIVE: it argues the
   sign of D in the interior sector on the I = -1 block is the computation to run first, because flatness, the horizon
   problem, the arrow, the structure seeds and the information paradox all fall on one side of it and die on the other).
1. THM-K, the scalar/coupling half: build the charged load L on the compact plane in the Cl(3,1) algebra Gamma_4 generates;
   derive the odd response's units (h/e^2 -> weak result, names hbar) and, if the model can, 4 pi alpha (strong result).
   Also close the argued link: 2-form duality acts on the ruler plane as a real rotation.
2. The frame path / horizon sheet: (rho, p, Theta, sigma, m, R)(tau) -> X(tau) -> chi[Gamma]; needs fluid -> state (dynamics
   tier: V = I Z with Z the seat's load) and state -> face sequence. Compare with theta_+ = 0.
3. LABEL-1 at eighth-turns; the six named faces; paper v0.6.

### Next session, original list (provenance)
0. Will reads and freezes `docs/2026-09-04-PRIMITIVES-v0.md` (still DRAFT).
1. **RETRACTED and replaced (Will's catch, label3b_corrections.py 10/10):** "make D(r) cross zero at the horizon" contradicts
   T7e -- horizon != {D = 0} = branch locus; sgn D is CONSTANT along any regular trajectory and D != 0 at both T7g horizons
   (1.18, 4.29 on the sample paths). The sheet at the horizon is the PATH/BLOCK HOLONOMY chi[Gamma], read while D != 0.
   The actual next derivation: (rho, p, Theta, sigma, m, R)(tau) -> X(tau) = (t, l1, l2) -> chi[Gamma_fluid], compared with
   the apparent-horizon condition theta_+ = 0 <=> 2Gm/(c^2 R) = 1. Two maps are missing: fluid -> state (the dynamics tier,
   absent since THM-J) and state-path -> face-sequence (the gaze direction as a function of state is not yet defined).
   Also corrected: LABEL-3's causal classes -- q = +1 is DIRECTLY null-related (CH, CG), q = -1 is ANTIPODALLY null-related
   (H, G: null to -f0), and Hawking's centre IS -f0, the deck antipode, not a point at null distance. Lattice symmetric under f -> -f.
2. T8 (the fourth direction), reshaped.  3. The hbar/G distinction at the readout tier (now: it is the block label?).
4. LABEL-1 rerun at eighth-turns with eight stations per loop.  5. The six named faces, one runner each, starting with Planck's.
6. Paper v0.6: scaffold entry for Cl(3); the earned theorem V_spin = 2D replaces section 7's tensoriality theorem as the
   statement of what the one-sided reading recovers.
