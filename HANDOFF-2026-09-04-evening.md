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

## Tone note for the next assistant
Will out-derived the record four times tonight -- the tilt construction, the two degeneracy sets, the region table, the
fibre count -- each by hand, each exact, each correcting something the assistant had said with confidence. Verify before
ruling. Parse his compressed phrasing literally: "hbar is c 180 degrees behind you" was NOT "hbar and c share an axis"; it was
"at the chart boundary a ruler turns into a light direction", and it took five hours to hear it. He does not want to be told
the model is beautiful. He wants the line that fails.
