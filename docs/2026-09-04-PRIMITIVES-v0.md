# PRIMITIVES v0 — the seated root without imports
Date: 2026-09-04. Status: DRAFT, uncommitted, Will's to correct before freeze.
Purpose: the ontology in words, with no algebra that already contains a target.
Rule: every mathematical object downstream must be REACHED from these, never REACHED FOR.

## What triggered this
The paper's pivot (BARE-1) was declared as a hyperbolic rotor in Cl(3). Cl(3) is the
Pauli algebra; its paravectors carry t^2 - |v|^2 structurally; its rotors ARE the Lorentz
group. Every "c-seat reproduces relativity" result was written in relativity's algebra
and could not have failed. SECT-1's trichotomy rests on the paravector norm the import
supplied. The three-axis geometry did not produce Lorentz. Lorentz was the ink.

## Primitives (words only)
P1. A point (the root).
P2. Three lines through the root. Each line has two ends (poles).
P3. Three planes, each spanned by two of the lines. Each plane has a fixed intrinsic
    metric that never deforms. The planes need not be alike.
P4. The state is the three angles between the lines. Nothing else.
P5. A seat is one line taken as reference. The seat cannot see its own line; both
    poles collapse to the root.
P6. A pivot is a rotation of the frame by an IMAGINARY angle. It is the operation that
    turns an angle into a rapidity. It is not a spatial motion.
P7. A presentation is what a seat reads after a pivot. Presentations are non-injective:
    distinct states can present identically to one seat.
P8. Axes cannot be manipulated in isolation. Moving one line moves the frame or shows
    the compression/stretching of the plane between converging lines.
P9. A constant is what a view collapses to a point. Constancy is a property of the
    view, not of the quantity.
P10. [Will's ruling, 2026-09-04 late; WORDING CORRECTED by Will same night] c is the
    NEGATIVE line of the frame's form (T7b2: derived, not declared). "Timelike" is the
    form's word for it; "time" is NOT. A (2,1) form contains a whole negative cone and
    identifies no vector with time, supplies no arrow, produces no evolution. What the
    seat constructs is a TILT: decompose c against the rulers' span P_c = span{hbar,G}
    (Gram S, |gamma| < 1), p = (hbar,G) S^-1 v, c_perp = c - p, eta = v^T S^-1 v;
    then q(p,p) = sinh^2 lambda, -q(c, n) = cosh lambda with n the unit normal to P_c,
    and sech(lambda) is the lapse -- the time-dilation factor between the seat's clock
    and the resolved one. Time is an ordered duration along a path through states;
    without dynamics T7 gives causal geometry and no time. Calling c "time" erases the
    tilt that produces dilation. (prim_t7d_tilt.py, 12/12.)
P11. [Will, 2026-09-04, RESOLUTION] A seat resolves the frame into one point, two
    lines, one plane. The point is the seated axis (its constant, both poles
    collapsed). The two lines are the other two axes, which are the edge-on traces of
    the two planes through the seated axis; the seat uses them as RULERS (unit
    references). They are imaginary lines: the edge of a plane you can enter only by
    an imaginary angle. The plane is the one not through the seated axis, seen
    face-on: unbounded, real -- what the seat calls space.
P12. [Will, 2026-09-04] Laws are seat-derived. The constants and laws a seat reads are
    derivable from that seat's resolution and no other's. Another seat derives
    another physics from the same frame.

CORRECTIONS forced by P11/P12:
 - P3 keeps "planes have intrinsic character (hyperbolic / compact), never deforming"
   and LOSES the word "metric": metric is a seat's word for its rulers (P11).
 - P4 is stronger than written: the three angles are ARBITRARY and unconstrained by
   the model. The elliptope (Delta >= 0, |gamma| <= 1) was the Euclidean seat's
   construction, not the model's.
 - There is NO bilinear form among the primitives. A form is what a seat's resolution
   constructs. T2's SO(2,1) is the c seat's construction (compact plane face-on, two
   hyperbolic edges). From hbar the same frame resolves with a hyperbolic plane
   face-on: hbar's space has the cone IN it. Same group, different resolution.
 - The Cl(3) crack and the "Euclidean then Wick" crack in T2 are the same crack: a
   form reached for. P11 says where forms come from.

## Geometry the rulings imply (T2 + P10), in words
The timelike line c pierces a two-sheeted surface: +c (light) on one sheet, -c
(temperature) on the other. The spacelike lines hbar and G pierce a one-sheeted surface
whose waist is a circle in the compact plane, with four marked points +-hbar, +-G.
Between them is the cone. Real boosts run down a meridian and never reach the cone
(the twist is logarithmic). The only crossing is the imaginary quarter-turn. Four
stations per meridian: c, iK, -c, -iK. Three closed journeys: pole-to-pole through the
cone twice (the thermal circle; holonomy 2pi on vectors, -1 on the cover); around the
waist (the seat's axial spin); down-along-up with a corner (the seat cycle; holonomy =
enclosed area). The black hole's horizon is the cone, its interior is the waist
(CONT-1's quarter turn), its centre is -c (CURV-1's spin -1). Matter at a Wick station
either crosses (continues the loop) or turns onto the waist (squeezes); the loop is a
divider; what completes is Hawking. [conjecture tier on top of proved stations]

## Banned as inputs (may appear only if DERIVED from P1-P9)
- Clifford algebras, Pauli matrices, spinors, any double cover
- Minkowski signature, the light cone, rapidity as a given
- Hilbert space, the Born rule, hbar-dependent angles or commutators
- The round metric on S^2, solid angles, the view sphere
- Positive-definiteness of the angle matrix (the elliptope), the Cayley cubic
- Any named physical theory as an arbiter

## First derivation targets, in order
T1. (2,2) rung: one plane, two lines, one angle. From P6 alone, what group do the
    pivots form, and what happens to the plane's quadratic invariant? Prediction, held
    out: the invariant changes signature (2,0) -> (1,1); the group is SO(1,1); rapidity
    adds. If so, 1+1 relativity is DERIVED from "imaginary angle", not imported.
T2. (3,3) rung: three lines, three planes. Which planes receive the imaginary angle?
    All three -> SO(3,C) as a real group = SO(1,3), Lorentz, with a fourth (scalar)
    direction forced. One plane -> SO(2,1), one line timelike. This is where P3
    ("planes need not be alike") becomes a choice the model must make or derive.
T3. Wigner rotation from T2, not from Cl(3). If two non-collinear pivots compose to a
    pivot times a real rotation, BARE-1's kill condition is re-earned with an edge.
    STATUS: T1-T3 RUN, 23/23, commit f87cf39. T1: (1,1) forced, SO(1,1), unique
    invariant form. T2: SO(2,1); boost-plane count in {0,2}, never 3; Lorentz ALGEBRA
    reachable, (1,3) SPACE not -- the model is 2+1 on its own primitives. T3: Wigner
    magnitude exact in SO(2,1), sign is orientation.
T4. Hawking's coefficient from SO(2,1): period of the Wick face on the vector reading
    (2pi) vs the 2-dim reading (4pi); the cover REACHED as Sym^2; the pinning's twist
    logarithmic with universal slope; kappa from the twist's scale = kappa from the
    seat's redshifted acceleration; T = (2pi) x kappa; first law and Euclidean
    regularity agree.  STATUS: RUN, 30/30. The tensoriality theorem is no longer
    load-bearing: the two-sided reading has no -1 to hide.
T5. Wigner is area: the T3 rotation equals the hyperbolic area of the triangle whose
    sides are the two boosts (Gauss-Bonnet on H^2 = SO(2,1)/SO(2), curvature -1).
    If so, the seat-cycle holonomy of section 9 is derived on the model's own surface
    and the Bargmann phase on S^2 was its spherical shadow. Kill: any discrepancy.
T6. The divider: the amplitude for a mode of winding w on the waist to cross the cone
    under a boost of rate a. Claim: e^{-pi w / a} per crossing, nothing put in;
    two crossings = Boltzmann. Kill: not exponential in w/a -> the divider picture
    is dead. (This IS the 09-03 hbar-read protocol; comparator tanh r = e^{-pi w/a}.)
T7. [REFRAMED by P11] The state space is three free angles; there is no constraint.
    What T7 now asks: (a) DERIVE the form a seat constructs from P11 -- one point,
    two imaginary rulers, one real plane -- and show it is (2,1) for the c seat
    without starting from a Euclidean R^3 and Wicking. (b) Do the same for the hbar
    seat and the G seat: what signature does each construct, which plane is its space,
    where is the cone in its view. (c) What must hold among the three arbitrary
    angles for a given seat's resolution to be well-defined (the seat's own branch
    locus, replacing the elliptope). Only then does SO(2,1) stand as REACHED.
    This is FOUNDATIONAL and runs BEFORE T5, T6, T8.
T8. The fourth direction. T2 says 3+1 is not on three lines. Fork: (a) 3+1 is what
    the c seat PRESENTS of a 2+1 frame -- the extra dimension is readout, not state;
    (b) the root needs extent as a primitive. Will's ruling owed.

## Scaffold register entry owed to the paper (v0.6)
BARE-1 / Cl(3): declared. Cargo: SL(2,C), the Lorentz group, the paravector Minkowski
norm, the spinor double cover, the Pauli representation. None of these are results.
Removal route: T1-T3 above. Until run, every c-seat reproduction is a band without an edge.
