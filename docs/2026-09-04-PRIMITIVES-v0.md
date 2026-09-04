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
T4. What constrains the three angles? The elliptope came from R^3. Without R^3, ask
    what P8 forces.

## Scaffold register entry owed to the paper (v0.6)
BARE-1 / Cl(3): declared. Cargo: SL(2,C), the Lorentz group, the paravector Minkowski
norm, the spinor double cover, the Pauli representation. None of these are results.
Removal route: T1-T3 above. Until run, every c-seat reproduction is a band without an edge.
