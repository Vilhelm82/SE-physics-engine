# CURV-1 — predictions written before any path runs
Not for the runner. Opened by the comparison stage after the suites are green.

## Will (2026-09-02, made in the previous session before the relevant prior work surfaced)
W-1  Omega goes imaginary inside, because det G' < 0 there.
W-2  The centre holonomy is finite and non-trivial: specifically the deck Z2.
W-3  The horizon crossing is a fixed quarter turn.

Will also holds earlier, independent work on this exact question from a different
corpus. It is not in this folder and enters only at the comparison stage.

## Claire (2026-09-02, hand derivation from D-19/D-20 and the centre Gram; not run; may be wrong)
Status of each: derived-by-hand | conditional on the section-2 conventions of PATH1-SPEC.

C-1  Reality inside splits by layer. Option A (K = root): the root's presented rod is
     rho_3 = r/(r - r_s) < 0 on the whole interior while the two visible rods are
     positive below their null radii, so the root's presented direction carries a
     factor i, S_eps is complex, V' is real, and Omega is complex with both parts
     nonzero — not purely imaginary. Option C, below all three null radii: all rods
     positive, det G' < 0, so V' is imaginary and S_eps real: Omega is purely
     imaginary or pi + i(real). W-1 holds on the world layer deep inside and fails
     on the seat layer.
C-2  Centre, option A: G'(0) = G - k k^T with k = (g13, g23, 1). The root's presented
     vector vanishes, Gam_13, Gam_23 -> 0 like sqrt(r), Gam_12 -> (g12 - g13 g23)/
     (sin th13 sin th23) = cos phi_3, the frame's dihedral angle at the root, and
     V' -> sin phi_3 (this is D-20). So the orbit -> {phi_3, phi_3, pi - phi_3, pi - phi_3}:
     finite, nontrivial, real. Picture: the root goes to the pole, the two visible
     axes to the equator; the presented triangle is the polar triangle whose excess
     is its apex angle. The seat's rod-free curvature reading at the centre is the
     frame's own dihedral angle at the root.
C-3  Centre, option C: G'(0) is the Gram of the projections p_i = a_i - k_i K onto the
     plane perpendicular to K, and sum_i p_i = 0 because K ~ a1+a2+a3. Three plane
     vectors summing to zero surround the origin, so S_(+++) = -4 prod cos(gap/2) < 0
     and V' = 0: the orbit -> {2 pi, 0, 0, 0}. Identity in the rotation, one member
     -1 in the spin lift: the deck Z2, on the world layer only. Agrees with W-2 for
     the hole; the infalling seat (C-2) sees an angle instead.
C-4  Horizon: when all k_i have one sign the (+++) member -> 0 from outside; any
     mixed-sign member is a 0/0 there and its limit is whatever the family gives.
     No value predicted.
C-5  (added 2026-09-02 afternoon, before any run) Layer A: the two visible presented directions
     collide, G'_12^2 = rho_1 rho_2, at r* = r_s Delta/(1 - g12^2), below both null radii.
     Half the orbit sees the pair antipodal there; those members' S + iV vanish.

## Outcome
See CURV1-RESULTS.md, "Predictions scorecard".  W-2, W-3, C-2, C-3, C-5 held; W-1 held on the world
layer only; C-1 held; C-4 made no claim.
