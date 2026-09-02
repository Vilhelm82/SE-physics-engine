# THM-TARGETS K, L, M — the pinning from temperature, rotation, and the swirl as the boost of the pinning
Written 2026-09-03, late.  `thm_k_clausius.py` **16/16** (`5c58e8e`), `thm_l_rotation.py` **10/10** (`d863efa`),
`thm_m_swirl.py` **9/9** (`68f3a9f`).  Read after `docs/2026-09-03-thm-target-J.md`.  Re-running is the verification.

## In one paragraph
THM-J closed the geometric route to the field law.  THM-K opened the thermal one and it went through: the seat's own Unruh
temperature, on a screen with an entropy proportional to its presented area, equilibrated as seen from infinity, returns
the pinning β² = r_s/r exactly with r_s = 2GM/c² falling out.  The exponent is not chosen — the seat's kinematics makes
the redshifted acceleration the gradient of sech²λ and nothing else — and the one fork inside (local vs redshifted
temperature) is decided by Mercury.  THM-L showed the transport law carries vorticity in the only way it can (boost by the
symmetric strain, rotate at half the curl) and reproduces the gyroscope and nodal frame-dragging numbers given the swirl.
THM-M derived the swirl: it is the Lorentz boost of the pinning, superposed over the source, coefficient 2 = 1 (clock)
+ 1 (rods); a scalar theory predicts half and Gravity Probe B kills it.  The gravitational sector now has no imported
profile.  It has one declared idea: a seat's presented sphere carries A/ℓ_P² thermal degrees of freedom.

## THM-K — the chain, with tiers
| step | statement | tier |
|---|---|---|
| K-1 | α = (c²/2)\|dβ²/dr\|/N for any profile; αN = (c²/2)\|d(sech²λ)/dr\| | derived \| KIN-1, RULE-2, H-5, E-4 |
| K-2 | KIN-2a: α = c²r_s/(2r²N); κ = c²/(2r_s) (schedule item 2 closed) | retrodiction |
| K-3 | T_∞ = N·ħα/(2πck_B) = (ħc/4πk_B)\|dβ²/dr\|; Hawking at the wall | derived \| D2 |
| K-4 | dE = T_H dS exact with S = k_B A/(4ℓ_P²): the ¼ and the 2π are one pair | derived \| SCREEN-1 |
| K-5 | equipartition on the presented sphere (area 4πr² proved), blind mass enclosed, Tolman equilibrium ⟹ dβ²/dr = −r_s/r², β² = r_s/r, r_s = 2GM/c² | derived \| SCREEN-1, THERM-1, MASS-1, EQ-1 |
| K-5x/y | local-temperature equilibrium ⟹ lapse harmonic ⟹ perihelion 2.5: dead | failure branch, decided by Mercury |
| K-6 | Poisson for sech²λ with the blind-mass density; Newton's a = GM/r² as N → 1 | derived \| same |

Declared: SCREEN-1, THERM-1, MASS-1, EQ-1.  Lineage at the comparison stage: Jacobson 1995, Verlinde 2010, Padmanabhan,
Tolman–Ehrenfest; the entropic-gravity objections (Kobakhidze 2011 among them) attach to this route and are not answered.

## THM-L — rotation
| step | statement | tier |
|---|---|---|
| L-1 | rule G′ (symmetric strain boost + rotation at ½∇×v⃗) conserves E and L_z identically for any β(r), ω(r); rule G fails; irrotational limit recovered | proved |
| L-2 | with ω = 2GJ/(c²r³): gyroscope at rest precesses at (G/c²r³)[3(J·r̂)r̂ − J]; vorticity a vacuum dipole; polar-orbit average GJ/(2c²r³) = 40.8 mas/yr (GP-B 37.2 ± 7.2) | derived \| KIN-3 (then derived in THM-M) |
| L-3 | node advance 2GJ/(c²a³): ratio 1.0000 to first order at a = 60 r_s, i = 60°; rule G gives zero | derived, numeric corroboration |

## THM-M — the swirl is the boost of the pinning
| step | statement | tier |
|---|---|---|
| M-1 | static presented field at O(r_s): clock h₀₀ = r_s/r, rods h_ij = (r_s/r)r̂_ir̂_j, no shift | derived \| THM-K, H-5 |
| M-2 | BARE-1 sandwich is a Lorentz map (LᵀηL = η); the presented metric is a rank-2 tensor; a moving element's shift is (r_s/r)[w⃗ + (w⃗·r̂)r̂], half clock, half rods | proved |
| M-3 | linear superposition over a rotating shell: far field (2G/c²)(J×r)/r³ to 5×10⁻⁵ at four points ⟹ ω = 2GJ/(c²r³), coefficient 2 = 1 + 1; clock-only gives half: GP-B 20.4 vs 37.2 ± 7.2, dead | derived \| THM-K's inputs, BARE-1, linear superposition |
| M-4 | the natural rotating-screen thermal principles give ∝ J/M with no G: the swirl is not a new thermal input | negative result, recorded |

Scope: first order in r_s and in J.  The full rotating solution is not reached.

## What is now derived, what is declared, what is owed
Derived, conditional on the four thermal declarations and BARE-1: the pinning, the swirl, the geodesic structure null and
timelike to second order, redshift, deflection (2 and 15π/16), perihelion (3), gyroscope dragging (Schiff), nodal dragging
(2), and the finite centre with its two readings.  Declared: SCREEN-1/THERM-1/MASS-1/EQ-1 (one idea), KIN-1, CONT-1,
RULING-2, TRANSPORT-1, STAT-1, ROT-1, PROP-1, linear superposition.  Owed: three rulings (the letter of the rod-free
definition; the semantics of "identity"; whether EQ-1 is a principle or Mercury's fact) and the comparison with the
entropic-gravity objections.  Next tests with kill numbers: rotation beyond first order (the full twisting river; the
Kerr quadrupole), the Hubble river (uniform source), and the neutron-interferometry objection to the thermal route.
