# The boundary note -- the readout tier, the horizon and the sky are one boundary

Written 2026-09-05 evening (Claire). Receipt: `bnd1_boundary_receipt.py` (10/10). Source of the identification:
MATH-FOUNDATION-AUDIT sec. 2.9 (`8771dda`), which named it and asked for no new computation; Codex's second audit sec. 3.9,
which supplies the 2+1 caveat. Tier: the identities are DERIVED; the reading of them is a NAMING, and one consequence is a
consistency the fourth direction must pass. Nothing here changes a number.

## 1. What the readout tier already proved, and what it is

The repo has these facts (thm_g, thm_rn, bridge_dbp, graph_cocycle, debt2b), all machine-verified before 09-04:

- the c-seat's reading family is s(u) = cosh lam + sinh lam u with the forced measure du/2 on the sphere of directions,
  and a pivot acts on directions by aberration with du'/du = s^-2 (RN-1; the transport exponent Q_c = 2);
- the hbar-seat's family is q(phi) = cosh 2 lam + sinh 2 lam cos 2 phi with dphi'/dphi = q^-1 (RN-2; Q_hbar = 1);
- the blind moment is the total mass of the source measure; the mirror is I(p, m) = I(-p, 1 - m); the moments are Legendre
  polynomials P_{n-1}(cosh 2 lam); the bath is a flat band and DEBT-2b forces c(omega)^2 g_0(omega) ~ omega^2.

**They are one object.** In the ball model of hyperbolic n-space, the Poisson kernel at hyperbolic distance lam from the
origin, evaluated at a boundary point at angle theta from the displacement, is

    P_n(lam, theta) = (cosh lam - sinh lam cos theta)^-(n-1),

and the boundary Jacobian of the isometry that moves the origin by lam is exactly P_n in the round measure (BND-1 b1a).
With u = -cos theta this is s^-(n-1). So:

- the c-seat's family is the reciprocal Poisson kernel of **H^3** (boundary S^2, the view sphere) and Q_c = 2 = n - 1 (b1b);
- the hbar-seat's family is that of **H^2** (boundary S^1, RP^1 in the doubled angle) and Q_hbar = 1 (b1c);
- the blind-mass theorem, integral of s^Q dmu' = 1, is the statement that the Poisson kernel integrates to one against the
  round measure -- harmonic measure (b2a, b2b);
- the Legendre polynomials are the zonal spherical functions of SO(2,1);
- the mirror (p, m) -> (-p, 1 - m) is the weight duality w <-> Q - w of the boundary (principal-series) representations,
  which graph_cocycle.py had already reduced it to;
- Q-RN ("why is a reading a power of a Radon-Nikodym derivative?") has its answer: a seat's readings are sections of
  homogeneous line bundles over the boundary of its symmetric space; the weight is the representation label; Q is the
  dimension of the boundary. The 08-28 handoff's "structural bet" that Q is a carrier dimension is confirmed.

Kill for the naming: if s^-Q were not the round-measure boundary Jacobian of the pivot on the boundary of H^{Q+1}, this dies.
RN-1 and RN-2 are exactly that check and they pass.

## 2. Two consequences the repo had not drawn

**(i) Q_c = 2 presupposes a SPHERE of directions -- i.e. 3+1.** On the primitives' own 2+1 frame the c-seat's null
directions form a circle and its Jacobian would be s^-1, Q_c = 1, like hbar's. So the transport census (2, 1, 2) is a
consequence of the fourth direction (T8, `d508e64`), not of the bare tier. That is a consistency Gamma_4 must pass, and the
repo can bank it: **any 3+1 completion in which the c-seat's carrier is not the boundary of H^3 fails RN-1.** Codex's second
audit (3.9) sharpens the 2+1 side: uniform du is not uniform angle on S^1 (weight 1/sqrt(1-u^2)) and the presented bath there
has an arcsine density, not a flat band -- so DEBT-2b's flat band is a 3+1 fact, to be recovered once Gamma_4 is earned, not
assumed before.

**(ii) The horizon is the state going to the boundary of the seat's own H^2.** The seat's tilt lam is hyperbolic distance
in its symmetric space: sinh^2 lam = eta (T7d) and q(nu, nu) = -eps/delta = -sech^2 lam (b3a) -- the dual normal's norm reads
the distance. The horizon eta -> oo is lam -> oo (b4a): the state leaves for the boundary, and nu goes null there. But the
seat's READOUT CARRIER is that same boundary (sec. 1). So the sky -- debt2's bath, the boundary of directions the seat reads --
and the horizon -- the boundary the state escapes to -- are **one boundary read two ways**: as the domain of the seat's
readings, and as the limit of its own position. That is the label-free form of "a ruler turns end-on and becomes light"
(P9 sharpened, `8715a9c`): the direction the seat can no longer construct is a point of the boundary it reads.

And Hawking's period lives there without a thermal label. The Wick face (T4, `a36c553`) is the compact real form's orbit inside
the complex group; on the boundary it acts by the complexified boundary action, and the period of that action in the
representation the seat reads (P4: the vector representation, period 2 pi) IS the period T4 computes. The thermal reading
(circumference hbar/k_B T) is a label the paper declares; the period is the model's.

## 3. What the pinch does to this picture

PINCH-1 (`67099dd`) in boundary language: the tilt map state -> H^2 is **discontinuous at the pinch**. Approached along the
horizon, eta = oo and the state is at the boundary (b5b). Approached transversally, eta -> sinh^2 l and the state sits at
hyperbolic distance l from the origin -- an interior point (b5a). One point of state space, two points of the symmetric space,
one interior and one at infinity (b5c). "The horizon has a hole at the pinch through which the seat survives" is this
discontinuity: at the pinch, and only there, a boundary point and an interior point of the seat's own space coincide.

That is also why PINCH-2's resolution parameter xi lives in K(sqrt delta) and not in K (`bab32b7` q5): the tilt lam is a
base-field function of the state (eta is rational), and a base-field function cannot have two limits at one point unless the
point is where the field's own coordinates fail -- the ramification point of the cover. The Kummer coordinate is what
separates the interior limit from the boundary limit.

## 4. What this note does not claim

It does not derive the readout tier from the primitives; it names the theory the readout tier is an instance of. It does
not select Gamma_4 or the 3+1 completion; it states a test any such completion must pass. It does not make the horizon a
gravitational horizon; that is D's job (Codex, rest-frame-trapping theorem). It does not touch the paper's physics sections.

## 5. What it makes available

- v0.6 can replace the cocycle/Mellin formulation of the readout tier by its name, at no cost in results, and gain the
  trichotomy (free / floored / pinned = the three positions of a positive cocycle's zero set relative to an orbit, a
  conjugation class and a base map -- the 08-28 item) and the aggregation-exponent question (a distinguished weight is a
  distinguished representation) as representation-theoretic statements.
- The fourth-direction work (Codex, frame-transport) has one more kill condition: the c-seat's carrier must come out as the
  boundary of H^3 with the round measure.
- The pinch's two limits have a home: interior vs boundary of the seat's H^2, separated by the cover.
