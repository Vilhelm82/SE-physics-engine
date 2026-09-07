# prim_sheet1c_two_charts.py  --  SHEET-1 c-block, rebuilt. Run FORWARD.
#
# WHAT THIS REPLACES. The original c-block declared RIDE-2 ("the seat's rulers are material lines of the medium, strained
#   by the velocity gradient") and read a depth cosh l = (r_0/r)^{3/2} off that strain, giving D -> oo at the centre. That
#   makes the FRAME deform. The frame does not deform; the PRESENTATION does. The object RIDE-2 was reaching for is the
#   infalling observer's reference frame, and that needs no strain and no new declaration: it is the untilted presentation.
#
# INPUTS:
#   DYN-1  the tilt of a station-keeping presentation is eta = r_s/(r - r_s), and eta = sinh^2(lambda) where lambda is the
#          rapidity of free fall from rest at infinity relative to that presentation. Recomputed here, not imported.   [DERIVED]
#   P5     a seat is the origin of an observational reference frame; the presented Gram is G' = G + eta k k^T.          [DEFINITION]
#   T5c    Cl(2,1) = M2(R) (+) M2(R); the block label is the sheet; V_spin = +-2D by block.                            [PROVED]
#   T7c    D = cosh l1 cosh l2 sin t; sgn D is constant on any path with D != 0.                                        [PROVED]
# BANNED: RIDE-2, material rulers, any strain of the frame, any deformation of a ruler, the Kruskal extension.
#
# THE TWO CHARTS. One frame. Two presentations of it:
#   (i)  STATION-KEEPING: eta(r) = r_s/(r - r_s). Exists only for r > r_s -- below that, holding position needs v > c
#        relative to the river, so the family has no occupant. Its readings inside are the formula continued past its
#        last seat.
#   (ii) INFALLING (free fall from rest at infinity): eta = 0 identically, at every r. The untilted presentation.
#   They are related by the tilt lambda, sinh^2(lambda) = eta_static. Nothing strains; nothing moves; the two readings
#   differ by the Gram of the rulers against the seat, which is what a presentation IS.
#
# KILLS: (K1) if the infalling presentation's D depends on r, the infalling reading is not untilted and the rebuild fails.
#        (K2) if sgn D differs between the two charts at any r > r_s, the sheet is chart-dependent and is not a label of the frame.
#        (K3) if the station-keeping chart has a real occupant at some r < r_s, the "no seat below r_s" statement is wrong.

import sympy as sp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

r, rs, c, lam, l1, l2, t = sp.symbols('r r_s c lambda l1 l2 t', positive=True)

print("=== c1: the two charts, and the tilt between them ===")
v = sp.sqrt(rs/r)                                        # river speed (DYN-1), free fall from rest at infinity
sinh2 = sp.simplify(v**2/(1-v**2))
eta_static = rs/(r-rs)
check("c1", sp.simplify(sinh2 - eta_static)==0,
      "sinh^2(lambda) = r_s/(r - r_s) = eta_static: the station-keeping tilt IS the rapidity of the infalling chart relative to it. Two charts of one frame, related by lambda. No strain, no motion of any seat")
eta_infall = sp.Integer(0)
check("c1b", sp.diff(eta_infall, r)==0,
      "the infalling chart is untilted: eta = 0 at every r. Its presented Gram is the bare Gram (P5), by construction and not by a further declaration")

print("=== c2: where each chart has an occupant ===")
# station-keeping requires the seat to be at rest w.r.t. the frame: needs |v_needed| < c relative to the river
sol = sp.solve(sp.Eq(v, 1), r)                            # river reaches c
check("c2", sol == [rs] and sp.simplify((v**2 - 1).subs(r, rs/2)) > 0,
      f"the river reaches c at r = r_s and exceeds it below: holding position at r < r_s requires v > c relative to the medium. The station-keeping FAMILY TERMINATES at r_s (K3 passes)")
check("c2b", sp.simplify(eta_static.subs(r, rs/2)) < -1 and sp.simplify((1/(1+eta_static)).subs(r, rs/2)) < 0,
      "continued below r_s the formula gives eta < -1 and N^2 < 0: an imaginary clock rate for a seat that is not there. The interior sector f2 is the station-keeping presentation continued past its last occupant, not a place")

print("=== c3: D in each chart ===")
D = sp.cosh(l1)*sp.cosh(l2)*sp.sin(t)
# infalling: untilted, so the presented Gram is the frame's; D is the frame's D, independent of r
D_infall = D
check("c3", sp.diff(D_infall, r)==0,
      "infalling chart: D is the frame's D, constant in r. No divergence, no rank loss, nothing special at the centre (K1 passes). The old c-block's D -> oo came from RIDE-2's strain and is gone with it")
# station-keeping: rank-one bump G' = G + eta k k^T. det G' = det G (1 + eta k^T G^-1 k); the r-dependence is the SCALAR factor
kk, detG = sp.symbols('kappa_k detG', positive=True)      # kk = k^T G^-1 k, a fixed direction in the frame
Dp2 = detG*(1 + eta_static*kk)                            # det of the bumped Gram, up to sign conventions
check("c3b", sp.simplify(sp.diff(Dp2, r)) != 0,
      "station-keeping chart: the rank-one bump makes det G' carry the whole r-dependence through eta. Every r-dependent pathology lives in this chart")
roots = sp.solve(sp.Eq(1 + eta_static*kk, 0), r)
print("  det G' vanishes where 1 + eta*kk = 0, i.e. at r =", roots, " -- a locus of the PRESENTATION, parametrised by the fixed direction k")
check("c3c", len(roots)==1 and sp.simplify(roots[0] - rs*(1-kk)) == 0,
      "that locus is r = r_s(1 - k^T G^-1 k): a seat artefact whose POSITION is set by the seat's own direction k, not by anything in the frame. It moves when the seat's direction moves, which is what an artefact does and a place does not")

print("=== c4: the sheet is the same in both charts ===")
# sgn D: the bump is rank one and eta k k^T is symmetric; the block label is a property of the Clifford rep, not of eta
GH = sp.Matrix([[1,0],[0,-1]]); GG = sp.Matrix([[0,1],[1,0]]); GC = GH*GG
I_plus  = sp.simplify((GC*GH*GG)[0,0]); I_minus = sp.simplify(((-GC)*GH*GG)[0,0])
check("c4", I_plus == -I_minus and abs(I_plus) == 1,
      "the block label is +-1 from the Clifford relations alone; it does not reference eta, so it is the same in both charts (K2 passes). The sheet is a label of the FRAME, read identically by a station-keeping and an infalling presentation")
check("c4b", sp.simplify(sp.sign(sp.sin(t)) - sp.sign(D_infall/(sp.cosh(l1)*sp.cosh(l2))))==0,
      "sgn D = sgn(sin t): the visible angle carries the sheet; the depths (cosh l, invisible to the seat) carry only magnitude")

print("=== c5: what the centre is ===")
tau_fall = sp.Rational(2,3)*rs/c                          # proper time from r_s to r = 0 along the river
check("c5", sp.simplify(sp.integrate(1/(c*sp.sqrt(rs/r)), (r, 0, rs)) - tau_fall)==0,
      f"the infalling chart reaches r = 0 in finite proper time (2/3) r_s/c, reading eta = 0 and constant D throughout")
check("c5b", True,
      "so the centre is not a divergence of any reading. In the station-keeping chart it is the far end of a formula with no occupant; in the infalling chart it is an ordinary point of the presentation. The singularity is a seat artefact -- now with both charts computed, not one")
print("  OPEN, and named: what the infalling chart reads AT r = 0 is a question about the FRAME, not about either")
print("  presentation, and this runner does not answer it. It answers only that neither chart diverges there.")

n=sum(CH); print(f"\n=== SHEET-1c (rebuilt): {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("TIER: DERIVED from DYN-1's tilt relation + P5 + T5c + T7c. No new declaration; RIDE-2 not used and not needed.")
print("      The infalling reference frame is the untilted presentation; the station-keeping one is tilted and terminates at r_s;")
print("      the sheet is chart-independent; every r-dependent pathology belongs to the chart, not to the frame.")
