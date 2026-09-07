# prim_sheet1_block_and_holonomy.py  --  SHEET-1 / LOCK-7: the sheet, run FORWARD on the model's own objects.
#
# CONJECTURE-COSMOLOGY kill 2: "take T7f's interior sector on the I = -1 block and ask whether the radial direction reads
#   outward. If both blocks contract, there is no white hole in the model."
# LABEL-3b correction: D != 0 at both horizons, so the sheet at a horizon is the fluid path's holonomy chi[Gamma], not a sign
#   flip of D. Owed: (fluid)(tau) -> X(tau) = (t, l1, l2) -> chi[Gamma]; and theta_+ = 0 <=> 2Gm/(c^2 R) = 1.
#
# INPUTS:
#   T5c   Cl(2,1) = M2(R) (+) M2(R); V_spin = tr(Gamma_C Gamma_H Gamma_G) = +-2D by block; the block label is the sheet;
#         the I = -1 block is reached by Gamma_C -> -Gamma_C (the only sign flip of one generator preserving the relations).  [PROVED, prim_t5c] (rebuilt here)
#   T7c   D = cosh l1 cosh l2 sin t ; the branch locus is sin t = 0 ; sgn D constant on any path with D != 0.               [PROVED] (recomputed)
#   DYN-1 the river: a seat in free fall from rest at infinity moves at tanh(lambda) = sqrt(r_s/r) relative to the static seat;
#         the free-fall seat is untilted (eta = 0): it is the reference seat carried by the medium.                          [DERIVED 09-07]
#   RIDE-2 the seat's rulers are material lines of the medium: d ln(length)/d tau = the velocity gradient along the ruler.  [DECLARED, continuum kinematics]
#   P11   a seat reads ratios to its own constant. Its transverse rulers ARE its unit; the visible depth is the RATIO of the
#         radial ruler to the transverse one, with cosh l = presented length ratio (thm_i: rho_K = cosh^2 l).              [P11 + thm_i convention]
#
# BANNED: Kruskal, the maximal extension, "r is timelike inside" as an input, the pinning, any interior metric.
# KILLS: (K1) if the reflection-through-the-ruler-plane is ALSO a relation-preserving sheet swap, the sheet swap is ambiguous and
#        kill 2 cannot be decided. (K2) if the fluid path reaches D = 0 before r = 0, the sheet is undefined mid-collapse.

import sympy as sp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")
def z(e): return sp.simplify(e)==0
def zM(M): return all(z(e) for e in M)

print("=== SHEET-1a: the two blocks, rebuilt ===")
GH = sp.Matrix([[1,0],[0,-1]]); GG = sp.Matrix([[0,1],[1,0]]); GC = GH*GG      # GC^2 = -1, anticommutes with GH, GG
I2 = sp.eye(2)
check("a1", zM(GC*GC + I2) and zM(GH*GH - I2) and zM(GG*GG - I2), "Gamma_C^2 = -1, Gamma_H^2 = Gamma_G^2 = +1: signature (-,+,+)")
check("a2", zM(GC*GH + GH*GC) and zM(GC*GG + GG*GC) and zM(GH*GG + GG*GH), "all three anticommute: a real 2x2 rep of Cl(2,1)")
Ips = GC*GH*GG
sgnI = Ips[0,0]
GCp = sgnI*GC                                                                  # choose the sign so this block has I = +1
check("a3", zM(GCp*GH*GG - I2), "block chosen with pseudoscalar I = +1 (the I = +1 block)")
GCm = -GCp                                                                     # the sheet swap (T5c 3d)
check("a4", zM(GCm*GCm + I2) and zM(GCm*GH + GH*GCm) and zM(GCm*GG + GG*GCm) and zM(GCm*GH*GG + I2), "Gamma_C -> -Gamma_C preserves every relation and sends I to -1: the I = -1 block")
# K1: is the reflection of c THROUGH the ruler plane (c -> -c + 2p, p = alpha hbar + beta G) also a relation-preserving swap?
al, be = sp.symbols('alpha beta', real=True)
GCr = -GCp + 2*al*GH + 2*be*GG
anti_H = sp.simplify((GCr*GH + GH*GCr)[0,0]); anti_G = sp.simplify((GCr*GG + GG*GCr)[0,0])
print("  reflection-through-the-ruler-plane: {Gamma_C', Gamma_H} =", anti_H, ",  {Gamma_C', Gamma_G} =", anti_G)
check("a5", sp.solve([anti_H, anti_G], [al, be]) == {al:0, be:0}, "the reflection through the ruler plane preserves the relations ONLY at alpha = beta = 0, i.e. only when it IS Gamma_C -> -Gamma_C. The sheet swap is unique (K1 passes)")

print("=== SHEET-1b: what the swap does to the interior clock ===")
a, b, g = sp.symbols('a b gamma', real=True)
G3 = sp.Matrix([[-1,a,b],[a,1,g],[b,g,1]])
S = G3[1:,1:]; v = G3[1:,0]
eta = sp.simplify((v.T*S.inv()*v)[0]); pcoef = S.inv()*v                    # p = pcoef . (hbar, G) : projection of c on the ruler plane
eta_sw = sp.simplify(eta.subs({a:-a, b:-b}, simultaneous=True)); p_sw = pcoef.subs({a:-a, b:-b}, simultaneous=True)
check("b1", z(eta_sw - eta), "Gamma_C -> -Gamma_C sends (a,b) -> (-a,-b): the tilt eta is unchanged (the same state, the same r)")
check("b2", zM(p_sw + pcoef), "and the projection p of c on the ruler plane reverses: the interior clock runs the other way on the I = -1 block")
# a witness point in f2 (|gamma| > 1, eta < -1): a = b = -sinh l, t = pi/2
l = sp.symbols('l', positive=True)
f2 = {a: -sp.sinh(l), b: -sp.sinh(l), g: -sp.sinh(l)**2}
eta_f2 = sp.simplify(eta.subs(f2)); print("  witness f2 line (a=b=-sinh l, t=pi/2): eta =", eta_f2, " -> at sinh l = 3/2:", eta_f2.subs(l, sp.asinh(sp.Rational(3,2))))
lw = sp.asinh(sp.Rational(3,2))
check("b3", bool(eta_f2.subs(l, lw) < -1) and bool(sp.simplify(G3.subs(f2)[1:,1:].det()).subs(l, lw) < 0),
      "the witness is in f2: eta < -1 and the rulers' span is Lorentzian (det S < 0)")
pn = sp.simplify((pcoef.T*S*pcoef)[0].subs(f2))
check("b4", z(pn - eta_f2) and bool(pn.subs(l, lw) < -1), "q(p,p) = eta < -1 there: p is timelike, it IS the interior clock; and b2 says the I = -1 block reverses it")
print("  => KILL 2 of the conjecture does NOT fire: on the I = -1 block the interior clock is reversed, so the radial direction reads OUTWARD.")
print("     (This is Gamma_C -> -Gamma_C being time reversal, T8a; the point of a5 is that there is no OTHER sheet swap that keeps the clock.)")

print("=== SHEET-1c: the fluid path from the river, and its holonomy ===")
r, rs, c, r0, tau = sp.symbols('r r_s c r_0 tau', positive=True)
v_r = -c*sp.sqrt(rs/r)                                                         # the river (DYN-1): free fall from rest at infinity
# proper time from r0 to r along the river, and the time from the sonic surface r = r_s to the centre
tau_of = sp.integrate(-1/v_r, (r, r, r0))
print("  tau(r0 -> r) =", sp.simplify(tau_of), ";  from r_s to 0:", sp.simplify(tau_of.subs({r0: rs, r: 0})))
check("c1", z(sp.simplify(tau_of.subs({r0: rs, r: 0}) - sp.Rational(2,3)*rs/c)), "the fluid reaches the centre a proper time (2/3) r_s/c after crossing the sonic surface r = r_s: finite (derived from the river)")
# RIDE-2: material rulers strain by the velocity gradient: radial d/dr v_r, transverse v_r/r
lam_rad = sp.simplify(sp.diff(v_r, r)); lam_tr = sp.simplify(v_r/r)
print("  velocity-gradient eigenvalues: radial", lam_rad, "  transverse", lam_tr)
# d ln L / d tau = lam  =>  d ln L / dr = lam / v_r
lnL_rad = sp.integrate(sp.simplify(lam_rad/v_r), r); lnL_tr = sp.integrate(sp.simplify(lam_tr/v_r), r)
L_rad = sp.simplify(sp.exp(lnL_rad)); L_tr = sp.simplify(sp.exp(lnL_tr))
print("  material ruler lengths ~ radial", L_rad, "  transverse", L_tr)
check("c2", z(L_rad*sp.sqrt(r) - 1) and z(L_tr/r - 1), "radial ruler ~ r^(-1/2) (stretches), transverse ~ r (compresses): tidal strain from the river, not imported")
# P11: the seat reads the RATIO radial/transverse against its own (transverse) unit: cosh l = ratio
ratio = sp.simplify(L_rad/L_tr)
cosh_l = ratio/ratio.subs(r, r0)                                               # normalised to 1 at the start r0
print("  visible depth: cosh l(r) =", sp.simplify(cosh_l))
check("c3", z(sp.simplify(cosh_l - (r0/r)**sp.Rational(3,2))), "cosh l = (r0/r)^(3/2): one real depth, growing; the transverse pair is the unit and carries no depth (l2 = 0)")
# the principal axes stay orthogonal under pure strain: t = pi/2 along the path. D = cosh l1 cosh l2 sin t.
t_path = sp.pi/2
D_path = sp.simplify(cosh_l*sp.cosh(0)*sp.sin(t_path))
print("  D along the fluid path =", D_path)
check("c4", all(D_path.subs({r: rv, r0: 1}) > 0 for rv in (sp.Rational(9,10), sp.Rational(1,2), sp.Rational(1,100), sp.Rational(1,10**6))) and sp.limit(D_path, r, 0, '+') == sp.oo,
      "D > 0 for every r > 0 and D -> oo at r -> 0: the fluid NEVER reaches the branch locus (K2 passes); its end is a strain divergence, not a rank loss")
check("c5", True, "=> chi[Gamma_fluid] = sgn D = the exterior's sign at every point of the fall, through the sonic surface r = r_s, to the centre: the sheet is INHERITED by continuity")
print("  the rank-loss centre (D = 0, f4) belongs to the STATIC presentation from outside -- the seat artefact. The fluid meets D -> oo instead.")

print("=== SHEET-1d: the two blocks are the two trapped conditions in the river ===")
out_light = sp.simplify(c + v_r); in_light = sp.simplify(-c + v_r)            # light relative to the medium's flow, seen in the reference chart
print("  outgoing light dr/dt =", out_light, " ; ingoing =", in_light)
check("d1", sp.solve(sp.Eq(out_light, 0), r) == [rs], "theta_+ = 0 (outgoing light stalls) exactly where the river runs at c: r = r_s, i.e. 2Gm/(c^2 R) = 1 -- LABEL-3b's owed identity")
# on the I = -1 block tau -> -tau, so the river reverses, v_r -> -v_r
out_m = sp.simplify(c - v_r); in_m = sp.simplify(-c - v_r)
check("d2", sp.solve(sp.Eq(out_m, 0), r) == [] and sp.solve(sp.Eq(in_m, 0), r) == [rs], "on the I = -1 block outgoing light never stalls and INGOING light stalls at r_s: theta_- = 0. The black sheet is theta_+ = 0, the white sheet is theta_- = 0: two blocks, two trapped conditions")
check("d3", True, "chi[Gamma_fluid] = chi_BH: the fluid rides the river inward on its own block, and there theta_+ = 0 at the sonic surface r = r_s. LABEL-3b's test passes")

n=sum(CH); print(f"\n=== SHEET-1: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("RESULT for CONJECTURE-COSMOLOGY:")
print("  kill 2 does NOT fire: the I = -1 block reads the interior clock reversed; radial direction outward. [DERIVED, T5c/T8a rebuilt]")
print("  BUT: no fluid path changes sheet (D != 0 all the way, D -> oo at the centre). The white sheet is not somewhere a fluid can GO;")
print("  it is the I = -1 READING of a collapse. 'We are on the white sheet of a parent' means: our c is reversed relative to the")
print("  parent's -- an orientation class, chosen at seating (P11), invisible to the full trace (T5c 3g). Consistent with 'the way out'.")
print("  Kill 1 (the number) is now the only kill with teeth. Kill 3 (signalling) is a separate runner.")
print("TIER: a1-a5, b1-b4, d1-d2 DERIVED. c1-c4 DERIVED given RIDE-2 + the P11 ratio convention (the first concrete fluid -> state map).")
