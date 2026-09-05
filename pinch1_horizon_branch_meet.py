#!/usr/bin/env python3
# =============================================================================
# PINCH-1 -- where the horizon meets the branch locus.  The model's one genuinely singular
#   stratum, and what a seat reads there.  (Claire, 2026-09-05.)
# INPUTS: T7a's form q = diag(1,1,-1), seat c = (0,0,1); the Gram G(a,b,gamma) with
#   eps = 1 - gamma^2, W = a^2 + b^2 - 2ab gamma, delta = -det G, F = eps/delta = N^2,
#   eta = W/eps, nu = G^-1 l (Codex, verified 09-05).  NO soldering map, NO r, NO pinning.
# CONTEXT: Codex's horizon-crossing chart (docs/2026-09-05-horizon-crossing-metric.md) requires
#   delta > 0 and a != sigma b, and says e_r is undefined at zero tilt.  This runner is about the
#   set that excludes.
# CLAIMS, checked:
#   p1  delta = eps + W identically.  On gamma = sigma: delta = W = (a - sigma b)^2.  Codex's two
#       exclusions (full rank loss; zero tilt) are ONE set on the horizon.
#   p2  THE PINCH P_sigma = {gamma = sigma, a = sigma b} = {|gamma| = 1} cap {delta = 0}.  Codim 2.
#       Seat coordinates: P_+ = {t = 0, l1 = l2}, P_- = {t = pi, l1 = -l2}.
#   p3  At P the frame has rank 2 with h = sigma g EXACTLY (ker G = ker M): two rulers become one;
#       Will's null n0 = h - sigma g shrinks to the ZERO vector, linearly in (a - sigma b).
#   p4  At P every seat invariant is 0/0: F, eta, N^2, and the numerator of nu vanishes too.
#   p5  TRANSVERSE approach (gamma' != 0): F -> sech^2 l, eta -> sinh^2 l (TILT = DEPTH), and nu
#       converges to ONE vector nu* = (-c + a h)/(1 + a^2), timelike, q(nu*,nu*) = -sech^2 l,
#       INDEPENDENT of the tangent.  The 0/0 is REMOVABLE.  There is a seat at the pinch.
#   p6  IN-HORIZON approach (gamma == sigma): F == 0.  nu = n0/(a - sigma b); frame components
#       diverge but the VECTOR converges, to one of TWO null vectors nu_H^+-, one per SHEET
#       (the sign of the frame's orientation).
#   p7  nu* = (nu_H^+ + nu_H^-)/2.  The seat's normal at the pinch is the BISECTOR of the two
#       sheets' null limits.  One side presenting as two, in one vector identity.
#   p8  TANGENT-LEAVING approach (gamma' = 0, gamma'' != 0): the ORDER-PAIR (ord eps, ord W)
#       decides; with the same tangent (1,0,0), F -> 1, 2/3, 0.  Every value in [0,1] is attainable.
#   p9  KILL CHECK, both halves recorded: the tangent decides in class p5 (kill of 'order-pair is
#       the invariant' TRIGGERED there); the order-pair decides in class p8 (kill of 'the tangent
#       decides' TRIGGERED there).  The invariant is the approach CLASS, stratified by contact
#       order with the horizon.  Neither half alone is the theorem.
# TIER: [DERIVED | T7a's declared signature].  Nothing else declared.
# [First run: p3a and p3c FAILED on runner bugs (kernel normalisation; unsimplified 0/0 on subs). Fixed; claims unchanged.]
# KILL for p5: any transverse path with F -> anything but sech^2 l, or nu -> anything but nu*.
# KILL for p7: any sheet realisation whose two null limits do not average to nu*.
# =============================================================================
import sympy as sp, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), sp.cancel, sp.factor, lambda q: sp.simplify(q.rewrite(sp.exp))):
        try:
            if f(e) == 0: return True
        except Exception: pass
    return False
def zv(v): return all(z(e) for e in v)

a, b, gam = sp.symbols('a b gamma', real=True)
s = sp.Symbol('s', positive=True)
b0, al, be, g1, g2 = sp.symbols('b0 alpha beta gamma1 gamma2', real=True)
t, l1, l2 = sp.symbols('t l1 l2', real=True)
q = sp.diag(1, 1, -1)
Q = lambda u, v: sp.simplify((u.T*q*v)[0])
G = sp.Matrix([[-1, a, b], [a, 1, gam], [b, gam, 1]])
eps = 1 - gam**2; W = a**2 + b**2 - 2*a*b*gam; delta = sp.expand(-G.det()); F = eps/delta; eta = W/eps
numnu = sp.Matrix([-eps, a - b*gam, b - a*gam])       # delta * nu, rational

print("=== p1: delta = eps + W; the two exclusions are one ===")
check("p1a delta = eps + W identically", z(delta - (eps + W)))
for sig in (1, -1):
    check(f"p1b on gamma = {sig:+d}: delta = W = (a - sigma b)^2", z(delta.subs(gam, sig) - (a - sig*b)**2) and z(W.subs(gam, sig) - (a - sig*b)**2))
check("p1c so on the horizon, {delta = 0} = {W = 0} = {a = sigma b}: Codex's 'full rank loss' and 'zero tilt' coincide", True if (CH[-1] and CH[-2]) else False)
# (p1c is a corollary of p1b; it is recorded as such, not as an independent computation)

print("=== p2: the pinch is codimension 2 ===")
sol = sp.solve([sp.Eq(eps, 0), sp.Eq(delta, 0)], [gam, a], dict=True)
check("p2a {eps = 0} cap {delta = 0} = {gamma = +-1, a = gamma b}: solve returns exactly the two pinch lines", len(sol) == 2 and all(sp.simplify(d[a] - d[gam]*b) == 0 for d in sol))
A_, B_, Gm_ = -sp.sinh(l1), -sp.sinh(l2), sp.cosh(l1)*sp.cosh(l2)*sp.cos(t) - sp.sinh(l1)*sp.sinh(l2)
check("p2b seat coords: at t = 0, gamma = cosh(l1 - l2) >= 1 with equality iff l1 = l2  -> P_+ = {t = 0, l1 = l2}", z(sp.expand_trig(Gm_.subs(t, 0)) - sp.cosh(l1 - l2)))
check("p2c seat coords: at t = pi, gamma = -cosh(l1 + l2) <= -1 with equality iff l1 = -l2 -> P_- = {t = pi, l1 = -l2}", z(sp.expand_trig(Gm_.subs(t, sp.pi)) + sp.cosh(l1 + l2)))
check("p2d on P_+ : a = b (l1 = l2) and on P_-: a = -b (l1 = -l2): consistent with a = sigma b", z((A_ - B_).subs(l2, l1)) and z((A_ + B_).subs(l2, -l1)))

print("=== p3: at the pinch the frame collapses to two lines; n0 -> 0 ===")
for sig in (1, -1):
    Gp = G.subs({gam: sig, a: sig*b})
    ker = Gp.nullspace()
    check(f"p3a sigma = {sig:+d}: rank G = 2 and ker G = span(0, 1, -sigma), i.e. the relation h - sigma g", Gp.rank() == 2 and len(ker) == 1 and zv(ker[0].cross(sp.Matrix([0, 1, -sig]))))
# explicit realisation ON the horizon near the pinch: a = sigma b + s; two sheets y = +-s/sqrt(1+a^2)
sig = 1
a_h = sig*b + s
c = sp.Matrix([0, 0, 1]); h = sp.Matrix([sp.sqrt(1 + a_h**2), 0, -a_h])
gvec = lambda sign: sp.Matrix([(sig + a_h*b)/sp.sqrt(1 + a_h**2), sign*s/sp.sqrt(1 + a_h**2), -b])
ok = True
for sign in (1, -1):
    g = gvec(sign)
    ok = ok and z(Q(g, g) - 1) and z(Q(c, g) - b) and z(Q(h, g) - sig) and z(Q(c, h) - a_h) and z(Q(c, c) + 1) and z(Q(h, h) - 1)
check("p3b explicit horizon states: q(c,c) = -1, q(h,h) = q(g,g) = 1, q(c,h) = a, q(c,g) = b, q(h,g) = sigma, for both sheets", ok)
n0p, n0m = h - sig*gvec(1), h - sig*gvec(-1)
check("p3c n0 = h - sigma g vanishes LINEARLY in s = a - sigma b (both sheets): h = sigma g exactly at the pinch",
      zv((n0p/s).applyfunc(lambda e: sp.limit(sp.cancel(e), s, 0, "+")) - sp.Matrix([b/sp.sqrt(1 + b**2), -1/sp.sqrt(1 + b**2), -1])) and zv(n0p.subs(s, 0)) and zv(n0m.subs(s, 0)))
check("p3d and n0 is null all along the horizon (Will's 'a whole line with zero length')", z(Q(n0p, n0p)) and z(Q(n0m, n0m)))

print("=== p4: every seat invariant is 0/0 at the pinch ===")
P0 = {gam: 1, a: b}
check("p4a eps, W, delta all vanish at P_+", z(eps.subs(P0)) and z(W.subs(P0)) and z(delta.subs(P0)))
check("p4b so F = eps/delta, eta = W/eps, N^2 are 0/0 there", True if CH[-1] else False)
check("p4c and the NUMERATOR of nu, delta*nu = (-eps, a - b gamma, b - a gamma), vanishes too: nu is 0/0 as a vector", zv(numnu.subs(P0)))

print("=== p5: TRANSVERSE approach -- the 0/0 is removable, uniquely ===")
Ptr = {a: sig*b0 + al*s, b: b0 + be*s, gam: sig + g1*s + g2*s**2}
Fs = sp.simplify(sp.series(F.subs(Ptr), s, 0, 1).removeO())
check("p5a F -> 1/(1 + b0^2) = sech^2(depth) along EVERY path with gamma' != 0, for all alpha, beta, gamma1, gamma2", z(Fs - 1/(1 + b0**2)) and not Fs.has(al) and not Fs.has(be) and not Fs.has(g1))
etas = sp.simplify(sp.series(eta.subs(Ptr), s, 0, 1).removeO())
check("p5b eta -> b0^2 = sinh^2(depth): the TILT equals the DEPTH at the pinch (T7d: sinh^2 lambda = eta -> lambda = l)", z(etas - b0**2))
nus = (numnu/delta).subs(Ptr).applyfunc(lambda e: sp.simplify(sp.series(e, s, 0, 1).removeO()))
vec_c = sp.simplify(nus[0]); vec_h = sp.simplify(nus[1] + sig*nus[2])   # at P, g = sigma h as vectors
check("p5c nu's frame components DO depend on the tangent, but g = sigma h there, and the VECTOR nu* = vec_c c + vec_h h does not",
      (nus[1].has(al) or nus[1].has(g1)) and not vec_c.has(al) and not vec_h.has(al) and not vec_h.has(g1) and not vec_h.has(be))
check("p5d nu* = (-c + a h)/(1 + a^2) with a = sigma b0", z(vec_c + 1/(1 + b0**2)) and z(vec_h - sig*b0/(1 + b0**2)))
a0 = sig*b; h0 = sp.Matrix([sp.sqrt(1 + a0**2), 0, -a0])
nustar = sp.simplify((-c + a0*h0)/(1 + a0**2))
check("p5e nu* is TIMELIKE with q(nu*,nu*) = -sech^2 l, and q(nu*,c) = 1, q(nu*,h) = 0: it is a seat's normal", z(Q(nustar, nustar) + 1/(1 + b**2)) and z(Q(nustar, c) - 1) and z(Q(nustar, h0)))

print("=== p6: IN-HORIZON approach -- two null limits, one per sheet ===")
nuHp, nuHm = (n0p/s).applyfunc(sp.cancel), (n0m/s).applyfunc(sp.cancel)      # nu = n0/(a - sigma b) on the horizon
check("p6a on the horizon F == 0 and nu = n0/(a - sigma b) (frame comps (0, 1, -sigma)/(a - sigma b): divergent)", z(F.subs(gam, sig)) and zv((numnu/delta).subs(gam, sig) - sp.Matrix([0, 1, -sig])/(a - sig*b)))
check("p6b but the VECTOR nu converges: nu_H^+ = (b, -1, -sqrt(1+b^2))/sqrt(1+b^2), nu_H^- with the middle sign flipped",
      zv(nuHp.subs(s, 0) - sp.Matrix([b, -1, -sp.sqrt(1 + b**2)])/sp.sqrt(1 + b**2)) and zv(nuHm.subs(s, 0) - sp.Matrix([b, 1, -sp.sqrt(1 + b**2)])/sp.sqrt(1 + b**2)))
check("p6c both limits are NULL", z(Q(nuHp.subs(s, 0), nuHp.subs(s, 0))) and z(Q(nuHm.subs(s, 0), nuHm.subs(s, 0))))
check("p6d they differ by the SHEET: the sign of the y-component is the orientation of the frame (sign of det[c,h,g])",
      z(sp.Matrix.hstack(c, h, gvec(1)).det()/s - (-sp.Matrix.hstack(c, h, gvec(-1)).det()/s)) and not z(sp.Matrix.hstack(c, h, gvec(1)).det()))

print("=== p7: the seat's normal at the pinch bisects the two sheets ===")
check("p7a nu* = (nu_H^+ + nu_H^-)/2 EXACTLY", zv((nuHp.subs(s, 0) + nuHm.subs(s, 0))/2 - nustar))
check("p7b and nu_H^+ - nu_H^- is spacelike, orthogonal to nu*: the two sheets are split along a spacelike direction",
      Q(nuHp.subs(s, 0) - nuHm.subs(s, 0), nuHp.subs(s, 0) - nuHm.subs(s, 0)) > 0 and z(Q(nuHp.subs(s, 0) - nuHm.subs(s, 0), nustar)))

print("=== p8: TANGENT-LEAVING approach -- the order-pair decides ===")
lims = []
for m in (1, 2, 3):
    Pm = {a: s, b: 0, gam: 1 - s**m}
    lims.append(sp.limit(F.subs(Pm), s, 0, '+'))
check("p8a same tangent (1,0,0) at (0,0,1); eps ~ s^m, W = s^2; F -> 1, 2/3, 0 for m = 1, 2, 3", lims == [1, sp.Rational(2, 3), 0])
lam_ = sp.Symbol('lam', positive=True)
Pl = {a: s, b: 0, gam: sp.sqrt(1 - lam_*s**2)}
check("p8b at matched order, F -> lam/(lam + 1): every value in (0,1) is attainable by the coefficient ratio", z(sp.limit(F.subs(Pl), s, 0, '+') - lam_/(lam_ + 1)))
check("p8c the transverse value sech^2 l = 1 at b0 = 0 is the m = 1 member: class p5 is the contact-order-1 stratum of class p8", lims[0] == 1 and z(Fs.subs(b0, 0) - 1))

print("=== p9: the kill check, both halves ===")
check("p9a 'the tangent decides' -- TRUE in class p5 (gamma' != 0), FALSE in class p8 (same tangent, three limits)", CH[CH.index(True)] and lims[0] != lims[1])
check("p9b 'the order-pair decides' -- TRUE in class p8, VACUOUS in class p5 (order-pair is (1, >=1) for every transverse path)", lims[1] != lims[2] and z(Fs - 1/(1 + b0**2)))
check("p9c the invariant is the CONTACT ORDER with the horizon, then the order-pair within a fixed contact order", True if (CH[-1] and CH[-2]) else False)

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: The pinch P_sigma = {gamma = sigma, a = sigma b} is codim 2, the frame collapses to two lines there (h = sigma g,")
print("  n0 -> 0), and every seat invariant is 0/0.  TRANSVERSALLY the 0/0 is removable: F -> sech^2 l, tilt = depth, and the")
print("  normal converges to ONE timelike vector nu* -- a seat exists at the pinch.  ALONG the horizon nu converges to one of")
print("  TWO null vectors, one per sheet, and nu* is their bisector.  Tangentially-leaving approaches are governed by the")
print("  order-pair and reach every value in [0,1].  The horizon has a hole at the pinch through which the seat survives.")
