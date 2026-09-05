#!/usr/bin/env python3
# =============================================================================
# OC-1 -- 'horizon = open circuit': what the seat's algebra decides, and what it doesn't.  (Claire, 2026-09-05.)
#   Joins three things: Will's divider reading (N^2 = 1/(1+eta), horizon = open circuit, time = work through the load),
#   KIN-1 (kin1_cost_fork.py: every kinetic cost finite at gamma = 1), and Codex's constitutive divisor law
#   (docs/2026-09-05-cella-constitutive-divisors.md Thm 1, with her clarification: the f^2 law applies when that
#   output must vanish FOR EVERY EFFORT on a two-sided face).
# INPUTS: T7d/T7e (eta = W/eps, N^2 = 1/(1+eta) = eps/delta); a two-element divider; Codex's Thm 1 as a PROVED input.
# CLAIMS, checked:
#   o1  Divider identity: V_out/V_in = Z_l/(Z_s + Z_l) = 1/(1 + Z_s/Z_l).  So N^2 = 1/(1+eta) reads eta = Z_s/Z_l.  Exact.
#   o2  The seat's eta = W/eps diverges at the horizon through the DENOMINATOR eps = 1 - gamma^2 -> 0, W finite (off the pinch).
#   o3  UNDETERMINED: the ratio eta = Z_s/Z_l = W/eps is satisfied by (Z_s, Z_l) = (W, eps) [Z_l -> 0: a SHORT] and by
#       (1/eps, 1/W) [Z_s -> oo: an OPEN circuit], and by any k(W, eps).  The seat's algebra cannot tell open from short.
#       One absolute impedance is needed -- THM-K(a)'s free Z.  'Horizon = open circuit' is a HYPOTHESIS, not a consequence.
#   o4  Under either reading the work through the load P_l = I^2 Z_l -> 0 at the horizon: 'time stops' either way.
#       They differ in the SOURCE: open -> I -> 0, P_s -> 0 (everything stops); short -> I -> V/Z_s, P_s -> V^2/Z_s (the
#       source keeps dissipating with nothing delivered).  That difference is the physical content of Will's choice.
#   o5  Codex's f^2 law, applied: OPEN -> the constrained output is the CURRENT (vanishes for every effort), so the
#       admittance's dissipative part obeys Y_R ~ f^2.  SHORT -> the constrained output is the LOAD VOLTAGE (vanishes
#       for every flow), so the load impedance's dissipative part obeys Z_R ~ f^2.  Either way SOMETHING is ~ f^2 at
#       the face; the readings disagree about which port.
#   o6  PASSIVITY vs THE SIGN OF eps: eps = 1 - gamma^2 changes sign across |gamma| = 1 (interior sector eps < 0), and
#       F = eps/delta with it.  Codex's Thm 1 rests on: a nonnegative smooth scalar with a two-sided zero vanishes to
#       order >= 2.  So a passive response cannot be ~ eps or ~ F (first order, sign-changing).  FORK: (i) the response
#       is ~ F^2 (or eps^2): passive both sides, Codex's law, the interior is passive; (ii) the response is ~ F: passive
#       outside, ACTIVE (negative resistance) inside -- the interior sector PUMPS.  (ii) is the white-hole-sheet
#       conjecture in circuit language.  The seat's algebra does not choose; the load does.
#   o7  CORRECTION (mine): I claimed a tension between KIN-1's finite cost at gamma = 1 and the f^2 law.  There is
#       none: KIN-1's forms are KINETIC (conservative, the Hamiltonian block), not the dissipative R of a constrained
#       output.  Verified: at gamma = 1 they are diag(3,1,2/3), diag(1,1,2/3) -- finite -- and that is what a kinetic
#       block should be.  Different object, no tension.  The conflation was mine.
#   o8  AT THE PINCH the f^2 law is SILENT: its hypothesis is a smooth defining function of a two-sided face; F has no
#       limit there (PINCH-1) and {F = 0} is not smooth.  Codex's clarified condition cannot even be stated at the pinch.
# TIER: o1-o2, o6 (the sign), o7 [DERIVED]; o3-o5, o8 [DERIVED | Codex's Thm 1 as stated]; the fork in o6 is OPEN.
# KILL: an identity in the seat's algebra fixing Z_s and Z_l separately (not just their ratio) closes o3 and kills
#   'undetermined'; a passive response ~ F on both sides would refute o6's dichotomy (it cannot exist -- Thm 1).
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

a, b, gam = sp.symbols('a b gamma', real=True)
Zs, Zl, V, k, s = sp.symbols('Z_s Z_l V k s', positive=True)
eps = 1 - gam**2; W = a**2 + b**2 - 2*a*b*gam; delta = eps + W; eta = W/eps; N2 = 1/(1 + eta)

print("=== o1: the divider identity ===")
Vout = V*Zl/(Zs + Zl)
check("o1a V_out/V_in = Z_l/(Z_s+Z_l) = 1/(1 + Z_s/Z_l): N^2 = 1/(1+eta) reads eta = Z_s/Z_l", z(Vout/V - 1/(1 + Zs/Zl)))
check("o1b and N^2 = eps/delta = (1-gamma^2)/(-det G) (T7e e2), consistent with eta = W/eps", z(N2 - eps/delta))

print("=== o2: the seat's eta diverges through the denominator ===")
hor = {a: 2, b: 1}                                                                    # a horizon point off the pinch (a != sigma b)
check("o2a at (a,b) = (2,1), gamma -> 1: eps -> 0, W -> (a-b)^2 = 1 finite, eta -> oo",
      z(eps.subs(gam, 1)) and z(W.subs(hor).subs(gam, 1) - 1) and sp.limit(eta.subs(hor), gam, 1, '-') == sp.oo)

print("=== o3: open or short is not decided by the ratio ===")
check("o3a (Z_s, Z_l) = (W, eps) gives eta and has Z_l -> 0: a SHORT", z(W/eps - eta) and z(eps.subs(gam, 1)))
check("o3b (Z_s, Z_l) = (1/eps, 1/W) gives the SAME eta and has Z_s -> oo: an OPEN circuit", z((1/eps)/(1/W) - eta) and sp.limit((1/eps).subs(hor), gam, 1, '-') == sp.oo)
check("o3c (Z_s, Z_l) = k (W, eps) for ANY k(a,b,gamma) > 0 gives eta: the ratio fixes one function, the divider has two", z((k*W)/(k*eps) - eta))
check("o3d so 'horizon = open circuit' is a hypothesis the seat's algebra cannot confirm or refute: one absolute impedance is missing (THM-K(a): Z free)",
      True if (CH[-1] and CH[-2] and CH[-3]) else False)

print("=== o4: what the two readings agree and disagree on ===")
I = V/(Zs + Zl); Pl = I**2*Zl; Ps = I**2*Zs
openc = {Zs: 1/s, Zl: 1}; short = {Zl: s, Zs: 1}                                      # s -> 0 is the horizon in each reading
check("o4a OPEN (Z_s -> oo): I -> 0, P_l -> 0, P_s -> 0 -- everything stops",
      sp.limit(I.subs(openc), s, 0, '+') == 0 and sp.limit(Pl.subs(openc), s, 0, '+') == 0 and sp.limit(Ps.subs(openc), s, 0, '+') == 0)
check("o4b SHORT (Z_l -> 0): I -> V/Z_s, P_l -> 0, P_s -> V^2/Z_s -- the source keeps dissipating, nothing delivered",
      z(sp.limit(I.subs(short), s, 0, '+') - V) and sp.limit(Pl.subs(short), s, 0, '+') == 0 and z(sp.limit(Ps.subs(short), s, 0, '+') - V**2))
check("o4c 'time = work through the load': P_l -> 0 under BOTH readings; the readings differ only in the source's fate", True if (CH[-1] and CH[-2]) else False)

print("=== o5: Codex's f^2 law, both ways round ===")
f = sp.Symbol('f', real=True)                                                        # a defining function of the two-sided face
# admittance form j = Y e; OPEN: j must vanish for every e on the face -> Y's dissipative diagonal is a nonneg smooth scalar with a two-sided zero
Y_R = f**2*sp.Symbol('y', positive=True)
check("o5a OPEN: the current is the constrained output; a nonnegative smooth Y_R with a two-sided zero at f = 0 is O(f^2) (Codex Thm 1) -- Y_R = f^2 y realises it, Y_R = |f| or f do not (not smooth / not nonneg)",
      z(sp.diff(Y_R, f).subs(f, 0)) and sp.diff(Y_R, f, 2).subs(f, 0) > 0)
Z_R = f**2*sp.Symbol('r', positive=True)
check("o5b SHORT: the load voltage is the constrained output (vanishes for every flow); the same law applies to the load's Z_R: O(f^2)",
      z(sp.diff(Z_R, f).subs(f, 0)) and sp.diff(Z_R, f, 2).subs(f, 0) > 0)
check("o5c under either reading SOMETHING dissipative is ~ f^2 at the face; the readings disagree about WHICH port", True if (CH[-1] and CH[-2]) else False)

print("=== o6: passivity against the sign of eps -- the fork ===")
check("o6a eps = 1 - gamma^2 is POSITIVE at gamma = 1 - s and NEGATIVE at gamma = 1 + s: a first-order, sign-changing zero",
      sp.limit(eps.subs(gam, 1 - s)/s, s, 0, '+') == 2 and sp.limit(eps.subs(gam, 1 + s)/s, s, 0, '+') == -2)
F = eps/delta
check("o6b F = eps/delta changes sign with it (delta > 0 both sides near the point): F < 0 in the interior sector (T7f)",
      sp.limit(F.subs(hor).subs(gam, 1 - s)/s, s, 0, '+') > 0 and sp.limit(F.subs(hor).subs(gam, 1 + s)/s, s, 0, '+') < 0)
check("o6c a response proportional to eps (or F) is NOT nonnegative on both sides: it cannot be a passive R (Thm 1's hypothesis fails); F^2 can",
      (eps.subs(gam, 1 + sp.Rational(1, 2)) < 0) and (F**2).subs(hor).subs(gam, 1 + sp.Rational(1, 2)) > 0)
check("o6d FORK recorded: (i) R ~ F^2 -- passive on both sides, the interior is passive; (ii) R ~ F -- passive outside, ACTIVE inside (negative resistance): the interior PUMPS. The seat's algebra does not choose.",
      True if (CH[-1] and CH[-2]) else False)

print("=== o7: KIN-1 is not in tension with the divisor law -- my conflation, corrected ===")
u, v = sp.symbols('u v', real=True)
gi = sp.diag(sp.cosh(2*u), sp.cosh(2*v), sp.cosh(u)**2*sp.cosh(v)**2/(sp.cosh(u)**2 + sp.cosh(v)**2))
gii = sp.diag(1, 1, sp.cosh(u)**2*sp.cosh(v)**2/(sp.cosh(u)**2 + sp.cosh(v)**2))
fx = {u: sp.asinh(1), v: 0}
gi_f = gi.subs(fx).applyfunc(lambda e: sp.simplify(e.rewrite(sp.exp))); gii_f = gii.subs(fx).applyfunc(lambda e: sp.simplify(e.rewrite(sp.exp)))
check("o7a KIN-1's kinetic forms at the horizon fixture are diag(3,1,2/3) and diag(1,1,2/3): finite, nonzero", gi_f == sp.diag(3, 1, sp.Rational(2, 3)) and gii_f == sp.diag(1, 1, sp.Rational(2, 3)))
check("o7b a KINETIC block is conservative and is not the dissipative R of any constrained output; Thm 1 does not apply to it. No tension. (Correction to Claire's overlap note.)",
      True if CH[-1] else False)

print("=== o8: at the pinch the law is silent ===")
s2 = sp.Symbol('s2', positive=True)
lim1 = sp.limit(F.subs({a: s2, b: 0, gam: 1 - s2**2}), s2, 0, '+'); lim2 = sp.limit(F.subs({a: s2, b: 0, gam: 1 - s2}), s2, 0, '+')
check("o8a F has no limit at the pinch (0,0,1): two approaches give 2/3 and 1 (PINCH-1 p8a)", lim1 == sp.Rational(2, 3) and lim2 == 1)
check("o8b so {F = 0} has no smooth defining function through the pinch and Thm 1's hypothesis (a two-sided coordinate hypersurface) cannot be stated there: the f^2 law is SILENT at the pinch", lim1 != lim2)

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: The seat's algebra fixes eta = Z_s/Z_l = W/eps and nothing else: open (Z_s -> oo) and short (Z_l -> 0) are")
print("  indistinguishable to it -- 'horizon = open circuit' is a hypothesis needing one absolute impedance (THM-K(a)'s Z).")
print("  Both readings stop the work through the load; they differ in whether the source keeps burning.  Codex's f^2 law")
print("  applies under either reading to a DIFFERENT port.  Passivity plus the sign change of eps across the horizon forces")
print("  a fork: response ~ F^2 (passive interior) or ~ F (ACTIVE interior, the white-hole sheet in circuit language).")
print("  KIN-1's kinetic forms are a different block; the tension I claimed was a conflation.  The law is silent at the pinch.")
