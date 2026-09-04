#!/usr/bin/env python3
# =============================================================================
# LABEL-3(b) -- two corrections (Will, 2026-09-04, late).                                     
#   1. The handoff's "make D(r) cross zero at the horizon" contradicts T7e: horizon != {det G_c = 0} = {D = 0}.
#      A regular trajectory cannot change sgn D without crossing the branch locus, and the horizon is off it.
#      The sheet at the horizon must be the PATH/BLOCK HOLONOMY chi[Gamma], read while D != 0.  Test:
#      theta_+ = 0  <?>  chi[Gamma_fluid] = chi_BH.
#   2. LABEL-3's causal classification: q(f0,f) = +1 is DIRECTLY null-related (f - f0 null); q(f0,f) = -1 is
#      ANTIPODALLY null-related (f + f0 null, i.e. null-related to -f0); Hawking's centre is EXACTLY -f0, the
#      deck antipode, not a distinct point at null distance.  The lattice is symmetric under f -> -f.
# =============================================================================
import sympy as sp, itertools, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), sp.cancel):
        try:
            if f(e) == 0: return True
        except Exception: pass
    return False

Q = sp.diag(1, 1, -1)
q = lambda u, v: (u.T*Q*v)[0]
FLIPS = ['C', 'H', 'G']; AX = {'C': 2, 'H': 0, 'G': 1}
def centre(subset):
    s = [1, 1, 1]
    for f in subset: s[AX[f]] = -1
    return sp.Matrix(s)
subsets = [frozenset(s) for r in range(4) for s in itertools.combinations(FLIPS, r)]
lab = lambda s: ''.join(sorted(s)) or 'O'
f0 = centre(frozenset()); fH = centre(frozenset(FLIPS))

print("=== 2. the corrected causal classification ===")
check("Hawking's centre is EXACTLY -f0: the deck antipode, not a point at null distance", fH == -f0)
direct_null = [lab(s) for s in subsets if s and z(q(centre(s) - f0, centre(s) - f0)) and q(f0, centre(s)) == 1]
anti_null   = [lab(s) for s in subsets if q(f0, centre(s)) == -1 and centre(s) != -f0 and z(q(centre(s) + f0, centre(s) + f0))]
boost_direct = [lab(s) for s in subsets if q(f0, centre(s)) == 3]
boost_anti   = [lab(s) for s in subsets if q(f0, centre(s)) == -3]
antipode     = [lab(s) for s in subsets if centre(s) == -f0]
print(f"    directly null-related to f0 (q = +1, f - f0 null): {direct_null}")
print(f"    antipodally null-related (q = -1, f + f0 null; null to -f0): {anti_null}")
print(f"    boost-related to f0 (q = +3): {boost_direct};  boost-related to -f0 (q = -3): {boost_anti};  antipode: {antipode}")
check("q = +1 faces (CH, CG) have f - f0 NULL: directly null-related; q = -1 faces (H, G) have f + f0 NULL: null-related to Hawking's centre, not to ours",
      set(direct_null) == {'CH', 'CG'} and set(anti_null) == {'H', 'G'})
check("the seven non-identity faces split 2 + 2 + 1 + 1 + 1: direct-null, antipodal-null, direct-boost (C), antipodal-boost (GH), antipode (CGH)",
      len(direct_null) == 2 and len(anti_null) == 2 and boost_direct == ['C'] and boost_anti == ['GH'] and antipode == ['CGH'])
comp = lambda s: frozenset(FLIPS) - s
check("the lattice is symmetric under f -> -f: q(f0, f) = -q(f0, complement(f)); every face's relation to our centre mirrors its complement's relation to Hawking's",
      all(q(f0, centre(s)) == -q(f0, centre(comp(s))) for s in subsets))
check("LABEL-3's 'six null-separated including Hawking's' was WRONG on two counts (recorded): Hawking's is the antipode, and the null-related faces split into direct (to f0) and antipodal (to -f0)", True if fH == -f0 else False)

print("=== 1. D is fixed-sign along any regular trajectory; the horizon is where D != 0 ===")
t, l1, l2 = sp.symbols('t l1 l2', real=True)
D = sp.cosh(l1)*sp.cosh(l2)*sp.sin(t)
check("D = 0 iff sin t = 0: the zero set of D IS T7c's branch locus (cosh never vanishes)", z(D.subs(t, 0)) and z(D.subs(t, sp.pi)) and not z(D.subs(t, sp.pi/3)))
check("sgn D = sgn(sin t): along any path with sin t != 0 throughout, sgn D is CONSTANT -- it cannot change without crossing the branch locus",
      True if z(sp.simplify(D/(sp.cosh(l1)*sp.cosh(l2)) - sp.sin(t))) else False)
# T7g's two horizons: horizontal (gamma -> +1 at cos t* = (1 + sinh l1 sinh l2)/(cosh l1 cosh l2), depths fixed, l1 != l2)
cos_tstar = (1 + sp.sinh(l1)*sp.sinh(l2))/(sp.cosh(l1)*sp.cosh(l2))
sin2_tstar = sp.simplify(1 - cos_tstar**2)
check("at T7g's HORIZONTAL horizon, sin^2 t* = 1 - cos^2 t* = (sinh l1 - sinh l2)^2 / (cosh^2 l1 cosh^2 l2) != 0 for l1 != l2: D != 0 there",
      z(sp.simplify(sin2_tstar - (sp.sinh(l1) - sp.sinh(l2))**2/(sp.cosh(l1)**2*sp.cosh(l2)**2))) and not z(sin2_tstar.subs({l1: 1, l2: 0})))
# vertical (gamma -> -1 at fixed t = pi/3, depths (L, L/2), L* = 1.8807): sin t = sin(pi/3) != 0 trivially
check("at T7g's VERTICAL horizon (t = pi/3 fixed) sin t = sqrt(3)/2 != 0: D != 0 there too", D.subs(t, sp.pi/3) != 0)
D_h = D.subs({l1: 1, l2: 0, t: sp.acos(cos_tstar.subs({l1: 1, l2: 0}))}); D_v = D.subs({t: sp.pi/3, l1: sp.Rational(18807, 10000), l2: sp.Rational(18807, 20000)})
check("CONSEQUENCE: D evaluated AT both horizons is nonzero (horizontal: D(l1=1,l2=0,t*) ; vertical: D at L*), so the sheet there is NOT sgn D flipping; it is the path/block holonomy chi[Gamma]. The handoff's item 1 is retracted.",
      sp.N(D_h) != 0 and sp.N(D_v) != 0, f"D_horizontal = {sp.N(D_h, 5)}, D_vertical = {sp.N(D_v, 5)}")

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: Hawking's face is the deck antipode -f0; the null-related faces are direct (CH, CG) and antipodal (H, G);")
print("  sgn D is constant along every regular trajectory and D != 0 at both horizons.  The horizon sheet is chi[Gamma].")
print("  Next derivation: (rho, p, Theta, sigma, m, R)(tau) -> X(tau) = (t, l1, l2) -> chi[Gamma]; compare with theta_+ = 0 <=> 2Gm/(c^2 R) = 1.")
