# prim_eq1_tolman_from_divider.py  --  EQ-1: thermal equilibrium across the potential ratio, run FORWARD.
#
# OPEN RULING (from earlier sessions): is EQ-1's Tolman condition a principle Will owns or a fact Mercury owns?
# LOCK REGISTER 09-07: neither -- it is a TARGET. This runner tries to derive it from the model.
#
# INPUTS:
#   DYN-1  A(r) = N^2 = W_in/(W_in+W_out) from two spreading resistances (recomputed here, not imported):
#          the lapse N is the ratio of the seat's proper time to the reference clock at infinity, dtau = N dt.   [DERIVED 09-07]
#   SI     hbar and k_B are EXACT by definition (2019). So a temperature is a rate: T = hbar*omega/k_B for a
#          thermal quantum; a thermal spectrum is characterised by ONE frequency scale.                            [GROUND: definitional]
#   COUNT  A rate omega_seat measured against the seat's clock is observed at the reference clock as
#          omega_inf = N * omega_seat: the SAME cycles counted against a clock that runs 1/N times as fast.       [DERIVED: pure counting on DYN-1]
#   EQUIL  Definition: two seats are in equilibrium iff no net exchange occurs in any mode. (This is what the word means.)  [DEFINITION]
#   NYQ    Johnson-Nyquist: a linear element at temperature T emits noise power k_B T per unit bandwidth (white).
#          Noise thermometry realises the kelvin in the 2019 SI.                                                   [GROUND]
#
# BANNED: Tolman-Ehrenfest as an input, any stat-mech equilibrium theorem, any metric, Hawking, Unruh.
# TWO ROUTES: (R1) mode-by-mode occupation balance; (R2) Nyquist power balance in a bandwidth. Must agree.
# KILLS: (K1) if R1 and R2 disagree the counting is wrong. (K2) if the result depends on which mode or which
#        bandwidth, equilibrium is not well-defined and the potential ratio cannot host a thermal state.

import sympy as sp, mpmath as mp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

r, rs, r1, r2 = sp.symbols('r r_s r_1 r_2', positive=True)
w, T, Tinf, df = sp.symbols('omega T T_inf Delta_f', positive=True)

print("=== EQ-1a: the lapse from the potential ratio (DYN-1, recomputed) ===")
RL = sp.integrate(1/r**2, (r, rs, r)); RS = sp.integrate(1/r**2, (r, r, sp.oo))
A = sp.simplify(RL/(RL+RS)); N = sp.sqrt(A)
print("  N^2 =", A, "   N =", N)
check("a1", sp.simplify(A-(1-rs/r))==0, "N^2 = 1 - r_s/r as output")

print("=== EQ-1b: R1 -- mode-by-mode balance ===")
# a mode at the seat has frequency w (seat clock). Seen from the reference clock it is N*w. The seat's thermal
# state assigns occupation n(hbar w / k_B T_seat). The reference reservoir assigns n(hbar (N w) / k_B T_inf).
# EQUIL: no net exchange in ANY mode  =>  the two arguments agree for every w:
Tseat = sp.symbols('T_seat', positive=True)
cond = sp.Eq(w/Tseat, N*w/Tinf)
sol = sp.solve(cond, Tseat)[0]
print("  no net exchange in mode w:  w/T_seat = N w/T_inf   =>  T_seat =", sol)
check("b1", sp.simplify(sol*N - Tinf)==0, "T_seat * N = T_inf : the product of a seat's temperature and its lapse is the reference temperature")
check("b2", w not in sol.free_symbols, "independent of the mode w (K2 passes: equilibrium is well-defined)")

print("=== EQ-1c: R2 -- Nyquist power balance in a band ===")
# seat linear element emits k_B T_seat * df_seat (NYQ). Each quantum arrives at the reference with energy scaled by N (COUNT)
# and the quanta arrive at a rate scaled by N (COUNT): power scales by N^2; the band df_seat maps to N*df_seat.
# The reference linear element emits k_B T_inf * (N df_seat) into that same band. EQUIL: equal.
P_seat_at_inf = (N**2) * Tseat * df           # k_B cancels
P_inf_in_band = Tinf * (N*df)
sol2 = sp.solve(sp.Eq(P_seat_at_inf, P_inf_in_band), Tseat)[0]
print("  N^2 T_seat df = T_inf N df   =>  T_seat =", sol2)
check("c1", sp.simplify(sol2-sol)==0, "R2 agrees with R1 (K1 passes)")
check("c2", df not in sol2.free_symbols, "independent of the bandwidth (K2 passes)")

print("=== EQ-1d: the law and its consequences ===")
Tr = Tinf/N
print("  T(r) =", sp.simplify(Tr))
check("d1", sp.simplify(Tr*sp.sqrt(A) - Tinf)==0, "T(r) sqrt(g_00) = const  -- this IS the Tolman-Ehrenfest condition, DERIVED, not input")
ratio = sp.simplify((Tinf/N.subs(r,r1))/(Tinf/N.subs(r,r2)))
print("  two seats: T(r1)/T(r2) =", ratio)
_tgt = sp.sqrt((1-rs/r2)/(1-rs/r1))
_ok = all(abs(sp.N((ratio-_tgt).subs({rs:1, r1:v1, r2:v2}), 30)) < 1e-25 for v1,v2 in ((sp.Rational(3,2),2),(2,10),(sp.Rational(11,10),1000),(5,sp.Rational(101,100))))
check("d2", _ok, "T1/T2 = N2/N1 (numeric at 4 seat pairs to 1e-25): the deeper seat is hotter, by exactly the clock ratio")
check("d3", sp.limit(Tr, r, rs, '+')==sp.oo, "at the lapse zero (N = 0) N -> 0 and T -> oo: a static thermometer at the horizon reads infinite temperature (no time to count the rate against)")
check("d4", sp.simplify(Tr.subs(r, sp.Rational(3,2)*rs) - sp.sqrt(3)*Tinf)==0, "at the light ring T = sqrt(3) T_inf")
Ain = A.subs(r, rs/2)
check("d5", Ain < 0, f"inside (r = r_s/2): N^2 = {Ain} < 0, no real N: NO static equilibrium exists in the interior. The interior is non-static as a THEOREM of the potential ratio.")
# reciprocity control: swap sigma- (the dual sheet)
Adual = sp.simplify(RS/(RL+RS)); Tdual = Tinf/sp.sqrt(Adual)
print("  dual sheet (involution sigma: eta -> 1/etaped): N^2 =", Adual, "  T_dual(r) =", sp.simplify(Tdual))
check("d6", sp.limit(Tdual, r, sp.oo)==sp.oo and sp.simplify(Tdual.subs(r,rs)-Tinf)==0, "on the dual sheet the divergence moves to infinity and the horizon is at T_inf: the swap exchanges which end is hot (FLAGGED, not claimed)")

print("=== EQ-1e: size of the effect on Earth (for the record; unmeasurable today) ===")
GM=mp.mpf('3.986004418e14'); c=mp.mpf('299792458'); R=mp.mpf('6.371e6')
dTT_per_m = GM/(c**2*R**2)                    # d ln T / dh = -d ln N / dh = -g/c^2
print(f"  d(lnT)/dh = -g/c^2 = {mp.nstr(-dTT_per_m,4)} per metre: a bath 1 m higher is colder by 1 part in {mp.nstr(1/dTT_per_m,3)}")
check("e1", abs(dTT_per_m - mp.mpf('1.09e-16')) < mp.mpf('2e-18'), "1.1e-16 per metre, the same number as the measured clock redshift gradient (optical clocks resolve cm)")

n=sum(CH); print(f"\n=== EQ-1: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("TIER: T*N = T_inf is DERIVED given DYN-1 + SI (hbar, k_B exact) + the definition of equilibrium. Two routes.")
print("RULING: EQ-1 is neither Will's principle nor Mercury's fact. It is a theorem of the potential ratio: temperature is a rate,")
print("        and a rate divides by the lapse. 'Time comes as a pair': (T, N) multiply to the reference node.")
print("GROUND STATUS: the clock-redshift gradient is measured (optical clocks); the thermal statement is one definitional step from it.")
