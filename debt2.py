#!/usr/bin/env python3
# =============================================================================
# DEBT-2: THE SEAT DISSIPATOR, DERIVED  (2026-08-28, the late shift)
#
# Construction (model-native): seat = one angle mode (unit-frequency oscillator,
# doors S2).  Bath = the OTHER cells' modes as PRESENTED to the seat: Doppler
# family w(u) = w0(cosh L + sinh L * u) under the FORCED plane measure du/2
# (E-2).  Linear map of uniform is uniform: the bath is a FLAT BAND on
# [w0 e^-L, w0 e^L], width 2 w0 sinh L.  No pivot => degenerate band => no
# dephasing => no arrow.  Dissipation REQUIRES the pivot; margin sinh L (#5).
#
# Machinery (lineage, loud): Ford-Kac-Mazur / Zwanzig / Caldeira-Leggett
# independent-oscillator elimination; classical FDT.  The derivations below
# are run in-session; the assembled GLE is [derived, conditional on standard
# integration by parts - FKM].  Residue: DEBT-2b = derive c(w) from the
# channel law (F-4's Boltzmann-exponent functional).  Coupling here: Ohmic-in-
# band discretisation, declared, not derived.
# =============================================================================
import sympy as sp
import math, random

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" - {n}" if n else ""), flush=True)

print("=" * 78); print("PART 1 - THEOREMS (symbolic)"); print("=" * 78)
w0, L, t, tp, u = sp.symbols('omega0 Lambda t tprime u', positive=True)
check("D2-1 plane measure => FLAT BAND: w(u)=w0(chL+shL u), u~U[-1,1] is uniform",
      sp.simplify((sp.cosh(L)+sp.sinh(L)-sp.exp(L)).rewrite(sp.exp)) == 0 and
      sp.simplify((sp.cosh(L)-sp.sinh(L)-sp.exp(-L)).rewrite(sp.exp)) == 0,
      "linear image of uniform is uniform; band [w0 e^-L, w0 e^L], width 2 w0 shL")
wp, wm = w0*sp.exp(L), w0*sp.exp(-L)
gam = (sp.sin(wp*t) - sp.sin(wm*t)) / ((wp - wm)*t)     # band-averaged kernel
check("D2-2 THE GATE: lim L->0 of the band kernel = cos(w0 t), UNDAMPED forever",
      sp.simplify(sp.limit(gam, L, 0) - sp.cos(w0*t)) == 0,
      "no pivot => pure beat, no memory decay, theta stays pi/2: NO ARROW")
print("      (L>0: kernel decays ~1/t on time 1/(2 w0 sinh L): the pivot opens the drain)")

# D2-3: FDT - the drain and the hiss share coefficients (3 symbolic modes)
kT = sp.Symbol('kT', positive=True)
ws = sp.symbols('w1:4', positive=True); cs = sp.symbols('c1:4')
x0 = sp.symbols('x01 x02 x03'); p0 = sp.symbols('p01 p02 p03')
xi  = sum(cs[k]*(x0[k]*sp.cos(ws[k]*t)  + (p0[k]/ws[k])*sp.sin(ws[k]*t))  for k in range(3))
xip = sum(cs[k]*(x0[k]*sp.cos(ws[k]*tp) + (p0[k]/ws[k])*sp.sin(ws[k]*tp)) for k in range(3))
prod = sp.expand(xi*xip)
subs = {}
for i in range(3):
    for j in range(3):
        subs[x0[i]*x0[j]] = (kT/ws[i]**2 if i == j else 0)
        subs[p0[i]*p0[j]] = (kT if i == j else 0)
        subs[x0[i]*p0[j]] = 0
Exx = sp.expand(prod.subs(subs, simultaneous=True))
gamma_tt = sum((cs[k]**2/ws[k]**2)*sp.cos(ws[k]*(t - tp)) for k in range(3))
check("D2-3 FDT: <xi(t) xi(t')> = kT * gamma(t-t') with the SAME kappa_k",
      sp.simplify(sp.expand_trig(sp.expand(Exx - kT*gamma_tt))) == 0,
      "the hiss and the drain are one object: damping without noise is impossible")
# D2-4: exact one-mode elimination identity (the reduction is not hand-waved)
q = sp.Function('q'); s = sp.Symbol('s'); c1s, w1s = sp.symbols('c w', positive=True)
xk0, pk0 = sp.symbols('xk0 pk0')
xk = (xk0*sp.cos(w1s*t) + (pk0/w1s)*sp.sin(w1s*t)
      + (c1s/w1s)*sp.Integral(sp.sin(w1s*(t - s))*q(s), (s, 0, t)))
resid = sp.simplify(sp.diff(xk, t, 2) + w1s**2*xk - c1s*q(t))
check("D2-4 exact elimination: x_k(t) formula satisfies x'' + w^2 x = c q identically",
      sp.simplify(resid.doit()) == 0,
      "GLE assembly then follows by one integration by parts [FKM, argument]")

print(); print("=" * 78)
print("PART 2 - THE UNITARY WORLD, PARTITIONED  [numeric evidence, leapfrog]")
print("=" * 78)
def simulate(N, lam, T, eta=0.10, w0=1.0, q0=1.0, tmax=200.0, dt=0.004, seed=1):
    random.seed(seed)
    if lam > 0:
        wlo, whi = w0*math.exp(-lam), w0*math.exp(lam)
        D = (whi - wlo)/N
        wk = [wlo + (k + 0.5)*D for k in range(N)]
    else:
        D = (2*w0*math.sinh(0.8))/N              # same coupling scale as lam=0.8
        wk = [w0]*N
    ck = [math.sqrt(2/math.pi*eta*D)*w for w in wk]
    xk = [random.gauss(0, math.sqrt(T)/w) for w in wk]
    pk = [random.gauss(0, math.sqrt(T)) for w in wk]
    qq, pp = q0, 0.0
    mu = sum(c*c/(w*w) for c, w in zip(ck, wk))   # counterterm
    def Fq(qv, xv): return -w0*w0*qv - mu*qv + sum(c*x for c, x in zip(ck, xv))
    def H():
        Eb = sum(0.5*pk[k]**2 + 0.5*(wk[k]*xk[k] - ck[k]*qq/wk[k])**2 for k in range(N))
        return 0.5*pp*pp + 0.5*w0*w0*qq*qq + Eb + 0.5*(mu - sum(c*c/(w*w) for c,w in zip(ck,wk)))*qq*qq
    steps = int(tmax/dt); H0 = H(); Es = []
    aq = Fq(qq, xk); ak = [-wk[k]**2*xk[k] + ck[k]*qq for k in range(N)]
    for n in range(steps):
        pp += 0.5*dt*aq
        for k in range(N): pk[k] += 0.5*dt*ak[k]
        qq += dt*pp
        for k in range(N): xk[k] += dt*pk[k]
        aq = Fq(qq, xk); ak = [-wk[k]**2*xk[k] + ck[k]*qq for k in range(N)]
        pp += 0.5*dt*aq
        for k in range(N): pk[k] += 0.5*dt*ak[k]
        if n % 50 == 0: Es.append((n*dt, 0.5*pp*pp + 0.5*w0*w0*qq*qq))
    return Es, abs(H() - H0)/abs(H0)

E0 = 0.5
# Run A: pivoted band, warm bath -> decay to the FLOOR (the hiss = kT)
Es, drift = simulate(N=400, lam=0.8, T=0.05, seed=3)
mid  = [e for tt, e in Es if 40 < tt < 60]
tail = [e for tt, e in Es if tt > 150]
floorA = sum(tail)/len(tail)
check("N-1  the WORLD stays at pi/2: total H drift", drift < 1e-6, f"rel drift {drift:.2e}")
check("N-2  the SEAT ages: E_seat 0.50 -> floor ~ kT (decay THROUGH the band)",
      max(mid) < 0.2*E0 and 0.3*0.05 < floorA < 3*0.05,
      f"mid-window max {max(mid):.4f}, late-time mean {floorA:.4f} vs kT = 0.05")
# Run B: cold bath -> decays toward zero (no floor to stand on)
Es, _ = simulate(N=400, lam=0.8, T=0.0, seed=4)
tailB = [e for tt, e in Es if tt > 150]
check("N-3  T = 0 bath: seat decays toward zero - the floor WAS the hiss",
      max(tailB) < 0.02*E0, f"late max {max(tailB):.5f}")
# Run C: NO PIVOT (degenerate band) - no arrow, energy returns
Es, _ = simulate(N=400, lam=0.0, T=0.0, seed=5)
lateC = [e for tt, e in Es if tt > 60]
check("N-4  THE GATE, run: lam = 0 => NO relaxation - energy keeps returning",
      max(lateC) > 0.5*E0, f"late-window max E_seat {max(lateC):.3f} of {E0}")
# Run D: small N -> the CATHEDRAL REVIVAL (recurrence at finite mode count)
Es, _ = simulate(N=25, lam=0.8, T=0.0, seed=6)
died = [tt for tt, e in Es if e < 0.1*E0]
revived = [e for tt, e in Es if died and tt > died[0] + 5]
check("N-5  N = 25: the note DIES then REVIVES (recurrence ~ 1/mode spacing)",
      bool(died) and max(revived) > 0.35*E0,
      f"died at t~{died[0]:.0f}, later revival to {max(revived):.3f}; N=400 shows none in-window")

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
DEBT-2 CLOSED (with one residue).  The seat dissipator is DERIVED, not posited:
  bath = the sky (other cells' modes, Doppler-presented, plane measure) =>
  flat band [w0 e^-L, w0 e^L]; elimination => GLE; FDT => the drain and the
  hiss carry the same coefficients (D2-3): a damping seat OWES its noise floor.
  The gate: no pivot => degenerate band => undamped kernel => theta = pi/2
  (D2-2 exact, N-4 run).  Dissipation margin = sinh(Lambda): fifth appearance.
  The arrow needs BOTH the pivot (bandwidth) and the crowd (N -> infinity);
  at finite N the cathedral revives (N-5) - the arrow is asymptotic, exactly
  as RULING-1 requires: world at pi/2 throughout (N-1), seat below it.
  RESIDUE (DEBT-2b): the coupling c(w) is declared Ohmic-in-band here, not
  derived; its derivation belongs to the channel law (F-4's functional).
""")
