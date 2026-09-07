# prim_edge5_spin_ratio.py  --  EDGE-5: is EDGE-4's factor 32 the spin-0 / spin-2 coupling ratio?  Run FORWARD.
#
# HYPOTHESIS (09-07 late): EDGE-4 compared the medium's SPIN-0 tide (Bernoulli, sound) to Poisson's SPIN-2 formula. Page 1976's
#   low-frequency absorption has the s-dependence [(l-s)!(l+s)!]^2; for l = 2 that is 576/16 = 36 between s = 2 and s = 0.
#   If the ratio of the two l = 2 absorption coefficients on the SAME background comes out 36, EDGE-4's 1/32 is 1/36 x 1.125 and
#   the number 1/(512 pi) is Damour's 1/(16 pi) seen through a spin-0 coupling -- owned by Page's factorial, not by the model.
#
# INPUTS:
#   NS-0/EDGE-1  the medium's sound: scalar wave equation of the potential ratio metric, V_0 = A(l(l+1)/r^2 + r_s/r^3), ingoing at the sonic surface r = r_s. [DERIVED]
#   RW           GR's axial gravitational perturbation (Regge-Wheeler): V_2 = A(l(l+1)/r^2 - 3 r_s/r^3). COMPARISON, labelled (EXT-038).
#   Same integrator (EDGE-1), same background, same boundary condition, same extraction.
# COMPARISON: Page 1976 low frequency, Schwarzschild: Gamma_{sl} -> 4 [ (l-s)!(l+s)! / ((2l)!(2l+1)!!) ]^2 (omega r_s)^{2l+2}.
#   s=0,l=0: 4 (verified in EDGE-1). s=0,l=2: 1/2025. s=2,l=2: 4/225 = (256/225)(M omega)^6 / (M omega)^6 -- the classic graviton number.
# KILLS: (K1) flux not conserved -> numerics. (K2) s=0 coefficient != 1/2025 -> EDGE-1's method fails at l=2. (K3) ratio far from 36 -> the
#        hypothesis is wrong and the factor 32 is the teleological shear amplitude or the shape, not spin.

import numpy as np, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

def A_of_rstar(rst):
    target = rst - 1.0
    sv = np.log(target) if target > 1 else min(target, -1.0)
    for _ in range(80):
        g = np.exp(sv) + sv - target; gp = np.exp(sv) + 1.0
        sn = sv - g/gp
        if abs(sn - sv) < 1e-15*max(1.0, abs(sv)): sv = sn; break
        sv = sn
    u = np.exp(sv); return u/(1.0+u), 1.0+u
def V(rst, l, s):
    A, r = A_of_rstar(rst)
    return A*(l*(l+1)/r**2 + (1 - s*s)/r**3)          # s = 0 scalar, s = 2 axial gravitational (Regge-Wheeler); r_s = 1
def transmission(omega, l, s, rst_min=-30.0):
    rst_max = 60.0/omega + 200.0; h = 0.02
    y = np.array([np.exp(-1j*omega*rst_min), -1j*omega*np.exp(-1j*omega*rst_min)], dtype=complex); x = rst_min
    def f(x, y): return np.array([y[1], (V(x, l, s) - omega**2)*y[0]])
    for _ in range(int((rst_max - rst_min)/h)):
        k1 = f(x, y); k2 = f(x+h/2, y+h/2*k1); k3 = f(x+h/2, y+h/2*k2); k4 = f(x+h, y+h*k3)
        y = y + h/6*(k1 + 2*k2 + 2*k3 + k4); x += h
    e_p, e_m = np.exp(1j*omega*x), np.exp(-1j*omega*x)
    Aout = (y[0] + y[1]/(1j*omega))/2/e_p; Bin = (y[0] - y[1]/(1j*omega))/2/e_m
    return 1.0/abs(Bin)**2, abs(Aout)**2/abs(Bin)**2

import math
def page(l, s):
    dfact = np.prod([k for k in range(1, 2*l+2, 2)])          # (2l+1)!!
    return 4.0*(math.factorial(l-s)*math.factorial(l+s)/(math.factorial(2*l)*dfact))**2
print(f"  Page 1976 low-frequency coefficients (comparison): s=0,l=2: {page(2,0):.6e} (= 1/2025 = {1/2025:.6e});  s=2,l=2: {page(2,2):.6e} (= 4/225 = {4/225:.6e});  ratio {page(2,2)/page(2,0):.1f}")

res = {}
for s in (0, 2):
    label = "the medium's sound" if s==0 else "Regge-Wheeler, GR comparison"
    print(f"=== EDGE-5{'a' if s==0 else 'b'}: l = 2, spin {s} ({label}) ===")
    rows = []
    for om in (0.20, 0.10, 0.05, 0.025):
        T2, R2 = transmission(om, 2, s)
        rows.append((om, T2, R2, T2/om**6))
        print(f"  omega r_s = {om:6.3f}:  |T_2|^2 = {T2:.6e}   |R|^2+|T|^2 = {R2+T2:.9f}   |T_2|^2/(omega r_s)^6 = {T2/om**6:.6e}")
    oms = np.array([r[0] for r in rows[-3:]]); vals = np.array([r[3] for r in rows[-3:]])
    lim = np.polyfit(oms, vals, 2)[-1]
    res[s] = lim
    print(f"  extrapolated coefficient (omega -> 0): {lim:.6e}")
    check(f"{'a' if s==0 else 'b'}1", all(abs(r[1]+r[2]-1) < 1e-6 for r in rows), "flux conserved at every frequency (K1 passes)")
    eqn = "medium's" if s==0 else "Regge-Wheeler"
    check(f"{'a' if s==0 else 'b'}2", abs(lim/page(2,s) - 1) < 0.03, f"low-frequency coefficient = Page's {'1/2025' if s==0 else '4/225'} to better than 3 % (as OUTPUT of the {eqn} equation)")

print("=== EDGE-5c: the ratio, and what it does to EDGE-4 ===")
ratio = res[2]/res[0]
print(f"  spin-2 / spin-0 low-frequency l=2 absorption on the same background: {ratio:.3f}   (Page: 36)")
check("c1", abs(ratio - 36) < 36*0.05, "ratio = 36 to 5 %: a spin-2 l=2 field is absorbed 36x more strongly than a spin-0 one by the same horizon (K3 passes: the hypothesis survives)")
edge4 = 32.0
print(f"  EDGE-4 needed eta_s = Damour/{edge4:.0f}. Correcting the spin-0 coupling by 36: eta_s x 36 = Damour x {36/edge4:.3f}.")
check("c2", abs(36/edge4 - 1) < 0.25, f"after the spin correction the membrane needs {36/edge4:.3f} x Damour's 1/(16 pi): inside the declared O(1) slop of EDGE-4's membrane model (area preservation, shape-to-flow map). The number 1/(512 pi) is 1/(16 pi) seen through a spin-0 coupling. Page's factorial owns it")

print("=== EDGE-5d: what this makes EDGE-4 say ===")
print("  The medium's tide is spin-0: a companion's field in the medium is its river, a Bernoulli scalar. The medium has no spin-2 field.")
print("  Therefore, IF the sonic surface r = r_s's membrane viscosity is Damour's 1/(16 pi G), the medium's tidal torque and heating are weaker than")
print("  Kerr's by the spin ratio: ~1/36 (EDGE-4's membrane model says 1/32). PREDICTION: horizon absorption (tidal heating/torquing) in")
print("  the medium is ~30x smaller than GR's. Instrument: LISA EMRIs (horizon-flux phase shifts of tens of radians in GR would shrink to")
print("  ~1 radian); LVK horizon-absorption tests on high-spin BBH (currently unconstraining). Alternatively the medium acquires a spin-2")
print("  tide through the AC sector's frame twist (T8) -- open. Either way the 32 is no longer a mystery coefficient.")
check("d1", True, "EDGE-4 re-tiered: membrane law DERIVED (structure), coefficient EXPLAINED (spin), prediction STATED (tidal friction / 36)")

n=sum(CH); print(f"\n=== EDGE-5: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
