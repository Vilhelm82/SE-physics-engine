#!/usr/bin/env python3
# =============================================================================
# FIGURES for PAPER-seated-root  (2026-08-28)
# Fig 1  the elliptope: state space, Cayley nodes, collision lines, strata
# Fig 2  the trichotomy: one rotor, three relationships (free/floored/pinned)
# Fig 3  the Kahler dial: theta interpolating Door C and Door R
# All three regenerate from the SAME mathematics the suites verify.
# =============================================================================
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

plt.rcParams.update({'font.size': 9, 'axes.linewidth': 0.8,
                     'figure.dpi': 200, 'savefig.bbox': 'tight'})
INK, ACC, WARN = '#1a1a1a', '#1f5c8b', '#a33b20'

# ---------------------------------------------------------------- FIG 1 -----
def fig1():
    fig = plt.figure(figsize=(9.6, 4.6))
    ax = fig.add_subplot(121, projection='3d')
    # boundary = rank-2 Grams: three unit vectors IN A PLANE at angles 0, a, b
    n = 260
    a = np.linspace(0, 2*np.pi, n); b = np.linspace(0, 2*np.pi, n)
    A, B = np.meshgrid(a, b)
    X, Y, Z = np.cos(A), np.cos(B), np.cos(A - B)
    ax.plot_surface(X, Y, Z, rstride=5, cstride=5, color=ACC, alpha=0.20,
                    linewidth=0, antialiased=True, shade=True)
    t = np.linspace(-1, 1, 60)
    for (lx, ly, lz) in [( np.ones_like(t),  t,  t), (-np.ones_like(t),  t, -t),
                         ( t,  np.ones_like(t),  t), ( t, -np.ones_like(t), -t),
                         ( t,  t,  np.ones_like(t)), ( t, -t, -np.ones_like(t))]:
        ax.plot(lx, ly, lz, color=WARN, lw=2.0, zorder=5)
    for (x, y, z) in [(1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1)]:
        ax.scatter([x],[y],[z], s=46, color=WARN, depthshade=False, zorder=10)
    ax.set_xlabel(r'$\gamma_{12}$', labelpad=-8); ax.set_ylabel(r'$\gamma_{13}$', labelpad=-8)
    ax.set_zlabel(r'$\gamma_{23}$', labelpad=-8)
    ax.set_xticks([-1,0,1]); ax.set_yticks([-1,0,1]); ax.set_zticks([-1,0,1])
    ax.tick_params(pad=-3); ax.view_init(elev=22, azim=42)
    ax.set_title('(a) $\\Delta=0$: the Cayley cubic', pad=0)

    ax2 = fig.add_subplot(122)
    g = np.linspace(-1.02, 1.02, 900); P, Q = np.meshgrid(g, g)
    shades = ['#cfe0ec', '#9dc0d8', '#5d94b8', '#1f5c8b']
    for s_, col in zip((0.0, 0.6, 0.9, 0.99), shades):
        Dl = 1 - P**2 - Q**2 - s_**2 + 2*P*Q*s_
        ax2.contourf(P, Q, Dl, levels=[0, 10], colors=[col], alpha=0.75)
        ax2.contour(P, Q, Dl, levels=[0], colors=[INK], linewidths=[0.6])
    ax2.scatter([1,-1], [1,-1], s=40, color=WARN, zorder=8)
    ax2.annotate('node', xy=(1,1), xytext=(0.30, 0.86), fontsize=8, color=WARN,
                 arrowprops=dict(arrowstyle='->', color=WARN, lw=0.9))
    ax2.set_xlabel(r'$\gamma_{12}$'); ax2.set_ylabel(r'$\gamma_{13}$')
    ax2.set_aspect('equal'); ax2.set_xlim(-1.1, 1.1); ax2.set_ylim(-1.1, 1.1)
    ax2.set_title(r'(b) slices $\gamma_{23}=0,\,0.6,\,0.9,\,0.99$', pad=6)
    fig.subplots_adjust(bottom=0.30, wspace=0.24, left=0.055, right=0.97)
    fig.text(0.5, 0.015,
        'Figure 1. The state space is the $3\\times3$ elliptope. Its boundary $\\Delta=0$ is the Cayley cubic:\n'
        'smooth points form the rank-2 coplanar stratum; the four red nodes are the rank-1 total collisions,\n'
        'which are also the cube\'s body diagonals, permuted as $S_4$ by the click group. Red lines are the six\n'
        'collision loci $\\gamma_{ij}=\\pm1$ (six of the surface\'s nine lines; the other three lie at infinity).\n'
        'In (b) the admissible region contracts as $\\gamma_{23}\\to1$, closing onto a node.',
        ha='center', va='bottom', fontsize=8)
    fig.savefig('fig1_elliptope.png'); plt.close(fig); print("fig1 ok")
fig1()

# ---------------------------------------------------------------- FIG 2 -----
def fig2():
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 4.2))
    lam = np.linspace(0, 3.0, 500)

    ax = axes[0]                      # c seat: FREE - the product is 1 on-state
    th = np.linspace(0, np.pi, 400)
    for L, al in ((0.4, .35), (0.9, .6), (1.6, 1.0)):
        ax.plot(th, np.cosh(L) + np.sinh(L)*np.cos(th), color=ACC, alpha=al,
                lw=1.3, label=f'$\\lambda={L}$' if L == 1.6 else None)
    ax.axhline(0, color=INK, lw=.7, ls=':')
    ax.set_xlabel(r'$\theta$'); ax.set_ylabel(r'$s(\theta)$')
    ax.set_xticks([0, np.pi/2, np.pi]); ax.set_xticklabels(['0', r'$\pi/2$', r'$\pi$'])
    ax.set_title(r'$c$ seat: FREE', color=ACC)
    ax.text(.05, .06, r'$s(0)\,s(\pi)=1$;  null rays exist ON the state',
            transform=ax.transAxes, fontsize=7.5)

    ax = axes[1]                      # hbar seat: FLOORED - e^{-2 lam}, never 0
    ax.semilogy(lam, np.exp(-2*lam), color=ACC, lw=1.5)
    ax.axhline(0, color=WARN, lw=1.2)
    ax.set_xlabel(r'$\lambda$  (free)'); ax.set_ylabel(r'$\min_\phi q\,/\,A$')
    ax.set_title(r'$\hbar$ seat: FLOORED', color=ACC)
    ax.text(.30, .70, 'zero only at\n$\\lambda=\\infty$:\norbit separation',
            transform=ax.transAxes, fontsize=7.5, color=WARN)

    ax = axes[2]                      # G seat: PINNED - wall gets an address
    r = np.linspace(1.0001, 6, 800)
    ax.plot(r, 1 - 1/r, color=ACC, lw=1.5, label=r'$N^2=\mathrm{sech}^2\lambda(r)$')
    ax.plot(r, 1 - 1/np.sqrt(r), color=INK, lw=1.0, ls='--', label=r'escaping pole $1-v$')
    ax.plot(r, 1 + 1/np.sqrt(r), color=INK, lw=1.0, ls=':', label=r'bound pole $1+v$')
    ax.axvline(1, color=WARN, lw=1.6)
    ax.set_xlabel(r'$r/r_s$'); ax.set_ylim(-0.05, 2.15); ax.set_xlim(0.8, 6)
    ax.set_title(r'$G$ seat: HORIZON-VANISHING', color=ACC)
    ax.legend(fontsize=6.6, loc='center right', frameon=False)
    ax.text(1.12, 0.12, r'wall at $r=r_s$', fontsize=7.5, color=WARN)
    fig.text(0.5, 0.015,
        'Figure 2. One rotor, three relationships between pivot parameter and state. The vanishing locus\n'
        'sits ON the state ($c$: free), at $\\lambda=\\infty$ with $\\lambda$ free ($\\hbar$: floored, unreachable),\n'
        'or at $\\lambda=\\infty$ with $\\lambda$ PINNED to position by $\\tanh\\lambda(r)=\\sqrt{r_s/r}$ ($G$: the wall\n'
        'acquires the address $r=r_s$). The three characters were frozen 2026-08-17, before any was computed.',
        ha='center', va='bottom', fontsize=8)
    fig.subplots_adjust(bottom=0.34, wspace=0.34)
    fig.savefig('fig2_trichotomy.png'); plt.close(fig); print("fig2 ok")
fig2()

# ---------------------------------------------------------------- FIG 3 -----
def fig3():
    import math
    def dot(u,v): return sum(a*b for a,b in zip(u,v))
    def sub(u,v): return [a-b for a,b in zip(u,v)]
    def add(u,v): return [a+b for a,b in zip(u,v)]
    def scl(s,u): return [s*a for a in u]
    def crs(u,v): return [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
    def nrm(u):
        n = math.sqrt(dot(u,u)); return [a/n for a in u]
    def gam(A): return (dot(A[0],A[1]), dot(A[0],A[2]), dot(A[1],A[2]))
    def Dv(g): return 1-g[0]**2-g[1]**2-g[2]**2+2*g[0]*g[1]*g[2]
    def grads(A):
        g = gam(A); D = Dv(g)
        c = (-(-2*g[0]+2*g[1]*g[2])/(2*D), -(-2*g[1]+2*g[0]*g[2])/(2*D),
             -(-2*g[2]+2*g[0]*g[1])/(2*D))
        return [add(scl(c[0],A[1]), scl(c[1],A[2])),
                add(scl(c[0],A[0]), scl(c[2],A[2])),
                add(scl(c[1],A[0]), scl(c[2],A[1]))]
    def rhs(A, th):
        Gr = grads(A); out=[]
        for i in range(3):
            Pg = sub(Gr[i], scl(dot(Gr[i],A[i]), A[i]))
            out.append(sub(scl(math.sin(th), crs(A[i],Gr[i])), scl(math.cos(th), Pg)))
        return out
    def run(th, N=6000, h=0.002):
        A = [[1.0,0,0],[0.55, math.sqrt(1-.3025),0],
             [-0.30,0.35,math.sqrt(1-.09-.1225)]]
        ts, Vs, g12 = [], [], []
        for k in range(N):
            k1=rhs(A,th); A2=[nrm(add(A[i],scl(h/2,k1[i]))) for i in range(3)]
            k2=rhs(A2,th); A3=[nrm(add(A[i],scl(h/2,k2[i]))) for i in range(3)]
            k3=rhs(A3,th); A4=[nrm(add(A[i],scl(h,k3[i]))) for i in range(3)]
            k4=rhs(A4,th)
            A=[nrm(add(A[i],scl(h/6,add(add(k1[i],scl(2,k2[i])),
                                        add(scl(2,k3[i]),k4[i]))))) for i in range(3)]
            if k % 6 == 0:
                ts.append(k*h); Vs.append(-0.5*math.log(Dv(gam(A)))); g12.append(gam(A)[0])
        return ts, Vs, g12

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
    cols = [(0.0,'#7a1f1f'), (math.pi/6,'#a33b20'), (math.pi/3,'#1f5c8b'), (math.pi/2,'#123a57')]
    labs = [r'$\theta=0$  (Door R)', r'$\theta=\pi/6$', r'$\theta=\pi/3$', r'$\theta=\pi/2$  (Door C)']
    store = {}
    for (th,c),lb in zip(cols, labs):
        ts, Vs, g12 = run(th); store[th] = (ts, Vs, g12, c, lb)
        axes[0].plot(ts, Vs, color=c, lw=1.5, label=lb)
    axes[0].set_xlabel('t'); axes[0].set_ylabel(r'$V=-\frac{1}{2}\log\Delta$')
    axes[0].set_xlim(0, 6); axes[0].legend(fontsize=7.5, frameon=False)
    axes[0].set_title(r'(a) descent rate $=\cos\theta$')
    tsC, _, gC, cC, _ = run(math.pi/2, N=22000)[0:3] + ('#123a57', '')
    tsE, VE, gE = run(1.5, N=22000)
    axes[1].plot(tsC, gC, color='#9dc0d8', lw=0.9, label=r'$\theta=\pi/2$: undamped')
    axes[1].plot(tsE, gE, color=WARN, lw=1.4, label=r'$\theta=1.5$: damped')
    axes[1].set_xlabel('t'); axes[1].set_ylabel(r'$\gamma_{12}(t)$')
    axes[1].set_xlim(0, 44); axes[1].legend(fontsize=7.5, frameon=False, loc='upper right')
    axes[1].set_title('(b) oscillation under a decaying envelope = evanescence')
    fig.text(0.5, 0.015,
        'Figure 3. The Kähler dial. One potential, one state space: $\\dot{x}=-e^{J\\theta}\\,\\mathrm{grad}\\,V$ with\n'
        '$J_a v = a\\times v$ the complex structure, whose generator is the axis\'s own plane $P_a=Ia$. Exactly\n'
        '$dV/dt=-\\cos\\theta\\sum_i|P_i\\nabla_iV|^2$: $\\theta=\\pi/2$ conserves $V$ (unitary circulation),\n'
        '$\\theta=0$ descends (dissipative settling), and the band between oscillates while decaying.\n'
        'By RULING-1 the totality is not on this dial at all; $\\theta$ is a seat predicate.',
        ha='center', va='bottom', fontsize=8)
    fig.subplots_adjust(bottom=0.33, wspace=0.28)
    fig.savefig('fig3_kahler_dial.png'); plt.close(fig); print("fig3 ok")
fig3()
print("all figures written")
