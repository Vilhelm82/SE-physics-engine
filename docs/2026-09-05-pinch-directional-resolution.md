# PINCH-1 continuation: retain the merging direction

5 September 2026

**PINCH-1 provides the missing input for extending the regular horizon work to rank loss.** Its timelike transverse limit and two null sheet limits sit inside one exact family. Retaining the direction of ruler merger gives a replacement frame with determinant exactly \(-1\), including at the pinch.

The algebra below is derived from PINCH-1's declared form and oriented frame. Treating the retained direction as part of the state at the pinch is a **proposed extension**. Selecting it from native dynamics remains a theorem target.

## Exact directional identity

Let \(d=\det[c,h,g]\) be the signed frame volume in the oriented ambient form \(q=\operatorname{diag}(1,1,-1)\). Then \(d^2=\delta\). For \(d\ne0\), define

\[
\xi=\frac{b-a\gamma}{d},\qquad v=\frac{g-\gamma h}{d}.
\tag{1}
\]

The letter \(v\) here denotes a vector, not the advanced coordinate in the metric note. Direct expansion gives

\[
\delta-(b-a\gamma)^2=(1+a^2)(1-\gamma^2).
\]

Therefore

\[
\boxed{F=\frac{1-\xi^2}{1+a^2},\qquad
\eta=\frac{a^2+\xi^2}{1-\xi^2}.}
\tag{2}
\]

Equation (2) is exact wherever the original ratios are defined. It makes the retained direction, including its coefficient rather than just vanishing order, explicit.

## A regular replacement for the collapsed frame

From the native pairings,

\[
q(h,v)=0,\qquad q(c,v)=\xi,\qquad q(v,v)=F.
\]

Thus the Gram matrix in \((c,h,v)\) is

\[
\boxed{
\widetilde G=
\begin{pmatrix}-1&a&\xi\\a&1&0\\\xi&0&F\end{pmatrix},
\qquad \det\widetilde G=-(1+a^2)F-\xi^2=-1.}
\tag{3}
\]

Its signature is \((2,1)\): the \((c,h)\) block has one eigenvalue of each sign, and its Schur complement is \(1/(1+a^2)>0\). This remains true at \(F=0\).

The continuation is realizable, rather than only a formal ratio. Write \(p=1+a^2\), choose any finite \(\xi\), and set

\[
c=(0,0,1),\quad h=(\sqrt p,0,-a),\quad
v=(a\xi/\sqrt p,1/\sqrt p,-\xi).
\]

For each \(\sigma=\pm1\), the exact local chart

\[
\boxed{
\gamma=\sigma\sqrt{1-Fd^2},\qquad
b=a\gamma+\xi d,\qquad g=\gamma h+dv
}
\tag{4}
\]

is real and smooth for sufficiently small \(d\). It reproduces the original Gram matrix and signed determinant \(d\). At \(d=0\), it gives \(\gamma=\sigma,b=\sigma a,g=\sigma h\), precisely PINCH-1's pinch. Meanwhile \(\det[c,h,v]=1\).

For \(d\ne0\), this is a change of basis in the original model. At \(d=0\), it retains the normalized difference discarded by the merged frame. It replaces each collapsed state by its finite merger directions. The chart covers finite \(\xi\); approaches with \(|\xi|\to\infty\) need separate treatment.

## PINCH-1's limits become one normal family

The dual normal satisfying \(q(\nu,c)=1\), \(q(\nu,h)=q(\nu,v)=0\) is

\[
\boxed{
\nu=-Fc+aFh+\xi v
=\left(\frac a{\sqrt{1+a^2}},
\frac\xi{\sqrt{1+a^2}},-1\right),\qquad q(\nu,\nu)=-F.
}
\tag{5}
\]

Consequently:

- \(\xi=0\) gives PINCH-1's transverse limit \(\nu_*=(-c+ah)/(1+a^2)\).
- \(\xi=\pm1\) gives its two null limits, with \(\nu_*=(\nu_{+1}+\nu_{-1})/2\).
- \(|\xi|<1\) gives the intervening timelike normals; \(|\xi|>1\) gives spacelike normals.

At a fixed pinch depth \(a_0\), the attainable finite lapse limits are

\[
F_0=\frac{1-\xi_0^2}{1+a_0^2}.
\]

The compact-side limits fill \([0,1/(1+a_0^2)]\), including the horizon endpoints; equation (4) realizes every value. This recovers PINCH-1's \([0,1]\) family at zero depth and extends it to arbitrary depth and the hyperbolic side.

For a transverse Gram path in PINCH-1, \(b-a\gamma=O(s)\) and \(\delta=O(s)\) with a nonzero leading coefficient on the admissible side, so \(\xi=O(\sqrt{s})\to0\). Along its horizons, (2) forces \(\xi^2=1\). The different limits thus correspond to different retained merger directions.

For smooth ambient frame curves passing through the pinch, \(\delta=d^2\) implies \(\delta'=0\). Differentiating the exact identity above then gives \(\gamma'=0\) there. This connects the contact-order classification to the regularity of the lifted frame: a transverse Gram parameter need not itself be a smooth parameter for that frame. The transverse limit remains realized after a suitable reparameterization.

## Extension of the horizon comparison

On this resolved state space, the two null directions are regular level sets:

\[
F=0\iff\xi=\tau,\qquad
dF|_{\xi=\tau}=-\frac{2\tau}{1+a^2}\,d\xi,
\qquad \tau=\pm1.
\tag{6}
\]

Here \(dF,d\xi\) are differentials; \(d\) in (1)–(4) is the signed volume. The label \(\sigma\) selects the original \(\gamma\) branch; \(\tau\) selects a null limit over its pinch.

With the fourth positive direction and a smooth invertible soldering of the **resolved** frame, the earlier metric construction applies near these horizons, including \(d=0\), using (2) as its scalar \(F\). This is a conditional regular continuation of that metric family.

For the stationary radial specialization, its surface-gravity expression becomes

\[
\kappa_H=\frac{A_H}{2B_H}\partial_rF\big|_H
=-\frac{\tau A_H}{B_H(1+a_H^2)}\partial_r\xi\big|_H,
\tag{7}
\]

which is finite for smooth radial fields and finite nonzero lapse factors. This establishes a regular coefficient, without selecting those fields dynamically.

The trapping comparison still reads

\[
F_{\rm trap}=\frac{1-\xi^2}{1+a^2}+\frac{2\theta_K}{\rho}.
\tag{8}
\]

Its hypotheses now concern the resolved normal and a smooth, cut-compatible soldering. The highest-value native continuation has two concrete parts: **derive the evolution of the merger direction \(\xi\), and derive its associated area transport \(\theta_K\)**. Together these determine whether evolution reaches a timelike, null or spacelike pinch direction and whether it coincides with gravitational trapping.

## Evidence

The input [PINCH-1 runner](/home/williaml/seated-root/pinch1_horizon_branch_meet.py) was read and run: **33/33 reported checks passed**. The independent [directional derivation runner](/home/williaml/seated-root/pinch_directional_resolution.py) passes **42 exact symbolic checks**, covering the exact identity, both branch realizations, the regular replacement frame, all three causal types and the surface-gravity coefficient. Its results are in [the check record](/home/williaml/seated-root/docs/pinch-directional-resolution-checks.json).

Related project results: [horizon-crossing metric](/home/williaml/seated-root/docs/2026-09-05-horizon-crossing-metric.md) and [rest-frame/trapping comparison](/home/williaml/seated-root/docs/2026-09-05-rest-frame-trapping-equivalence.md).

The next target is developed in [pinch selection and area transport from a load](/home/williaml/seated-root/docs/2026-09-05-pinch-load-transport.md): native rate selection, a Riccati law on the resolved pinch, passive response conditions, and a solved finite-time crossing.
