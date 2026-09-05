---
title: "Seated Root: mathematical foundation audit"
date: "5 September 2026"
lang: en-AU
---

# Mathematical foundation audit

**Purpose:** assess the framework against its intended destination: a physics model capable of earning novel predictions, especially about singularities and black holes. The current construction is a candidate design, and every mathematical choice is open to replacement if a stronger construction does its job better.

**Scope:** the live `/home/williaml/seated-root` checkout at `45c744365ace464bd3e486a1d1b4282e75e35516`, including its untracked curvature and inertia documents. The current primitives and latest handoff take precedence over paper v0.5 where they conflict. This is a mathematical and architectural audit, not a complete formal verification of every runner or an observational review of the cosmology conjecture. No model files were changed.

## 1. Verdict

**The model can become a useful framework without physics labels. It is not yet one simply by deleting the labels.** Its most reusable idea is that a state, a reference choice, and the presentation available from that reference are different objects. Gram invariants, observation fibres, transformations, orientation-sensitive operations and transport can all be studied independently of light, mass, temperature or quantum mechanics.

The current construction nevertheless selects several special mathematical worlds: quadratic forms, particular realifications, three-dimensional determinant identities, particular covering representations, round angular measures, and a restricted gravitational geometry. These choices support genuine conditional theorems. They do not all follow from the verbal primitives, and they cannot all survive unchanged in a framework for arbitrary use cases.

**The strongest route is a general mathematical core with explicitly selected realizations, followed by a physical reconstruction and a closed dynamics.** Existing mathematics already supplies most of that core. The promising original work is in the reconstruction, compatibility, singularity and dynamics theorems that connect its parts. There is no demonstrated obstruction requiring a replacement for linear algebra, Clifford algebra or differential geometry.

The audit also found concrete errors of domain or inference, not merely opportunities for generalization:

- Generic coplanar Gram-coordinate blowup is not an intrinsic curvature blowup of the oriented frame quotient.
- Zero lapse can coincide with rank loss along a fully realizable finite path.
- A globally claimed coboundary fails at the fixed poles, and a closed-loop factor of two cannot be removed by a genuine vertex coboundary.
- The four-radical Kummer extension is closed under axis permutations, but not the full signed click group.
- Pin kernel sign, geometric orientation and Clifford block label are distinct. Their identification needs maps that the current equations do not supply.
- Legacy measure, thermal and quantum-readout results have not all been transferred to the rebuilt primitives.

These findings change the immediate priority. Before interpreting a block as an expanding interior, repair the domains and specify the dynamics and observation map that would make “expanding” a calculable property.

## 2. What counts as a stronger foundation

An input is not circular merely because it restricts the model. A predictive theory must exclude possibilities. The relevant questions are whether the restriction has an independent reason, whether the result could have failed under the permitted alternatives, and whether the chosen representation already carries the physical conclusion being advertised.

This audit uses five criteria: explicit assumptions; preservation of the intended observations and operations; ability to represent alternatives that can fail; a manageable path to well-posed dynamics; and the capacity to produce an uncalibrated observable. Greater generality alone is not an improvement if it only introduces unrestricted functions that can reproduce any answer.

Three distinctions organize the findings:

1. **Representation:** an efficient way to calculate a structure already justified. Clifford algebra is excellent at this job.
2. **Reconstruction:** a theorem showing that independently stated operational assumptions force that structure. Several of these theorems are still missing.
3. **Physical realization:** a map from mathematical states, operations and parameters to preparations, measurements, clocks, sources and events. An isomorphism of algebras is not this map.

The current project is strongest at representation and conditional identities. Its main opportunity is reconstruction and physical closure.

## 3. Assessment of the mathematical bodies

### 3.1 Incidence, reference choices and imaginary-angle pivots

**Verdict: retain the relational idea; make the analytic and dimensional choices explicit.**

P1–P9 provide a root, lines and planes, reference choices, partial presentations and imaginary-angle pivots. Incidence and reference selection are reusable outside physics. Exactly three lines, a particular meaning of “angle”, and analytic continuation of ordinary rotations are additional structure.

T1 correctly turns a continued ordinary rotation into a real hyperbolic matrix:

\[
S^{-1}R(i\lambda)S=
\begin{pmatrix}\cosh\lambda&\sinh\lambda\\\sinh\lambda&\cosh\lambda\end{pmatrix},
\qquad S=\operatorname{diag}(1,i).
\]

Its invariant symmetric forms are proportional to \(\operatorname{diag}(1,-1)\). This is sound given the rotation representation, analytic continuation, real slice and quadratic-invariant question. Bare lines through a point do not choose those ingredients. [T1–T3 runner](/home/williaml/seated-root/prim_t1_t3_pivot_group.py:43).

T2's count of zero or two boost planes holds for its diagonal realifications of an orthogonal coordinate frame. It does not classify arbitrary pair-spans of nonorthogonal state lines. The Gram matrix

\[
G_* = \begin{pmatrix}1&-2&-2\\-2&1&-2\\-2&-2&1\end{pmatrix}
\]

has eigenvalues \((3,3,-3)\), while every pair restriction has determinant \(-3\). Three unit spacelike lines can therefore span three hyperbolic planes. This does not invalidate the coordinate-generator count; it prevents promoting that count to a larger claim. [T2 enumeration](/home/williaml/seated-root/prim_t1_t3_pivot_group.py:111).

For a general core, allow a declared collection of reference objects and admissible transformations. Choose continuous Lie groups, analytic continuation or a fixed dimension only when the use case justifies them. There is no need to invent a new group theory to make this separation.

### 3.2 Bilinear gluing, Gram geometry and Schur complements

**Verdict: keep the linear algebra; qualify the foundational derivation.**

T7 starts with a symmetric bilinear form and supplies its three plane restrictions, including the diagonal normalizations \((-1,+1,+1)\). These restrictions uniquely determine its six coefficients. The proved result is a gluing theorem **inside the category of bilinear forms**. The words “compact”, “hyperbolic”, “reference” and “collapse” have not independently selected that category. [T7 gluing](/home/williaml/seated-root/prim_t7_seat_form.py:44); [primitives corrections](/home/williaml/seated-root/docs/2026-09-04-PRIMITIVES-v0.md:132).

This is a substantive gap rather than a demand for assumption-free mathematics. For example,

\[
L_\epsilon(x,y,z)=x^2+y^2-z^2+
\epsilon\frac{x^2y^2z^2}{(x^2+y^2+z^2)^2}
\]

on nonzero vectors, extended by zero at the origin, has exactly the same restrictions on all three coordinate planes. It is homogeneous of degree two but nonquadratic for \(\epsilon\ne0\): the parallelogram defect for \(u=(1,1,0),v=(0,0,1)\) is \(2\epsilon/9\). It is a counterexample to selection by those plane restrictions alone, not to an additional full linear-isometry axiom. A polarization law or a sufficiently strong operational linearity/isometry principle could exclude it.

Cone and Lorentz–Finsler geometry are established comparison classes if quadratic response cannot be justified. They allow causal cones and anisotropic response without presupposing a quadratic metric. They also add freedom; adoption is warranted only if that freedom serves an identified measurement or dynamics problem. [Javaloyes–Sánchez](https://arxiv.org/abs/1805.06978).

Once the form is justified, Schur complements are an excellent fit. With

\[
G=\begin{pmatrix}-1&a&b\\a&1&\gamma\\b&\gamma&1\end{pmatrix},
\quad S=\begin{pmatrix}1&\gamma\\\gamma&1\end{pmatrix},
\quad v=(a,b)^T,
\]

the compact-ruler domain \(|\gamma|<1\) gives

\[
\eta=v^TS^{-1}v\ge0,\qquad
N^2=\frac1{1+\eta}=\frac{1-\gamma^2}{-\det G}.
\]

These identities should be retained. Their physical use as a clock lapse requires the associated observer/clock construction. [T7d](/home/williaml/seated-root/prim_t7d_tilt.py:45).

The projected ruler plane \(c^\perp\) and the plane \(\operatorname{span}\{\hbar,G\}\) are different when the rulers have depth. The current primitives and runners sometimes describe both as the seat's real face-on plane. Resolve that distinction in the definitions: it affects which plane is always compact and which can become Lorentzian. [T7 plane choice](/home/williaml/seated-root/prim_t7_seat_form.py:68).

### 3.3 Euclidean elliptope, Cayley surface and signed state space

**Verdict: retain as exact realizations, not a universal state space.**

For three Euclidean unit vectors modulo \(O(3)\), the unit-diagonal PSD Gram matrix is the right state description. Its determinant

\[
\Delta=1-x^2-y^2-z^2+2xyz
\]

and rank strata have their usual algebraic meaning. Correlation-matrix quotient geometry already provides a general framework for arbitrary numbers of vectors and bounded ranks. [Chen](https://arxiv.org/abs/2401.03126); [Cayley runner](/home/williaml/seated-root/cayley.py:17).

The rebuilt Lorentzian frame is a different real sector. Its realizability condition is

\[
|\gamma+ab|\le\sqrt{(1+a^2)(1+b^2)}.
\]

Thus three free angle parameters do not mean three independent unrestricted Gram entries. The parameterization by \((t,\ell_1,\ell_2)\) enforces the condition. In it,

\[
\det G=-D^2,
\qquad D=\cosh\ell_1\cosh\ell_2\sin t.
\]

The real rank-one nodes of the Euclidean sector are absent here because a negative unit line cannot coincide with a positive one. This is a real-domain change, not removal of all algebraic singularities: \(\det G(a,b,\gamma)=-\Delta(ia,ib,\gamma)\), so the complexified determinant surfaces are equivalent. [T7 state realization](/home/williaml/seated-root/prim_t7_seat_form.py:84).

For a general framework, maintain configuration space, its symmetry quotient, its orientation-forgetting quotient, and its chosen coordinate chart separately. Their singular sets need not coincide.

### 3.4 Observation maps and loss of information

**Verdict: promote this to the central abstraction, with explicit domains and fibres.**

The current unpivoted normalized-angle readout is \(p_c(t,\ell_1,\ell_2)=\cos t\). Both depths are invisible throughout the regular domain. Noninjectivity therefore occurs generically, not only on the determinant branch locus. At \(\sin t=0\), the differential of this scalar readout also loses rank; that is a further event. [T7 readout](/home/williaml/seated-root/prim_t7_seat_form.py:96).

This is the right place to generalize the model. Ask which states are indistinguishable under a specified family of operations and observations, and which additional readings recover the missing data. Noninjectivity alone supplies neither quantum superposition nor a rule for evolution.

An observation map does not automatically inherit autonomous dynamics. On \(X=\mathbb R^2\), take \(p(x,h)=x\) and \(\dot x=h,\dot h=0\). The states \((0,1)\) and \((0,-1)\) look identical under \(p\), but their immediate observed velocities differ. This simple countermodel applies to the proposed seat-derived laws: either the hidden data must be recovered, or the presentation needs memory or additional state.

### 3.5 Clifford algebra, Pin/Spin, Wigner rotation and the sheets

**Verdict: retain as conditional representation machinery; separate the covering maps.**

Clifford algebra is the universal associative algebra for a specified quadratic relation \(v^2=q(v)1\). Its use after deriving \(q\) is not inherently circular. It does not independently derive \(q\), identify a spinor with matter, or supply probabilities. Pin/Spin is similarly well suited to reflection lifts and rotational holonomy. [Meinrenken, §3](https://www.math.toronto.edu/mein/research/keio_06.pdf).

The Wigner/area identity is a good fit for the chosen connection on \(H^2\). It remains a theorem about that connection; a dynamical state path must be shown to use it. [T5 area runner](/home/williaml/seated-root/prim_t5_area_and_lifts.py:5).

Four signs currently sit too close together: geometric orientation \(D\); the Pin kernel \(u\mapsto-u\); permutation parity of orthogonal reflection lifts; and the central Clifford block \(\chi=\pm1\). The exact trace is

\[
\mathcal V_\chi=2\chi D.
\]

It reads their product. It does not identify the factors. With central unit pseudoscalar \(I\), every element of \(Cl(2,1)\), including every ordinary Pin product, commutes with \(e_\pm=(1\pm I)/2\). Multiplying by the kernel sign cannot exchange the blocks. Pin reflection of vectors uses the twisted adjoint; block exchange is a different operation. [T5c](/home/williaml/seated-root/prim_t5c_corrections.py:47); [T8a](/home/williaml/seated-root/prim_t8a_fourth_direction.py:38).

The Boolean flip/parity picture also needs its domain: mutually orthogonal reflection lifts anticommute. For nonorthogonal vectors the relation is \(uv+vu=2q(u,v)\), and their reflection products generally do not commute. Antisymmetry of one trace is weaker than equality of whole operators up to sign. [T5b general correction](/home/williaml/seated-root/prim_t5b_parity.py:104).

The octahedron is useful combinatorics for three binary choices. Its named faces are a proposed physical interpretation of that combinatorics. Arbitrary numbers of references lead to different combinatorics and different Clifford decompositions; neither the octahedron nor its two-block realization is universal.

### 3.6 The fourth generator and its adjoint

**Verdict: a useful minimal graded extension, with three unfinished physical identifications.**

T8a correctly solves the stronger problem of finding an operator that anticommutes with all three Clifford vectors. Its general form is

\[
\Gamma_4=\begin{pmatrix}0&\alpha J\\\beta J&0\end{pmatrix},
\qquad J^2=-1,
\qquad\Gamma_4^2=-\alpha\beta\,\mathbf1_4.
\]

Both signs exist. But block exchange alone does not force this anticommutation condition: \(E=\begin{pmatrix}0&1_2\\1_2&0\end{pmatrix}\) exchanges the central blocks while commuting with two existing vector generators. “An internal block-coupling operator” and “an additional spacetime vector direction” therefore require a further identification. [T8a solution](/home/williaml/seated-root/prim_t8a_fourth_direction.py:48).

T8b's invariant real spinor form is symplectic. Its adjoint calculations are valid, but symplectic self-adjointness is not a positive measurement rule. For commuting real spinor components,

\[
B^T=-B\Rightarrow u^TBu=0,
\qquad A^{\dagger_B}=A\Rightarrow(BA)^T=-BA\Rightarrow u^TBAu=0.
\]

An operational state/effect pairing is still needed. The parameter identification “one seated line, therefore equal relative spinor-form sign” is also an extra bridge: the runner checks the consequences of \(s=+1\), not its derivation from the number of clocks. [T8b form](/home/williaml/seated-root/prim_t8b_adjoint.py:57); [T8c sign table](/home/williaml/seated-root/prim_t8c_reciprocity.py:94).

### 3.7 Galois/Kummer extensions, invariants and monodromy

**Verdict: retain the algebra; repair the base field and acting group.**

Over \(F=k(x,y,z)\), \(k=\mathbb Q(i)\), the four independent square classes \(\Delta,1+x,1+y,1+z\) generate a degree-16 extension. Over the permutation-invariant field \(F^{S_3}\), the corresponding group has order 96. These are useful exact algebraic results. [Latest descent receipt](/home/williaml/seated-root/thm_galois_deck_descent.py:6).

The full signed click group also sends \(x\) to \(-x\). It therefore sends \(\sqrt{1+x}\) to a square root of \(1-x\), whose divisor is absent from the four-class extension. The field is not closed under all those clicks. The smallest editorial repair is to name the acting subgroup \(S_3\). If full signed-click closure is required, adjoining all six \(\sqrt{1\pm x},\sqrt{1\pm y},\sqrt{1\pm z}\), together with \(\sqrt\Delta\), gives a generic rank-seven extension. The resulting group over the signed-click invariant base is \(C_2^7\rtimes S_4\), of order 3072. This completion is a mathematical option, not a new physical claim.

A Bargmann expression can descend to an orientation subfield without being invariant under every physical relabelling. Specify the field lift and the transformation of the observable. Complex branch monodromy also does not automatically imply a real path obstruction: the interior of the Euclidean elliptope admits a globally chosen positive \(\sqrt\Delta\). [Monodromy runner](/home/williaml/seated-root/thm_b_monodromy.py:14).

### 3.8 Symplectic/Kähler geometry and the proposed dissipator

**Verdict: retain on the explicit lifted phase space; do not silently put it on the three-angle quotient.**

The product \((S^2)^3\cong(\mathbb{CP}^1)^3\) is six-dimensional and can carry the declared product Kähler structure. A three-dimensional state quotient cannot itself carry a nondegenerate symplectic two-form. The inherited structure must be identified through a specified reduction, Poisson bracket or symplectic leaves. Those choices may also require extra phase variables beyond P4. [Kähler construction](/home/williaml/seated-root/kahler.py); [paper §3](/home/williaml/seated-root/PAPER-seated-root-v0.5.md:208).

The bath and channel calculations are conditional on a spectral density, coupling law, bath statistics and golden-rule conventions. DEBT-2b's functional equation can select the coupling inside those premises; it does not derive irreversible statistical mechanics from the static frame. The inspected files do not implement a general contact-geometry or Lindblad/GKSL foundation. Those should remain proposed or queued tools, not inventory items already earned. [DEBT-2](/home/williaml/seated-root/debt2.py); [DEBT-2b](/home/williaml/seated-root/debt2b.py).

**Gradient/Hamiltonian flows and log-determinant energies.** The repository already supplies a revealing comparison: `doors.py` chooses different dynamics and potentials on the same frame geometry. The barrier \(-\tfrac12\log\Delta\) diverges at the boundary; \(-\Delta\) is smooth there. Boundary exclusion or crossing therefore depends on the energy, kinetic structure and evolution law. It is not a consequence of the incidence data. Calling the log determinant a Gaussian information quantity additionally assumes a Gaussian statistical model. Retain these flows as comparators and possible independently motivated laws, rather than treating an attractive invariant as the uniquely selected energy. [Declared dynamics](/home/williaml/seated-root/doors.py:5).

**A useful reduction already available.** If the product-sphere Poisson structure is accepted, the oriented shape variables obey \(\{x,y\}=V\), \(\{x,z\}=-V\), \(\{y,z\}=V\), with \(V^2=\Delta\). The function \(x+y+z\) is a Casimir. This gives a generically rank-two Poisson structure and two-dimensional symplectic leaves, rather than a Kähler structure on the full odd-dimensional quotient. It still needs a Hamiltonian. This is a practical existing mathematical route for the Euclidean realization.

**Symplectic squeezes and the positive floor.** For a declared positive matrix \(C_0=A1_2\), squeezing preserves \(\det C=A^2\). The real symmetric-matrix representation and its doubled angles are reusable without quantum labels. But arbitrary \(A>0\) does not establish a universal nonzero quantum floor: choosing \(A=\hbar/2\) and interpreting the matrix as quantum covariance requires a physical rule. The full covariance expression \(q q_\perp-\mathrm{cov}^2=A^2\) is invariant under rotations; only the weaker marginal product \(q q_\perp\ge A^2\) saturates solely on principal axes. These two meanings of minimum uncertainty should be distinguished. [Squeeze construction](/home/williaml/seated-root/thm_g.py:56); [comparison language](/home/williaml/seated-root/thm_g.py:216).

### 3.9 Measures, cocycles and the trichotomy

**Verdict: retain the transformation identities; rederive carriers, measures and boundary claims.**

The exponents depend on what is being integrated. The legacy c/G construction uses round \(S^2\) solid angle, with a \(s^{-2}\) aberration factor. A native \(2+1\) construction has null directions on \(S^1\), with angular factor \(s^{-1}\). Writing \(u=\cos\theta\) still gives \(du'/du=s^{-2}\), but uniform \(du\) is not uniform angular measure on \(S^1\): the latter has weight \(1/\sqrt{1-u^2}\). The old census and derived bath exponent therefore need a carrier/measure derivation after the dimensional rebuild. A subsequently justified physical \(3+1\) realization might recover them. [G runner](/home/williaml/seated-root/thm_g.py); [one-norm census](/home/williaml/seated-root/sect1_one_norm.py).

There is also an exact domain correction. For

\[
g_L(u)=\frac{u+\tanh L}{1+u\tanh L},\qquad
\sigma_L(u)=\cosh L+u\sinh L,
\]

the expression \(\sigma_L(u)=b(g_Lu)/b(u)\), \(b(u)=(1-u^2)^{-1/2}\), is valid on \(-1<u<1\). At the fixed poles, \(\sigma_L(\pm1)=e^{\pm L}\); any finite nonzero coboundary at a fixed point equals one. It is therefore not a global trivialization on the carrier including the poles. [PRED-1 cocycle](/home/williaml/seated-root/pred1_cross_seat.py:43).

Likewise a genuine multiplicative vertex coboundary telescopes to one on every closed loop. The runner's closed-loop factor two can be a covering transfer/normalization issue, but cannot be absorbed as such a coboundary. Its second inverse-cover branch is also incorrect: adding \(\pi/2\) to \(\tfrac12\arccos u\) produces \(\cos2\phi=-u\). The other branch modulo \(\pi\) is \(\pi-\tfrac12\arccos u\). [Cover claims](/home/williaml/seated-root/pred1_cross_seat.py:77).

Groupoid cocycles, density bundles and pushforward/transfer operators are appropriate general tools. They preserve exactly the distinctions between observable values, measure Jacobians, multiplicities and path holonomy that these findings require.

**Radon–Nikodym, moment and harmonic-analysis machinery.** If \(d\mu'/d\mu=f^{-n}\), then \(\int f^n d\mu'=\mu(\Omega)\). Keep this general identity. The separate question is why the physical readout is that power of that density. The Doppler-family uniqueness of blind moments is conditional on its carrier and measure, not a theorem for arbitrary response families. Legendre/Mellin moment identities and power means remain useful analytical tools after the detector weighting is chosen. [THM-RN's own boundary](/home/williaml/seated-root/thm_rn.py:6).

**Exact sequences and gauge language.** Splitting logarithmic readings into their sum and difference is valid linear algebra. The kernel of the sum map is not automatically a physical gauge freedom: the individual readings can still detect the difference. Gauge identification depends on the complete admitted observation algebra. [Accounting bridge](/home/williaml/seated-root/bridge_dbp.py:56).

### 3.10 Thermal periods, quantum readouts and operational predictions

**Verdict: retain periodicity and trace identities; keep their statistical and physical bridges explicit.**

T4's vector/spinor periods and logarithmic asymptotics are useful algebra. Temperature additionally uses the identification between imaginary-time period and \(\hbar/(k_BT)\), a physical clock/radial reconstruction and, for the first-law route, area entropy. Agreement of routes sharing these ingredients is an internal consistency result. It is not several independent derivations of the ingredients. [T4 temperature assembly](/home/williaml/seated-root/prim_t4_hawking_period.py:155). Fermionic antiperiodicity at a fixed inverse temperature prevents interpreting a spinorial minus sign alone as half the ordinary temperature. A thermal claim also needs a state and time evolution with the appropriate correlation-function boundary condition; a spin representation's period alone does not supply them. [Laine–Vuorinen, thermal field theory](https://arxiv.org/abs/1701.01554).

The operational PRED-1 runner directly supplies qubit registers \(P_i=(1+a_i\cdot\sigma)/2\). It establishes what a joint directed operation on those registers can read. It does not compute the entire state-to-presentation-to-preparation channel that would prove the same experiment is available from the rebuilt primitives. The older Pauli/Bloch carrier must therefore be reconstructed or retained as a declared realization. [Operational chain](/home/williaml/seated-root/pred1_operational_chain.py:11); [physical protocol](/home/williaml/seated-root/pred1_physical_protocol.py:109).

The ordinary Bargmann/Pancharatnam and graph-cycle calculations remain strong invariant tools on that carrier. Binet–Cauchy relates products of signed volumes to cross-Gram minors; anchored cycle decompositions work on charts with the required nonzero overlaps and reference volume. These identities describe what a representation forgets. They do not select physical effects. Use \(\operatorname{Im}B\) for a continuous signed algebraic separator; \(\arg B\) may be \(\pi\) on negative-real coplanar cases and is undefined when \(B=0\). The project's all-cycle quadratic result also has known mathematical comparison work, already acknowledged by its runner. [Graph receipt](/home/williaml/seated-root/graph_cocycle.py:7); [projective-frame invariants](https://arxiv.org/abs/1312.5393).

For a physics-neutral statistical layer, states, effects, transformations and system composition should be specified separately. General probabilistic theories provide an existing comparison class containing classical, quantum and other possibilities; they allow the model to ask which further assumptions select a quantum realization. Merely having indistinguishable presentations does not select it. [Barrett](https://arxiv.org/abs/quant-ph/0508211).

**Oscillator elimination and thermal baths.** The generalized-Langevin calculation is appropriate for the declared oscillator bath, coupling and initial Gaussian thermal moments. Those moments supply its fluctuation–dissipation relation. Finite baths recur; an appropriate continuum/weak-coupling limit is a further step. Under the native isotropic \(S^1\) alternative, the presented monochromatic bath has an arcsine frequency density rather than a flat band, so its memory kernel changes. Likewise the golden-rule functional equation with angular exponent \(Q\) forces the product \(c(\omega)^2g_0(\omega)\propto\omega^Q\) under its uniform-response and resonance premises. It does not independently determine density and coupling separately, the overall strength, or Bose statistics. [Bath construction](/home/williaml/seated-root/debt2.py:27); [coupling theorem](/home/williaml/seated-root/debt2b.py:43).

One further inference in paper §4 needs a dynamics premise: a Lorentz-vector current with fixed invariant norm need not be conserved. In units \(c=1\), \(n^\mu=J(\cosh at,\sinh at,0,0)\) has constant norm but \(\partial_\mu n^\mu=Ja\sinh at\ne0\). Covariance does not supply a continuity equation. [Paper current claim](/home/williaml/seated-root/PAPER-seated-root-v0.5.md:272).

### 3.11 Hodge duality, constitutive response and impedance

**Verdict: retain the conditional signature theorem; generalize the starting response.**

T8c's chain is valid under its stated constitutive condition:

\[
\text{compact duality on two-forms}
\Rightarrow\star^2=-1
\Rightarrow n_-\text{ odd}
\Rightarrow\Gamma_4^2=+1
\]

for the two four-generator candidates it compares. The unfinished step is that the seat's channels realize this duality. Moreover \(n=4\) follows from choosing degree-two fields and requiring Hodge duality to act internally; in general that condition is \(n=2p\). This is a meaningful conditional selection, not an unconditional dimension theorem. [T8c](/home/williaml/seated-root/prim_t8c_reciprocity.py:32).

The current construction \(K=Z^{-1}\star_q\) builds the metric cone into the response. It correctly checks compatibility and leaves \(Z\) free. If the aim is to earn the cone from a load, the stronger comparison is a general constitutive map first, then its characteristic polynomial. In four dimensions a local linear response has 36 components; metric response is a restricted sector. Premetric electrodynamics derives a general Fresnel polynomial and states conditions, including symmetry and closure, under which a conformal metric emerges. Those conditions must be derived or declared, not hidden. [Rubilar](https://arxiv.org/abs/0706.2193); [Hehl and collaborators](https://arxiv.org/abs/gr-qc/0506042); [current response runner](/home/williaml/seated-root/thm_k_response_map.py:14).

Even a recovered cone cannot select all coupling data. Locally at the principal-symbol level,

\[
H=\lambda\star_qF+\vartheta F
\]

has the same metric characteristic cone for nonzero scalar \(\lambda\): \(k\wedge F=0\) removes the \(\vartheta\) term from the other wave equation. Normalization and an axion-like component therefore require additional response information. Interface/source effects can distinguish data invisible to homogeneous bulk characteristics. [Obukhov–Hehl](https://arxiv.org/abs/physics/0504172).

HUNCH-Z0 correctly removes dependence on electrical units, but its numerical check inserts \(\alpha\) before constructing \(Z_0\). It verifies an identity, not an independent value of \(\alpha\). A prediction needs a dimensionless output obtained without that target input. [HUNCH-Z0](/home/williaml/seated-root/hunch_z0_impedance.py:47).

Finally, \(N^2=1/(1+\eta)\) can be written as a voltage-divider ratio on its positive domain. That useful analogy does not by itself specify a physical impedance, effort/flow pairing, clock or propagator. A passive real divider is also not automatically a model for the sector where \(\eta\le-1\). A constitutive law must say whether it is scalar/tensorial, local/nonlocal, frequency-dependent and dissipative. \(V=IZ\) becomes dynamics only when these objects and the derivative/current are defined.

### 3.12 Mechanical connection, quotient curvature and rank strata

**Verdict: retain the connection and inertia detector; correct the intrinsic-divergence interpretation.**

For the Euclidean product-round configuration \(E=(S^2)^3\), the repo derives

\[
\mathbb I=3\mathbf1_3-AA^T,
\qquad
F(X,Y)=2\mathbb I^{-1}\sum_iX_i\times Y_i
\]

for horizontal vectors. The locked-inertia theorem \(\det\mathbb I=0\iff\operatorname{rank}A=1\) is sound. These are established geometric-mechanics objects and well suited to the problem. [Connection proof](/home/williaml/seated-root/docs/GRAM_SUBMERSION_CURVATURE/GRAM_SUBMERSION_CURVATURE.md:27); [rank detector](/home/williaml/seated-root/docs/INERTIA_NODE_DETECTOR/INERTIA_NODE_DETECTOR.md:23).

At a rank-two coplanar frame the \(SO(3)\) action is still free and \(\mathbb I^{-1}\) is smooth. For horizontal unit vectors,

\[
\left\|\sum_iX_i\times Y_i\right\|\le1,
\qquad\|F(X,Y)\|\le2\|\mathbb I^{-1}\|.
\]

Thus the invariant connection curvature is locally bounded there. The code instead uses the lifts of coordinate-unit changes in the three Gram entries. Those lifts diverge at the orientation-forgetting fold. Summing their coefficients without the inverse quotient-metric contractions produces a coordinate-response norm, not an intrinsic tensor norm. Dividing a loop angle by coordinate area validates that same response. [Norm implementation](/home/williaml/seated-root/docs/GRAM_SUBMERSION_CURVATURE/GRAM_SUBMERSION_CURVATURE.md:186).

An exact witness makes the distinction concrete. At the coplanar trine with angles \(0,120^\circ,240^\circ\), \(\mathbb I=\operatorname{diag}(3/2,3/2,3)\). Choosing the horizontal orthonormal pair specified in the verification receipt gives

\[
\|F(X,Y)\|^2=\frac89,
\qquad K_{E/SO(3)}(X,Y)=\frac43.
\]

Both are finite on the coplanar stratum. A singular response to prescribed overlap changes may still be an interesting observable, but it must retain that protocol. It cannot establish a spacetime curvature singularity.

Rank-one nodes are different: the stabilizer changes and \(\mathbb I\) degenerates. They need a slice/stratified analysis, invariant norms and declared path classes. The present exponent argument establishes particular cancellation and big-O behavior, not a universal nonzero asymptotic coefficient for every approach. Nor can the positive-definite quotient proof be transferred to \(SO(2,1)\) without examining null directions, noncompact actions and the chosen metric on variations.

### 3.13 Gravitational profiles, transport, variational laws and Gauss–Codazzi

**Verdict: preserve the conditional successes and precise failures; replace the restricted carrier as the general dynamics.**

THM-H's rank-one presentation update, THM-I's spherical transport and THM-K's thermal derivation of the pinning remain useful conditional results. Their inputs include the physical radial/clock interpretation, null or mass-shell readings, source/entropy/equilibrium rules and transport. Two computational realizations of the same Lorentz scaffold do not independently derive that scaffold. [Paper gravitational sector](/home/williaml/seated-root/PAPER-seated-root-v0.5.md:336).

THM-J makes a valuable negative statement about its candidate energy. For

\[
E[\lambda]=\int r^2F(\lambda)(\lambda')^2dr,
\qquad s(\lambda)=\int\sqrt{F(\lambda)}d\lambda,
\]

the field equation becomes \((r^2s')'=0\) where the change of variable is regular. The required profile is obtained in this class only when \(s\) is affine in the pinned variable. Choosing \(F\) to enforce this would insert the target; the specified quotient energy fails. This does not rule out other fields, derivative laws or source closures. [THM-J](/home/williaml/seated-root/thm_j_dyn.py:8).

THM-M earns a first-order rotating response under its declared boost/superposition rule. THM-N kills that particular second-order continuation through its forbidden anisotropic mass aspect. It does not kill all rotating extensions, and a material shell's quadrupole is not automatically a Kerr-hole comparison. [THM-N scope](/home/williaml/seated-root/thm_n_kerr_quadrupole.py:1).

THM-O's identities are sound for the explicit metric

\[
ds^2=-dt^2+\delta_{ij}(dx^i-v^idt)(dx^j-v^jdt).
\]

With its sign convention, \(G_{nn}=e_2(K)\) and \(G_{ni}=-\tfrac12(\nabla\times\nabla\times v)_i\). NORMAL-1 then imposes four normal constraints. The six tangential/evolution equations are not supplied. Furthermore, a pointwise orthonormal set of rods does not prove spatial flatness: a curved spatial metric also has local orthonormal coframes. Integrability and transport between frames are additional data. [THM-O](/home/williaml/seated-root/thm_o_strain_law.py:140).

Its multipole obstruction is a useful no-go inside the stationary perturbative, unit-lapse, flat-slice sector around the spherical profile. It argues for generalizing the spatial metric and connection. It is not a reason to discard Gauss–Codazzi, whose role is to relate constraints and curvature once the carrier is specified. [Gourgoulhon](https://arxiv.org/abs/gr-qc/0703035).

## 4. Singularity and black-hole consequences

Four questions must remain distinct: does the observation map lose information; does the state/orbit space change stratum; does physical causal transport become incomplete; and does the evolution equation lose existence, uniqueness or regularity? A positive answer to one does not answer the others.

The claimed universal separation between zero lapse and rank loss is false. In the compact-ruler domain take \(0<s<1\),

\[
a=s,\qquad b=0,\qquad\gamma=1-s^4.
\]

Then

\[
\eta=\frac1{s^2(2-s^4)}\to\infty,
\quad N^2\to0,
\quad-\det G=s^2+2s^4-s^8\to0.
\]

The finite endpoint \((0,0,1)\) is on the rank-loss locus, and every preceding state has signature \((2,1)\). The generic fixed-depth path staying off rank loss remains correct. What fails is the quantifier “requires”. The right next theorem classifies allowed approaches by their vanishing orders and distinguishes finite endpoints from escape to infinity. [T7e overclaim](/home/williaml/seated-root/prim_t7e_two_degeneracies.py:83).

If the model retains its effective Schwarzschild/river metric and physical free-fall transport into the centre, the curvature invariant is \(12r_s^2/r^6\), and a radial rain trajectory reaches \(r=0\) in finite proper time:

\[
\Delta\tau=\frac{2r_0^{3/2}}{3\sqrt{r_s}}\quad(c=1).
\]

Finite angular observables or deck holonomy cannot establish singularity resolution while that metric and tidal detector law remain physical. The constructive opportunity is to derive where this effective reconstruction stops applying and what well-defined transport replaces it. Modern singularity theorems concern causal geodesic incompleteness, with explicit geometric hypotheses; they do not reduce to divergence of a coordinate or a chosen scalar. [Senovilla–Garfinkle](https://arxiv.org/abs/1410.5226).

A zero in the static observer's lapse is also not a ban on infall. The existing river geometry has unit ADM lapse and freely falling trajectories cross its horizon in finite proper time. Its spherical future-null expansions illustrate what the white/black distinction actually needs:

\[
\theta_\pm=\frac2r(v_r\pm1).
\]

Inflow and outflow \(v_r=\mp\sqrt{r_s/r}\) have the same even pinning but opposite trapping character. An internal block label has to be connected to this future-oriented propagation law; its sign alone cannot supply the connection. Local trapping horizons require null normals, expansions and area; event horizons additionally require global causal information. [Hayward](https://arxiv.org/abs/gr-qc/9303006).

Accordingly, the white-hole-sheet conjecture remains an admissible research hypothesis. Its “direction” test needs an event/clock map and a dynamical vector field or propagator before it can return a physical verdict. The present audit neither establishes nor rules out a future realization of that conjecture.

## 5. A general framework that preserves the useful idea

The following is a **proposal for organizing the mathematics**, not a replacement model already adopted.

Start with a configuration/state space \(X\), a collection of references \(s\), admissible operations \(T\), and observation maps

\[
p_s:U_s\subseteq X\longrightarrow Y_s.
\]

The domain \(U_s\) records where that reference's reading is defined. All compositions below are restricted to their common domains. Do not initially require \(X\) to be a sphere, an elliptope, a Hilbert space or a particular Clifford module. Its structure is supplied by a chosen realization. A fibre of \(p_s\) is the set of states indistinguishable under that reading. An algebra of accessible observations consists of functions factoring through \(p_s\).

For a proposed transition from one presentation to another, the exact condition for a well-defined map \(\bar T\) is

\[
p_s(x)=p_s(x')\quad\Longrightarrow\quad
p_t(Tx)=p_t(Tx').
\]

If it holds, \(p_tT=\bar T p_s\). If it fails, the operation needs information absent from the starting presentation. This condition is useful in geometry, control, coarse-grained mechanics and statistical models, without physics labels.

For a flow \(\Phi_\tau\), an autonomous seat law requires the analogous factorization

\[
p_s\Phi_\tau=\psi_{s,\tau}p_s.
\]

On regular smooth domains, this is the familiar projectability/lumpability problem. The pushed-forward vector field must take the same value on every point of a fibre. Existing work supplies necessary and sufficient conditions for smooth surjective submersions. The model's rank-changing observation maps are where a specific extension or compatibility theorem could be valuable. [Horstmeyer–Atay](https://arxiv.org/abs/1607.01237).

Reversible operations may be organized as a groupoid. Irreversible preparations, measurements and coarse-graining need a more general process category or explicit transition kernels. This is bookkeeping for actual operations, not a requirement to adopt category theory as a new physical substance.

Add further structures only for the jobs they perform:

- A quadratic response, if independently justified, produces Gram matrices and the corresponding Clifford algebra.
- A carrier measure produces density/Jacobian laws; it is not inferred from an exponent selected in advance.
- A phase or spin bundle represents specified transport and coverings, with each sign map explicit.
- A statistical realization adds states, effects, probabilities and composition rules.
- A flow or field equation produces dynamics; its principal symbol determines propagation when applicable.
- A physical reconstruction assigns events, rods, clocks, sources, areas and detector readings.

This core is general enough for non-physics applications. The current model becomes one family of realizations inside it. Its low-dimensional coincidences remain useful theorems of those realizations, not universal conclusions.

## 6. Which alternatives are worth pursuing

**First choice: stratified observation geometry plus constrained dynamics.** It directly addresses the present errors and the user's singularity focus. Begin with the signed-frame domain and derive or classify an equivariant vector field preserving its admissibility conditions. Prove existence, uniqueness and access to its boundaries. Add a local hyperbolic PDE when propagation between events is required. This can proceed without first deriving electromagnetism or the fine-structure constant. General PDE theory supplies the relevant initial-value questions; the model must supply the coefficients. [Geroch](https://arxiv.org/abs/gr-qc/9602055).

**For the cone/response problem: premetric constitutive reconstruction.** It permits non-metric responses to survive long enough to fail a physical criterion. Its advantage over starting with \(\star_q\) is that the characteristic geometry can be an output. Its cost is additional constitutive freedom and explicit requirements for causality, stability and channel identification. The existing metric response should be the comparison case.

**For spacetime dynamics: a general coframe/connection, with matter-first gravitational closure as a bounded candidate.** A coframe and connection allow nonflat spatial geometry while preserving the frame viewpoint. They do not themselves choose a gravitational law. If an actual local matter/load equation is obtained, gravitational-closure methods offer a conditional route to compatible geometry actions by requiring matter and geometry to evolve on shared Cauchy surfaces independently of intermediate foliation. They require physicality conditions and solving closure equations; they do not manufacture the matter law or guarantee uniqueness. [Düll and collaborators](https://arxiv.org/abs/1611.08878). A weakly birefringent example retains eleven parameters, illustrating the remaining freedom. [Schneider and collaborators](https://arxiv.org/abs/1708.03870).

**Reserve more elaborate frameworks for a demonstrated need.** Lorentz–Finsler geometry is useful if quadratic response fails; metric-affine geometry is useful if independent torsion or nonmetricity is physically selected; singular symplectic reduction is useful after a Hamiltonian lift exists. Each changes the admissible model class and adds obligations. None is automatically stronger merely because it is broader.

## 7. Where original mathematics could genuinely help

No finite audit can prove that a construction is historically novel or that no future framework could improve it. The sources examined already cover most of the general language. The present evidence supports original **theorems connecting the model's structures**, rather than a new algebra chosen to make the old desired outputs inevitable.

The strongest bounded targets are:

1. **Operational selection of quadratic response.** State laws of ruler composition, scaling and reference change, then prove that they force a bilinear form—or classify the nonquadratic alternatives. Failure is informative: two admissible responses with different cones disprove uniqueness.
2. **Observation and singularity classification.** Classify rank/stabilizer strata and observation-critical sets, with invariant norms and allowed approach classes. Distinguish a divergent instrument response from incomplete physical evolution. The coplanar correction is an immediate prototype.
3. **Compatibility of the different covers.** Construct the maps between signed orientation, Pin lifts, permutation operations and central blocks. Prove when they agree and exhibit the obstruction when they do not. An identification by shared cardinality is insufficient.
4. **Descent of dynamics through a seat.** Determine which candidate flows preserve observation fibres, what extra readings close the reduced law, and what happens when fibre rank changes. The regular case has established machinery; the model-specific boundary problem may require new results.
5. **Reconstruction and closure of the response law.** Classify admissible local laws at an explicitly chosen derivative order and prove whether the physical cone, coupling and evolution are fixed. If free parameters survive, show their effects rather than selecting their measured values and calling them derived.
6. **A singularity-extension or obstruction theorem.** With a specified physical curve and detector class, prove unique continuation and finite physical response through the proposed transition, or identify exactly where it fails. Only after that theorem does a black/white orientation test have a definite meaning.

These are original research opportunities, not promised discoveries. A successful framework might combine established mathematics in a new way and derive a new physical law without creating a new branch of mathematics. Conversely, a new mathematical vocabulary with no stronger observable or closure would not advance the goal.

## 8. Recommended next work

**First, repair the mathematical contracts.** Separate observation from state, real sector from complexification, the distinct sign structures, and coordinate response from intrinsic curvature. Correct the explicit domain/group errors identified above. Preserve successful conditional theorems and retain failed candidates with their exact scopes.

**Second, choose one reconstruction question.** For the stated black-hole focus, the best immediate question is: what minimal evolution law on the admissible signed-frame state space permits a physically interpretable path through its observation and rank boundaries? Keep the law independent of the desired collapse/expansion answer. Work through the definitions and derivation before a new numerical suite.

**Third, establish the physical comparison.** Specify the event carrier, clock, causal propagation, area and detector map. Show precisely which part of the effective Schwarzschild construction is recovered, and where any replacement changes it. A boundary in an internal presentation is not yet a horizon.

**Fourth, earn one uncalibrated consequence.** After the law and remaining parameters are fixed by independent inputs, derive a collapse criterion, a detector response near the transition, or a source-matched rotational observable that was not used to choose the law. A novel prediction may be qualitative or quantitative, but it must distinguish admissible outcomes rather than relabel a known one.

The recommendation is therefore to generalize the model's **method of construction**, while keeping each physical realization sharply specified. Its best asset is the reference-dependent presentation idea. Making its observation and dynamics maps explicit gives that idea a chance to produce physics that the mathematics did not select in advance.

## 9. Verification and limits

The audit combined three independent scoped inspections and a synthesis of the constitutive and generalization questions. Current source files, mathematical definitions and primary literature were inspected. The original THM-E/F runners cited by older material were not present in this checkout; those families were assessed through the available formulas and later receipts, without claiming to reverify their complete original proofs. Exact counterexamples were derived before computation. A short independent SymPy receipt confirmed six decisive examples: the joint zero-lapse/rank-loss path; three hyperbolic pair-spans; finite trine quotient curvature; block exchange without fourth-vector anticommutation; vanishing symplectic quadratic expectations; and the incorrect inverse-cover branch.

The receipt is [audit-counterexamples.py](</mnt/files/My Files/Documents/Codex/2026-09-05/ca/outputs/audit-counterexamples.py>), with [results](</mnt/files/My Files/Documents/Codex/2026-09-05/ca/outputs/audit-counterexamples.json>). Existing PASS counts were not treated as new validation, and no broad model suite was rerun. The report is not a proof of physical novelty, observational viability, or the sufficiency of any proposed replacement. It is a foundation audit with explicit retained results, counterexamples and a recommended direction.
