# Working instructions

## Research dependencies

Maintain `docs/EXTERNAL-MATHEMATICS-DEBT.md` as you work.

- Retain standard mathematics that is independently reproducible from its stated first principles. A derivation is worth
  repeating when it removes a dependency, not when it removes a label.
- Before adopting external mathematics, check the Cella DAG claims and read their backing proofs, including applicable DIS
  derivation notes. The topic library uses symlinks; follow them or use canonical source paths. Record existing solutions
  and concrete extensions first; a keyword miss establishes nothing.
- Record dated search queries, the sources actually read, and whether a dependency was adopted.
- For every adopted result, record its source, exact use, assumptions, and whether it is retained standard mathematics or
  leaves a specific native construction to derive.
- Preserve provenance after rederivation: state clearly whether a proof is conditional on imported structure, or derives
  that structure from the native primitives.
- Prioritise Cella extensions that add mathematical content and complete their natural construction beyond a single
  application. Defer purely external items until a native connection emerges.
- Foundational debt closes when the construction is native, not when a check passes.

## Naming

Files and runners carry exact physics labels — headers, input lists, check strings, print statements, results pages,
ledger entries and commit messages name the mathematical object. State what a quantity *is*: `N² = 1 − φ(r)/φ(r_s)`,
`the sonic surface r = r_s`, `the involution η → 1/η`.

Analogy belongs in conversation, and in handoffs when framed as a conjecture or framing tool.

## Results

- Verify before ruling: read the file, run the computation. Two independent paths where it matters.
- Label every claim — proved / derived-given-X / retrodiction / conjecture / hunch — and name X out loud.
- Give every result a kill condition. If none can be stated, say so; that is information.
- Commit computational results with their logs. Commit messages carry check counts and what was held out.
