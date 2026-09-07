# Research dependency ledger

Will's instructions, 5 September 2026: register each online search and mathematical import, extend the native Cella work, and retain strong standard mathematics without ceremonial rederivation.

Maintain `docs/EXTERNAL-MATHEMATICS-DEBT.md` during research:

- Retain standard mathematics that is independently reproducible from its stated first principles; do not repeat a derivation merely to remove an external label. Prioritize Cella extensions that add mathematical content, complete their natural construction beyond a single application, and defer purely external items unless a native connection emerges.
- Before adopting external mathematics, check the relevant Cella DAG claims and read their backing proofs, including applicable DIS derivation notes. The topic library contains symlinks; use canonical source paths or follow the links. Record existing solutions and concrete extensions before declaring a gap. Keyword misses alone do not establish absence.
- Record dated search queries, sources actually read, and whether a dependency was adopted.
- For every adopted external mathematical result or framework, record its source, exact use, assumptions, and whether it is retained standard mathematics or leaves a specific native construction to derive.
- Preserve provenance after rederivation. Distinguish a proof conditional on imported structure from derivation of that structure from the native primitives.
- Do not mark foundational debt closed merely because a symbolic or numerical check passes.
- This ledger supports continuing the authorized research; it does not introduce an approval gate.

> Standing rulings: `docs/RULINGS.md`. Handoff index: `handoffs/README.md`. Runner era index: `suites/README.md`.

## Analogy containment (Will's ruling, 2026-09-07)

**Analogy stays in chat.** Only exact physics labels go into files or runners: headers, input lists, check strings, print
statements, results pages, ledger entries and commit messages name the mathematical object, not the picture that found it.

**Analogy may enter handoffs**, and only there, when explicitly framed as a conjecture or framing tool.

Reason, from the session that produced the rule: the circuit vocabulary was an excellent heuristic — it found the lapse-as-divider,
the matched termination, the negative bulk viscosity and the load/source involution before the mathematics did — and a poor
notation. The notation leaked twice in one day.

1. `prim_dyn1_ohmic_divider.py` was headed "derived from Ohm's law". The spreading resistance 1/a − 1/b *is* the Newtonian
   potential difference; Ohm's law does no mathematical work anywhere in the chain. The actual content is: **the exact lapse is
   linear in the Newtonian potential to all orders**, N² = 1 − φ(r)/φ(r_s) normalised where escape velocity = c. First order is
   automatic for anything matching Newton; the *tail* is the claim, and that is what β = 1 and Mercury's 42.9807″ tested.
2. "Electrode" appears 59 times across seven runners as a name for the inner boundary at r_s. No runner uses an electrode
   property — no injection, source, charge or emf — so no mathematics is contaminated. But it puts the emphasis on the wrong
   terminal and implies the medium is driven from r_s, when the river flows inward and r_s is a sink (and, per EDGE-1, a matched
   termination). The faithful description is a **potentiometer**: r_s one terminal, infinity the other, **the seat is the tap**.

**Vocabulary edit owed** across the DYN/EDGE chain, before the AC sector (where L, C and impedance make naming load-bearing):
r_s → "the inner terminal" / "the termination"; seat → "the tap"; DYN-1's header → Newton, not Ohm. Vocabulary only — re-run every
affected runner to demonstrate that no output changes.
