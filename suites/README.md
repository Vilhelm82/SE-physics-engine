# suites/ — live work and scaffold

**The problem this file exists to fix.** Cl(3)-era runners sit beside the post-pivot rebuild with nothing distinguishing
them. A session reading `thm_g.py` finds "the hbar-seat pivot is a SQUEEZE" and concludes a seat can pivot, which
`docs/RULINGS.md` R1 forbids. Pre-pivot files are **provenance**: they record what was declared at the time and must
**not** be edited to match current rulings — they are superseded, not rewritten.

**Classification is Will's.** Below is a machine pass on the file headers only, offered as a starting point, not an
answer. The test used: a header that lists Cl(3) or BARE-1 as an INPUT/DECLARED/GIVEN is scaffold; a header that lists
it as BANNED is post-pivot rebuild; anything else is unclassified.

## Confirmed scaffold — do not cite as live inputs
- `thm_g.py`, `thm_g2.py` — both assert an "X-seat pivot"; contradicts RULINGS R1. Left as written (provenance).
- Any runner whose header takes Cl(3), BARE-1 or the elliptope as an input.

## Machine pass

**Declares Cl(3)/BARE-1 as an input (2):**
- `e8x_cross_pairing.py`
- `sect1_one_norm.py`

**Declares Cl(3)/BARE-1 BANNED — post-pivot rebuild (5):**
- `prim_dyn3_stokes_swirl.py`
- `prim_t1_t3_pivot_group.py`
- `prim_t4_hawking_period.py`
- `prim_t7_seat_form.py`
- `prim_t7b_labelling.py`

**Unclassified by this test (89)** — needs Will, or a header that says which era it belongs to:
- `bargmann.py`
- `bnd1_boundary_receipt.py`
- `branched_resolved_response.py`
- `bridge_dbp.py`
- `cayley.py`
- `cella_constitutive_divisors.py`
- `cella_frame_transport.py`
- `cella_horizon_reuse.py`
- `cella_normal_geometry.py`
- `cella_predictive_state.py`
- `cella_tensor_valuation.py`
- `census_c.py`
- `correspondence_july.py`
- `curv1_path1.py`
- `curv1_path23.py`
- `debt2.py`
- `debt2b.py`
- `doors.py`
- `galois.py`
- `graph_cocycle.py`
- `hg1_hbar_G_is_the_sheet.py`
- `horizon_crossing_metric.py`
- `hunch_z0_impedance.py`
- `kahler.py`
- `kin1_cost_fork.py`
- `label1_wick_loops.py`
- `label2_octahedron.py`
- `label3_lattice.py`
- `label3b_corrections.py`
- `oc1_open_circuit.py`
- `pinch1_horizon_branch_meet.py`
- `pinch2_deck_and_blowup.py`
- `pinch3_cella_valuation.py`
- `pinch_directional_resolution.py`
- `pinch_load_transport.py`
- `pivot_map.py`
- `pivot_map_closed_forms.py`
- `pred1_cross_seat.py`
- `pred1_deck_separator.py`
- `pred1_operational_chain.py`
- `pred1_physical_protocol.py`
- `pred1_protocol.py`
- `prim_dyn1_lapse_from_potential.py`
- `prim_dyn2_horizon_period.py`
- `prim_edge1_horizon_absorption.py`
- `prim_edge2_tensor_love.py`
- `prim_edge3_rotating_dissipation.py`
- `prim_edge4_tidal_friction.py`
- `prim_edge5_spin_ratio.py`
- `prim_edge6_membrane_viscosity.py`
- `prim_eq1_tolman_from_lapse.py`
- `prim_fold1_lightring_dressing.py`
- `prim_gud1_complement.py`
- `prim_gud2_inversion_map.py`
- `prim_ns0_river_euler.py`
- `prim_ns1_swirl_inertia.py`
- `prim_sheet1_block_and_holonomy.py`
- `prim_t5_area_and_lifts.py`
- `prim_t5b_parity.py`
- `prim_t5c_corrections.py`
- `prim_t7d_tilt.py`
- `prim_t7e_two_degeneracies.py`
- `prim_t7f_regions.py`
- `prim_t7g_two_horizons.py`
- `prim_t8a_fourth_direction.py`
- `prim_t8b_adjoint.py`
- `prim_t8c_reciprocity.py`
- `quantum_prediction_probes.py`
- `rest_frame_trapping_equivalence.py`
- `thm_b_monodromy.py`
- `thm_d2_unruh.py`
- `thm_g.py`
- `thm_g2.py`
- `thm_galois_deck_descent.py`
- `thm_h.py`
- `thm_h2_d1.py`
- `thm_i_field.py`
- `thm_i_pre.py`
- `thm_i_transport.py`
- `thm_j_dyn.py`
- `thm_k_clausius.py`
- `thm_k_response_map.py`
- `thm_l_rotation.py`
- `thm_m_swirl.py`
- `thm_n_kerr_quadrupole.py`
- `thm_o_strain_law.py`
- `thm_rn.py`
- `view1_incidence.py`
- `xlink1_cella_coupling_form.py`

## Rule going forward
Every new runner header states its era and its inputs. A runner that does not say what it takes as given cannot be
cited as a live input by anything else.
