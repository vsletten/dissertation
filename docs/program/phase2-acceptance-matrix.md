# Phase-2 scientific-output acceptance matrix

**Standing rule:** change this matrix first, then make the driver and tests match.

This contract governs every Phase-2 ladder route before thermochemistry or
`results.json`/`store.sqlite` publication. Optimizer completion, geometry
identity, minimum character, saddle index, and reaction connectivity are
separate gates; none substitutes for another.

## Inputs and precedence

1. Revoke prior canonical outputs before a new live attempt.
2. Accept a cached production complex only when its input geometry, atom order,
   electronic state, frozen shell, production settings, endpoint hash, and
   chemistry-gate receipt all match.
3. After production optimization, reapply the same identity, state, shell,
   collision, and exact oxygen-proton-owner gates before checkpointing.
4. Reject the reactant if any PHVA imaginary mode exceeds 30 cm^-1.
5. Reject the saddle unless exactly one PHVA imaginary mode exceeds 30 cm^-1.
6. Reject quick-IRC unless both endpoints have unambiguous oxygen-proton owners,
   the endpoints have distinct connectivity signatures, one is hydrolyzed
   product `(M-Ow=True, M-Obr=False, H-Obr=True)`, and the other retains the
   attacked `M-Obr` bond.
7. Publish the store first and atomically publish `results.json` last.

The connectivity signature uses the same Phase-1 covalent cutoffs: M-O < 2.3 A
and O-H < 1.25 A. The non-product endpoint may be the separated reactant or the
verified associative intermediate because the pilot established a sequential
mechanism.

## Quick-IRC signature matrix

| M-Ow | M-Obr | H-Obr | Role / verdict |
|---|---|---|---|
| false | false | false | reject: no attacked bridge |
| false | false | true | reject: protonated but non-associative fragment |
| false | true | false | accepted non-product: separated reactant |
| false | true | true | accepted non-product: protonated bridge intermediate |
| true | false | false | reject: cleaved bridge without transferred proton |
| true | false | true | required hydrolyzed product |
| true | true | false | accepted non-product: associative reactant |
| true | true | true | accepted non-product: associative proton-transfer intermediate |

A valid endpoint pair contains exactly two distinct rows: the required product
row and any accepted non-product row. Endpoint order is irrelevant. Any
unassigned or ambiguous proton rejects the pair before this table is evaluated.

## Test anchors

- `test_existing_reactant_complex_without_receipt_is_rejected`
- `test_verified_reactant_complex_skips_both_optimization_rungs`
- `test_production_endpoint_rejects_changed_proton_owner`
- `test_load_xyz_rejects_changed_atom_order`
- `test_reactant_minimum_rejects_significant_imaginary_mode`
- `test_quick_irc_must_span_intact_bridge_and_hydrolyzed_product`
- `test_new_attempt_quarantines_stale_canonical_outputs`

## Explicitly deferred to A3a

A3a owns persistence of a nonconverged production endpoint, one and only one
fresh-optimizer continuation, and an independent final production-gradient
receipt. Until that card passes, A3 remains blocked and no fourth Osa-neutral
n=1 calculation or n=2..4 launch is authorized.
