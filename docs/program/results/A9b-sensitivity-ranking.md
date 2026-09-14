# A9b sensitivity ranking — corrected approximate-rate ensemble

## Verdict

The independently replayed survey-tier campaign is complete, reproducible, and
scientifically negative. All 232 required trajectories (29 scenarios × 8 fixed
seeds) completed at 298 K and pH 4; every scenario passed the full per-kind
state-distribution stationarity gate. Nevertheless, no trajectory observed a
biased original-lattice Si or Al release. The unnormalized importance-sampling
estimator therefore returns typed
`censored_zero_biased_observations`, not zero rates and not finite upper bounds.

Consequently:

- original-lattice Si rate: **unavailable / censored**, mol m^-2 s^-1;
- original-lattice Al rate: **unavailable / censored**, mol m^-2 s^-1;
- Si:Al stoichiometry: **undefined**;
- ordinal sensitivity ranking: **not scientifically accepted**;
- rankable barrier families: **none**;
- top-k survey targets: **none**;
- families demonstrated irrelevant at this tier: **none**.

Undefined responses remain unranked. This is a survey-tier platform-test result,
not a production calculation.

## Design and evidence

The campaign applies the merged reachability/biased-sampling implementation and
the hash-bound pH 3–5 open-flow reservoir/origin contract. It evaluates nominal
plus ±3 kcal mol^-1 barrier and ×0.1/×10 prefactor perturbations for all seven
families:

1. `siloxane-neutral`
2. `sioal-si-neutral`
3. `sioal-al-neutral`
4. `connectivity-ladder`
5. `al-o-al-analogue`
6. `adsorption`
7. `cation-desorption`

Every replica preserves its command receipt, complete event stream, likelihood
increments, checkpoints, PGIF lineage, stopping evidence, and typed censoring.
The bundle contains 232 receipts and 232 each of event, likelihood, and checkpoint
streams. The public launcher admitted the run through the fixed bounded systemd
unit with 4 workers × 4 threads, niceness 10, and `RuntimeMaxSec=86400`; observed
runtime was 19.101 s, peak memory 106.0 MiB, and no swap.

Durable local evidence (not committed because the bundle is 68 MiB):

- bundle: `/mnt/data/vsletten/run-outputs/A9b-sensitivity-ranking-20260914`
- operator receipt: `/mnt/data/vsletten/run-outputs/A9b-sensitivity-ranking-20260914-operator/launch-attempt-1.json`
- systemd invocation: `003188e2e63545af8dd2111887353089`
- manifest SHA-256: `b224a5a48729024d41c70a87cb185fb0c875496df4209f40c127a689d1d07ac2`
- preliminary analysis SHA-256: `7830c5eb46ccd12dca111c686a01a79124c7adb1a37234b18aab987ead236fb2`
- independent verification SHA-256: `66fdbde89522a00e612007d924ff3bed489bc15c9bd84267fd1f14d2a264b04d`

The independent verifier regenerated all 29 scenario decks, replayed all 232 raw
replicas, rehashed source inputs, recomputed likelihood increments, conversions,
state stationarity, lineage, estimator outputs, censoring, and ranking, and then
proved its derived bytes equal the analyzer's finalized result.

## Laboratory comparison

The contract's Palandri–Kharaka acid-plus-neutral 298 K approximation gives
kaolinite formula-unit rates of:

| pH | laboratory rate (mol m^-2 s^-1) | simulated gap |
|---:|---:|---|
| 3 | 8.892533283451393e-14 | unavailable; comparison bound only |
| 4 | 6.988878750916011e-14 | unavailable; simulated rate censored |
| 5 | 6.670760828695020e-14 | unavailable; comparison bound only |

The pH 3 and pH 5 values are laboratory-comparison bounds only; only pH 4 was
simulated. Because both simulated species rates are censored, a kaolinite
formula-unit rate and log10 gaps cannot be formed without fabricating information.

## Reproduction

```bash
petra/scripts/launch_a9b_sensitivity_campaign.sh \
  /mnt/data/vsletten/run-outputs/A9b-sensitivity-ranking-20260914 4
```

The launcher self-dispatches to the bounded transient unit when called outside it.
The verifier is the final authority; a successful process exit never upgrades a
censored response into a numerical rate or rank.
