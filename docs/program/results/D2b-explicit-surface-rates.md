# D2b — explicit-surface LH rates: results and verdict

*Campaign run 2026-08-27 on the 4090 workstation by fable
(hermes-custom-build-001; profile=workstation). Driver:
`qm/scripts/surface_rate_protocol.py` (final commit `87ea5eb`). Full
artifacts (geometries, scans, CC receipts, rate tables):
`qm/runs/D2b-explicit-surface-rates/` in the campaign worktree,
gitignored per big-output policy; logs `qm/runs/D2b-campaign{,-stdout}.log`.
Two passes: 5/10 clean, then four driver-hardening fixes from live
failures and a catch-up pass → **10/10 reactions, 0 failures**.*

## The protocol (three rungs, each with receipts)

D2a compared quarry's gas-phase supermolecule rates against Song &
Kästner's explicit-ASW instanton fits and found the addition channel
4–5 orders low — the explicit NO-GO this card resolves. The protocol:

- **Rung A — surface partition functions** (Meisner/Lamberts/Kästner
  implicit-surface convention): unimolecular Eyring rates from the
  adsorbed pre-reactive complex, vibration-only partition functions
  (`quarry.rates.surface_thermo_from_frequencies`). Measured effect:
  **factor ~2** on the addition channel. Not the gap.
- **Rung B — coupled-cluster barriers**: in-driver open-shell
  UHF-UCCSD(T) (A2's engines are closed-shell only); CBS(cc-pVTZ,QZ)
  on gas-phase systems, direct TZ on 1-water clusters, additive
  gas-phase ΔΔ transfer to 2-water clusters. Corrections: +1.8 kJ/mol
  (H+CO), −0.1 (addition), **+3.3 (abstraction)** — the direction and
  channel asymmetry Simons+20 documented against S&K's DFT rates.
- **Rung C — explicit clusters + site sampling**: every channel on
  (H₂O)₁ and (H₂O)₂ in distinct binding orientations (3 wet sites for
  H+CO, 2 each for the H2CO channels), substrate atom order invariant.

Anchors, with the caveats the field itself documents: S&K 2017 fits
(explicit-ASW QM/MM instanton, valid ≥59 K, *documented overestimates* —
Simons+20 §3.2), Simons+20 Table 2 adopted plateaus (2·10⁵ s⁻¹ H+CO,
1.5·10⁵ s⁻¹ per H2CO channel — the values a KMC benchmarked against the
Chuang/Watanabe experiments actually uses), Fuchs+09 experimental
effective barriers (12–16.5 K, compared only in their own ν=2·10¹¹ s⁻¹
convention), Andersson+11 T_c = 79 K for H+CO.

## What the numbers say

**Barriers (the validated tier).** ZPE-corrected, CC level:

| Channel | gas | 1w | 2w | literature |
|---|---|---|---|---|
| H + CO → HCO | 1692 K | 1900 / 1994 K | 1697 K | 1700–2500 K effective band |
| H + H2CO → CH3O | 2260 K | 2470 K | 2017 K | S&K DFT 1900 K (documented low) |
| H + H2CO → H2+HCO | 2640 K | **3028 K** | 2977 K | **S&K surface 3030 K** |

The 1-water abstraction barrier reproduces S&K's published surface
value to 2 K. Water raises both H2CO barriers at 1w (H-bonded TS) and
relaxes at the 2w site — a real ~450 K site spread, receipts for why
site sampling is mandatory. Worst UHF ⟨S²⟩ 0.889 (addition TS; 0.75
exact) — moderate contamination, noted; UCCSD(T) is tolerant but this
channel's reference is the shakiest.

**Rates vs anchors** (k_surface_cc, s⁻¹; * = Eckart extrapolation
below the validity floor, reported only for anchor comparison):

| Channel (site) | 60 K | vs S&K fit | 15 K* | vs Simons plateau |
|---|---|---|---|---|
| H+CO (gas) | 1.8·10⁵ | — | 3.7·10³ | 54× low |
| H+CO (1w sites) | 1.5–2.5·10⁴ | — | 3.0–5.5·10² | ~500× low |
| CH3O (gas) | 5.6·10⁴ | 3.7·10⁻⁵ | 2.8·10³ | 54× low |
| CH3O (1w / 2w) | 3.9·10³ / 2.2·10⁵ | 2.5·10⁻⁶ / 1.4·10⁻⁴ | 1.3·10² / 6.3·10³ | 1100× / 24× low |
| H2+HCO (gas) | 8.2·10⁶ | 22× | 1.1·10⁶ | **7× high** |
| H2+HCO (1w / 2w) | 4.5·10⁶ | 12–14× | 3.7–4.5·10⁵ | **within 2.5–3×** |

Fuchs-convention effective barriers from our extrapolated rates run
150–250 K below their experimental fits (H+CO: 218–305 K vs 390–520 K)
— our extrapolated plateau rates are *fast* in their convention while
being *slow* against the intrinsic Simons values; the two conventions
bracket, as expected for effective-vs-intrinsic comparisons.

**The decomposition of D2a's 4–5 orders**, now itemized: ~2× partition
functions (rung A) + 0–3.3 kJ/mol barrier error (rung B; worth up to
~10² at 60 K) + **the S&K fits being documented overestimates** (up to
~10³·⁵ of the "gap" was in the anchor, not the calculation) + the
remainder, 1.5–3 orders in deep tunneling, is **Eckart-vs-instanton
barrier-shape error** — confirmed independently by the crossover
mismatch on H+CO: our PWB6K imaginary mode (645–678i) gives T_c ≈
148–155 K vs Andersson's 79 K on a CC-fitted PES; PWB6K's barrier is
too stiff, and Eckart inherits the wrong width exactly where width is
everything.

## Verdict (the card's gate)

- **NO-GO — the D3b rate table is not emitted.** At 10–20 K the
  addition and H+CO channels sit 25–1600× below the corrected
  literature anchor with a ~60× site spread, and the implied
  abstraction:addition branching (~70:1 at 15 K) contradicts the
  near-equal branching the experiment-benchmarked network requires.
  Feeding these into the D3-lineage deck would poison it. Per the
  card's acceptance, no machine-readable table ships.
- **GO — the protocol tier itself.** Surface partition functions,
  open-shell CC barrier machinery, tiered wet-site gates, and
  orientation sampling are validated and reusable: the barrier tier
  reproduces the field's benchmark surface barrier to 2 K, and the
  abstraction channel's full k(T) lands inside the documented spread
  at every anchor temperature.
- **Validity floor, stated:** with CC barriers, rates at T ≥ 100 K
  carry ≲1 order of method uncertainty; 59–100 K ≲2 orders
  (channel-dependent, abstraction best); below 59 K asymmetric-Eckart
  is systematically low by 1.5–3 orders against instanton-class
  anchors and must not be published as absolute surface rates.
- **The missing piece is named and bounded:** instanton or
  small-curvature tunneling along the existing IRC machinery — the
  scoping doc's stretch goal — is the single upgrade that flips this
  gate. Carded as D2c-instanton-tier.

## Method caveats

Geometries/frequencies at pwb6k/def2-SVP/D3BJ/DF; CC single points
frozen-core UHF-UCCSD(T) (no F12); 2w barriers carry gas-phase ΔΔ
transfer (direct-TZ transferability check at 1w agrees to ~1 kJ/mol
for abstraction, ~2.7 for addition — the transfer is the weakest rung-B
link). Wet-site gates tolerate spectator imaginary modes <100 cm⁻¹
(disclosed per-structure in `frequency_data`); the chemical mode gate
(≥200 cm⁻¹, exactly one) is unchanged from D2a. Site sampling covers
binding orientation, not an ASW ensemble — a Grassi-style distribution
over an amorphous slab remains future (Phase-3) work.
