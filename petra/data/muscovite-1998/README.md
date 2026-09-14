# Sletten & Onstott (1998) comparison data

These CSV files are **manual image digitizations**, not recovered instrument files.
They support qualitative model discrimination only.

## Source

V. W. Sletten and T. C. Onstott (1998), “The Effect of the Instability
of Muscovite During In Vacuo Heating on 40Ar/39Ar Step-Heating Spectra,”
*Geochimica et Cosmochimica Acta* 62, 123–141,
[doi:10.1016/S0016-7037(97)00323-2](https://doi.org/10.1016/S0016-7037(97)00323-2).

The publisher marks the version of record open archive / CC BY-NC-ND 4.0.
The plotted source rasters used for digitization are the publisher's figure
images:

- Figures 3 and 5:
  `https://ars.els-cdn.com/content/image/1-s2.0-S0016703797003232-gr3.gif`
  and `...-gr5.gif`.
- Figure 4:
  `https://ars.els-cdn.com/content/image/1-s2.0-S0016703797003232-gr4.gif`.
- Figures 6–11:
  `https://ars.els-cdn.com/content/image/1-s2.0-S0016703797003232-gr{6..11}.gif`
  (replace the brace expression with the individual figure number).

Downloaded-source SHA-256 receipts (13 September 2026):

- Figure 3: `1b9ed5b7010eba4502d8d09ec3e9b6f27a3dec6d536404f12c8506c80bbbbd31`
- Figure 4: `6bcc205753a98391c6c6706efe32216ceaeee0ec4ac486b1449b77378e225462`
- Figure 5: `b9a3770ff2f30ef552ca9572244d280533000e0609d5ecee899c158e2d58b30c`
- Figure 6: `b007a3788d0d5f4000d056872b63945c69e5e98468e8a78e8702c6b98766be4a`
- Figure 7: `8c090d7aa4e2310bde154508c4111a5c31a6b26939ea3b19e8d183713bebfe50`
- Figure 8: `b79aafd8fb2c3c68aa35ac536ef4a27b3a41ee57a792c173ac57851a8ea436bd`
- Figure 9: `0cf7714444445c5b3e48a999f26716f47d656249279c82b8064c056acad5ebd9`
- Figure 10: `88caaba89104f722fe28cd2e369b055602d219c537a0358cc852fb0a7f62f78b`
- Figure 11: `c58f3c7472a9aae6562bea6a49c27d16683c45e446173d2435b820132a2899ce`

The raster images are not redistributed here.  The CSVs are new point estimates
read from their plotted axes.

## Method and limits

Values were read at visible step boundaries, markers, and curve extrema.  The
`estimated_uncertainty_*` columns are conservative plot-reading uncertainty,
not analytical 1-sigma errors. Grain-size ranges are represented by their
midpoint (`150–180 µm` → `165`; `52–64 µm` → `58`). Curves are intentionally
sparse; interpolation must not be presented as recovered raw measurements.

Figure 4 contributes the square-marker 40Ar* release-rate series. Figure 3
contributes the data staircases for the 40 h / 500 °C, 64 h / 700 °C, and
45 d / 700 °C / 2 kbar treatments. Figure 5 contributes the magnified
low-release-fraction steps used to test the anomalously old initial-step
signature. `figures6-11-time-series.csv` adds sparse cumulative-loss and
cylinder-D/a² points for the 500 and 700 °C isothermal runs (Figures 6, 7, 9,
10), the Figure 8 dehydroxylation/Ar comparison, and the Figure 11a
³⁶Ar/⁴⁰Ar diffusivity-ratio traces. The Figure 6/7/9/10 lower panels contain
several geometry inversions; E4b digitizes and overlays the paper's labelled
cylinder branch only because the model analysis deliberately applies the same
infinite-cylinder inversion.

The original experiment's exact incremental dwell-time table was not present in
the A8 thesis archive. E4 therefore uses an explicit **comparison schedule** of
500–1200 °C in 50 °C, 600 s increments: it matches the Figure 4 temperature
support and preserves E2's 600 s dwell convention, but it is not represented as
the original furnace controller log. Any quantitative fit requires recovery of
the original GD150/source data or the full methods table.

## E3b calibration decision fragment

`e3b-calibration-fragment.toml` records all three coordinate-matched E3b
`disagrees` verdicts and the value-selection decision for E4. The PBE-D3 and
classical numbers are frozen-path rises, not relaxed barriers, so E4 applies
none of them as corrections. Its existing local-dehydroxylate rule retains the
64.095991 kcal/mol E3a incomplete-convergence value strictly as a bounded
sensitivity; the reconstructed-replication and Xe E3b numbers are not consumed
by the current Ar-only E4 decks.
