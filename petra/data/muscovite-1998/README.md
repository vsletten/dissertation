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

Downloaded-source SHA-256 receipts (13 September 2026):

- Figure 3: `1b9ed5b7010eba4502d8d09ec3e9b6f27a3dec6d536404f12c8506c80bbbbd31`
- Figure 4: `6bcc205753a98391c6c6706efe32216ceaeee0ec4ac486b1449b77378e225462`
- Figure 5: `b9a3770ff2f30ef552ca9572244d280533000e0609d5ecee899c158e2d58b30c`

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
signature.

The original experiment's exact incremental dwell-time table was not present in
the A8 thesis archive. E4 therefore uses an explicit **comparison schedule** of
500–1200 °C in 50 °C, 600 s increments: it matches the Figure 4 temperature
support and preserves E2's 600 s dwell convention, but it is not represented as
the original furnace controller log. Any quantitative fit requires recovery of
the original GD150/source data or the full methods table.
