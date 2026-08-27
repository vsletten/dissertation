# E4-1998-comparison — synthetic spectra over the founding data

- status: blocked
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (CPU ensembles; needs the E3a barrier set)
- depends: E3a-classical-neb-barriers, E2b-grain-size-sweep, A8-thesis-archive-intake
- blocked-on: E3a-classical-neb-barriers, E2b-grain-size-sweep, A8-thesis-archive-intake
- claimed-by:

## Objective

The point of the track: integrate E3a's computed barriers (with E3b
calibration verdicts where available) into the full-mechanism deck, run the
laboratory T(t) schedules from Sletten & Onstott (1998) at E2b-validated
volumes, and lay the synthetic ⁴⁰Ar/³⁹Ar age spectra, cumulative release
curves, and ³⁶Ar/⁴⁰Ar correlations directly over the 1998 target figures
(Figs. 3–11) and Villa's isotope-correlation test. Digitized 1998 release
data comes from the thesis archive intake (A8).

Testable claims on the line (scoping doc §5): staircase spectra from a
homogeneous-age lattice with defect + alteration reservoirs (no fossil
⁴⁰Ar* gradient); the rise-then-fall D/a²(t) signature; species-resolved
reservoir release ordering.

## Constraints

- Every plotted synthetic curve states its barrier provenance (computed vs
  calibrated vs residual proxy) — the comparison is only as honest as its
  labels.
- No fitting knobs beyond the physically carded parameters; if the spectra
  don't match, that result is reported as-is (NO-GO is a finding).
- This card produces the paper-1 figure set; manuscript work is a separate
  future card.

## Acceptance

- Synthetic-vs-1998 overlay figures for the untreated and hydrothermally
  treated splits, with ensemble bands, landed under
  `docs/program/results/` with the generating decks/schedules committed.
- An explicit verdict per §5 claim: reproduced / partially / not — with the
  mechanism attribution for each.
- Petra tests, replay gates, lint green; card/PLAN/STATUS bookkeeping.

## Progress

- 2026-08-27 — filed by Fable. Blocked on the computed-barrier campaign
  (E3a), the statistics ladder (E2b), and thesis-archive intake (A8),
  which supplies the digitized 1998 release data.
