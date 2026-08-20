# D2a-astro-rate-reproduction — reproduce 5 known astrochemical rates

- status: done
- track: D (astrochemistry)
- priority: P1
- machine: workstation (GPU preferred; systems are small — CPU fallback
  acceptable for a subset)
- depends: —
- claimed-by: (see agents/D2a-astro-rate-reproduction branch)

## Objective

The de-risking gate for the whole astro track: before computing any
*new* rate, reproduce ~5 published tunneling-dominated reaction rates
with quarry and quantify agreement. Target set (final pick justified in
the run log, from docs/scoping/astrochemistry.md §validation): gas-phase
H + H2CO → HCO + H2; H abstraction OH + H2 → H2O + H; one or two steps
of the CO hydrogenation ladder (H + CO → HCO; H + H2CO → CH3O) against
published instanton/experimental values; one barrierless control (should
give κ≈1 sanity). For each: barrier + imaginary frequency + Eckart κ(T)
over 50–300 K, tabulated k(T) vs the literature values, with the
comparison honest about where Eckart degrades (expected: low-T
underestimate vs instanton — characterize, don't hide).

Deliverable: `qm/runs/` outputs + a summary doc
`docs/program/results/D2a-rate-reproduction.md` (tables + verdict:
is quarry's 50–300 K accuracy sufficient for KIDA-grade submissions?)
+ any quarry fixes the campaign forces (each its own tested commit).

## Context

- docs/scoping/astrochemistry.md (validation anchors, Eckart caveat)
- qm/quarry/{ts,rates,pipeline}.py; scripts/phase1_xiao_lasaga.py is
  the campaign-driver pattern to follow (checkpointing, tee-to-log,
  structural gates)
- .claude/learnings/debugging.md — the saddle-search traps apply
  verbatim to these reactions

## Acceptance

- Summary doc on main with ≥5 reactions, each: our barrier vs
  literature, our k(T) at 3+ temperatures vs literature, deviation
  stated; a clear go/no-go verdict for Phase D2b (new rates).
- All campaign logs teed to files; results carry provenance.
- Compute etiquette respected (GPU-first, thread caps).

## Progress

- 2026-08-19 — local Hermes worker (gpt-5.6 sol) claimed via branch push
  — the protocol's first non-Claude runtime. Built the checkpointed
  campaign driver + open-shell GPU Hessian fallback + first-crest
  seeding + per-reaction failure isolation (commits 689fbed..b03d599);
  ran the full 6-reaction campaign on the GPU (etiquette followed,
  logs teed). Hit the tool-iteration ceiling during closeout.
- 2026-08-20 — fable — closeout: summary doc written from the durable
  artifacts, verdict recorded. DEVIATION: card scoped "~5 reproduced
  rates"; delivered 4 completed + 2 honestly-rejected-with-receipts —
  for a de-risking gate the rejections are results.

## Result

Verdict: **GO for gas-phase rate classes, NO-GO for surface (LH) rates
from the cheap gas-phase protocol** (4–5 orders low vs surface
literature for H+H2CO→CH3O; abstraction channel within 3–50× across
50–300 K; barrierless control κ=1 exact). Full analysis:
docs/program/results/D2a-rate-reproduction.md. Artifacts:
qm/runs/D2a-astro-rate-reproduction/ (workstation). Follow-ups worth
carding: OH+H2 TS guess strategy; CH2OH constrained approach; explicit-
surface machinery before any KIDA surface submission.

- 2026-08-18 — card created (fable).
