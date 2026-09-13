# A3j-calc005-si-n1-verification — evidence-only ruling that A3i's receipt stands

- status: done
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: A3i-calc005-si-n1-pilot
- claimed-by:

## Objective

Record Victor's 2026-09-13 mission-control TASK-298 evidence-only ruling that
A3i's hash-bound terminal-failure receipt stands, and that this card does not
perform the independent optimizer-free verifier originally specified for the
CALC-005 Si `i=1`, `(x=1,y=0)`, `Osa.sih -> Osa.albr` pilot. This closeout is a
separately authorized board ruling, not a typed `verified-pass` /
`verified-fail` receipt. It releases parent A3 from this gate; it does not
re-establish component identity, stationary minima, thermodynamic arithmetic,
Store provenance, or the no-emission boundary.

## Constraints

- Follow fleet `POLICY.md` and `docs/program/PROTOCOL.md`.
- Do not treat this ruling as a substitute for the original hash, artifact,
  energy, gradient, thermochemical, or no-emission checks, or as a typed
  verification receipt.
- Do not launch the original verifier, replay A3i, retry, `n=2..4`, Petra
  emission, Store republish, or CALCULATIONS update from this card.
- A3i's receipt remains the standing hash-bound terminal-failure evidence; this
  card does not overwrite or re-open it.

## Acceptance

- Card status is `done` by the separately authorized evidence-only ruling, not
  by verifier PASS.
- Objective, Constraints, and this Acceptance section document that the original
  independent-verifier campaign and typed verification receipt are not the
  done-gate.
- No `verified-pass` / `verified-fail` receipt is required or produced to mark
  this card done.
- Parent A3 is released from the A3j gate by that ruling; the retired verifier
  campaign is not restarted here.

## Result

Closed by ruling 2026-09-13 (mission-control TASK-298); A3i's hash-bound
terminal-failure receipt stands. The original independent verifier was not run
and no typed verification receipt was emitted. (hermes-custom-build-001;
profile=workstation)
