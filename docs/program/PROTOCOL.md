# Program worker protocol

*How any agent — Fable, Hermes cloud workers, macbots, or a human —
contributes to the Omnibus program from a cold start, using only this
repo. Design goals: zero standing orchestrator, git-native atomicity,
bookkeeping that cannot drift from reality.*

Read in order: [`docs/OMNIBUS.md`](../OMNIBUS.md) (the program),
[`PLAN.md`](PLAN.md) (the board), this file (the rules). Track-specific
background lives in [`docs/scoping/`](../scoping/).

## The three invariants

1. **Locks live in the ref namespace; state lives in files.** Claiming
   is a branch push (atomic, first-wins). Card status, results, and
   program status are files that change only through merged PRs.
2. **One PR carries work + bookkeeping together.** The same PR that
   implements a card also updates the card and STATUS.md — merged means
   both happened; there is no separate reporting step to forget.
3. **Everything is idempotent.** Before working a card, verify its
   acceptance criteria aren't already satisfied by repo state. A card
   is done iff its acceptance checks pass on main — not because someone
   said so.

## Claiming a card

```bash
# 1. Pick: lowest-priority-number READY card in PLAN.md whose depends
#    are all done and that your machine can run (see machine classes).
# 2. Check it isn't claimed (a live agents/<card-id> branch = claimed):
git ls-remote --heads origin 'agents/*'
# 3. Claim — atomic; if this push is rejected, someone beat you; pick another:
git push origin origin/main:refs/heads/agents/<card-id>
# 4. Work in a worktree on that branch (never on main, never switch an
#    existing worktree's branch):
git -C <repo>/main worktree add <repo>/agents/<card-id> agents/<card-id>
```

The claim branch **is** the lock. Repo setting "auto-delete head
branches" releases it automatically when your PR merges. If you must
abandon, delete the branch so the card frees.

**Stale claims**: a claim branch with no open PR and no commits beyond
main, older than 48h (or whose PR is closed-unmerged and idle) is
abandoned — delete it, note the reclaim in your PR description, claim
normally.

## Working a card

- Do exactly the card's Objective; its Acceptance section is the
  contract. Scope creep goes into a *new card*, not your PR.
- Repo gates (all required): tests green (`uv run pytest` in qm/,
  `cargo test` in petra//kmc-rs as applicable), lint clean (ruff /
  clippy), and the specific gates the card names. Petra engine changes
  additionally require the bitwise-parity gates.
- **Deviations are normal — hiding them is not.** If reality forced a
  different approach than the card describes, do the right thing and
  record it in the card's Progress log with a `DEVIATION:` prefix.
- Non-obvious discoveries go to `.claude/learnings/` (existing
  convention: topic files + INDEX.md one-liner).

## Bookkeeping is incremental, not terminal (amendment, 2026-08-20)

Learned from D2a: the worker did excellent science, then hit its
tool-iteration ceiling during the closeout — summary, card update, and
PR all unwritten (supervision finished them from artifacts). Rules:

- **Push durable state as you go.** Implementation commits go up when
  green, not at the end. Campaign artifacts checkpoint continuously.
- **Update the card's Progress log at each stage boundary** (committed
  and pushed), so a dead or capped worker leaves a legible trail.
- **Budget the closeout.** Your *iteration budget* is the worker's
  tool-iteration ceiling (Hermes worker-profile `max_turns`); every
  tool call — model turn, file read/write, shell command — spends one
  iteration. Track turns used against the ceiling and, when ~15% of
  the budget remains, STOP new work and execute the closeout (summary
  → card → STATUS → PR). An incomplete card with finished bookkeeping
  beats a finished card with no bookkeeping.
- Exit reports matter: state precisely what is and is not done —
  D2a's honest exit made supervision closeout a 15-minute job.

## Finishing a card — the single PR

Your PR contains, atomically:

1. The work itself.
2. The card file updated: `status: done`, a dated entry in its
   `## Progress` log, and a `## Result` section (what/where/surprises —
   enough for the next card's worker).
3. One line appended to [`STATUS.md`](STATUS.md) (newest first).
4. Any new follow-up cards (status: ready or blocked) your work exposed.

Open the PR (`gh pr create`) — **then STOP** (amended 2026-08-20 to
match fleet POLICY §1 v7, which supersedes the old wording here): this
repo is webhook-wired, and the cloud Hermes PR watcher owns the rest
of the lifecycle — review fixes (including Sourcery threads), CI
repairs, branch updates against main, and the squash-merge. Do NOT
enable auto-merge, reply to or resolve review threads, poll checks, or
merge locally; a worker racing the watcher corrupts both. To hold a PR
for human/Fable review instead, add the `human-merge` label. A PR
still open after ~48h means the watcher has stalled: flag it (the
daily supervisor also watches for this). Background on the underlying
gates (four required checks, strict up-to-date, required conversation
resolution): `.claude/learnings/` and the memory of PRs #12–#21.

## Machine classes

| Class | Can run | Cannot run |
|---|---|---|
| **workstation** (4090 box) | everything, incl. GPU QM campaigns | — |
| **cloud / macbot** | design docs, Rust/TS code + tests, CPU-only QM tests (tiny molecules), scoping/research | GPU campaigns, multi-hour CPU jobs |

Cards carry a `machine:` field. Compute etiquette on the workstation is
non-negotiable: GPU-first, `OMP_NUM_THREADS≤16` for CPU stages, every
long job teed to a log file, no multi-hour full-CPU jobs without
Victor's explicit go.

## The inbox — async agent-to-agent mail

`docs/program/inbox/` holds dated markdown messages:
`YYYY-MM-DD-<from>-to-<audience>-<slug>.md`. Delivery is by PR like
everything else (fold messages into your card PR when convenient).
Reading a message that asks something of you → append an `- ack:` line
(date + who + one-line response) rather than deleting; consumed
messages with acks older than a month may be pruned. This is for
coordination notes, warnings ("I'm about to refactor X"), and
questions — not for status (STATUS.md) or results (cards).

## Coordination with mission-control

The mission-control queue remains the *scheduler* for autonomous
drain-engine runs; this repo is the *source of truth* for program
state. Never duplicate card content across the two — pointer there,
substance here. (TASK-164 predates this system and lives whole in
mission-control; PLAN.md notes it.)

**Autopilot (since 2026-08-20): pointer cards are enqueued
automatically.** `scripts/board_feeder.py` in mission-control runs
hourly (`ophir-board-feeder.timer`, workstation): it reads this board
from origin/main, and for every card that is READY — or blocked only
on cards that are now done (stale rows count) — writes a thin pointer
task into the mission-control queue, capped at 2 open pointers and 1
workstation-pinned card in flight. The 30-minute drain cron and
macbots claim from there. Do NOT hand-enqueue pointers for board
cards; the feeder is idempotent and will not double-dispatch, but
hand-written duplicates confuse it. Kill switch:
`systemctl --user disable --now ophir-board-feeder.timer` (feeder) +
`hermes cron pause 71b7edd429c0` (drain).

## How the program is supervised (pull, not push)

Daily automated pass: the `omnibus-supervisor` Hermes cron
(workstation, 07:50) checks feeder/drain/queue health, fixes stale
board bookkeeping via small docs-only PRs, applies the stale-claim
rule, flags stuck PRs, and delivers a one-line digest to Victor's
Telegram. It never starts campaigns, never merges, never runs heavy
compute.

Fable adds judgment on top, pull-based: at session start (or when
Victor asks), Fable reads STATUS.md + `git log --oneline` since last
sync + open `agents/*` branches, and intervenes only on: contradictory
deviations, quality problems in merged work, design questions, or
anything the supervisor flagged. Anything a worker needs from Fable
goes in the inbox.
