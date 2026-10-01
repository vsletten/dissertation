---
name: pr-merge-gauntlet
description: Land a PR in this repo — auto-merge, required checks, up-to-date branch, and mandatory Sourcery thread resolution. Use whenever opening a PR here or when a PR sits OPEN/BLOCKED with green checks.
---

# Landing a PR in vsletten/dissertation

Branch protection on `main`: 4 required checks (Analyze c-cpp, Analyze
javascript-typescript, CodeQL, Sourcery review), **strict up-to-date**,
and **required conversation resolution**. The last two cause silent
stalls; this skill exists because they burned us (PR #12, #14).

## Standard flow

```bash
gh pr create --title "..." --body "..."   # body: ## Summary / ## Test plan / ## Breaking changes
gh pr merge --auto --squash
```

Then poll until MERGED — do not assume:

```bash
gh pr view <N> --json state,mergeStateStatus --jq '.state+" "+.mergeStateStatus'
```

## Diagnosing a stuck PR

| Symptom | Cause | Fix |
|---|---|---|
| `BEHIND` | another PR landed first (concurrent workers are normal here) | `git fetch origin main && git merge origin/main && git push` (this box's gh lacks `pr update-branch`) |
| `BLOCKED`, all checks green | **unresolved Sourcery threads** | see below |
| `gh pr merge` says "base branch policy prohibits" | same | same |
| checks "expected", not starting | CodeQL lag | wait; never hand-retry |

## Sourcery threads (the one everyone forgets)

List unresolved:

```bash
gh api graphql -f query='{repository(owner:"vsletten",name:"dissertation"){pullRequest(number:<N>){reviewThreads(first:40){nodes{id isResolved path line comments(first:1){nodes{body}}}}}}}' --jq '[.data.repository.pullRequest.reviewThreads.nodes[] | select(.isResolved|not)]'
```

Policy: **fix legitimate points** (they usually are — validate inputs,
edge cases, tests), push, then resolve each thread:

```bash
gh api graphql -f query='mutation{resolveReviewThread(input:{threadId:"<ID>"}){thread{isResolved}}}'
```

Never resolve-without-reading to force a merge. GitHub 500s
intermittently — retry with backoff.

## House rules riding along

Commits end with the Co-Authored-By line your harness specifies. If the
PR implements a program card, it must also carry the card update +
`docs/program/STATUS.md` line (see docs/program/PROTOCOL.md).
