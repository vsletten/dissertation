---
name: learn
description: Record a learning, decision, or note to the project's learnings directory
user-invocable: true
disable-model-invocation: true
---

# /learn — Record a project learning

The user wants to record a specific learning for the current project. The argument contains what they want to remember.

## Steps

1. **Find the learnings directory**: Look for `.claude/learnings/` in the current project (walk up from CWD if needed). If it doesn't exist, create it by running `~/.claude/scripts/init-learnings.sh` on the project root.

2. **Categorize the learning**: Based on the content, determine the best topic file:
   - `gotchas.md` — bugs, pitfalls, things that break in non-obvious ways
   - `architecture.md` — design decisions, rationale for choosing approach A over B
   - `patterns.md` — effective patterns, reusable approaches that work well
   - `debugging.md` — debugging approaches, "when X happens check Y"
   - `performance.md` — performance findings, benchmarks, optimizations

   If unclear, default to `gotchas.md`.

3. **Write the entry** to the appropriate topic file using this format:
   ```
   ### <short descriptive title> (YYYY-MM-DD)
   <The learning, written as 2-3 clear sentences>
   ```

4. **Update INDEX.md**: Add a one-line summary under the appropriate section with a pointer to the topic file.

5. **Confirm** to the user what was recorded and where.
