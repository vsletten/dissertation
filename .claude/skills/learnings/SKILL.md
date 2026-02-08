---
name: learnings
description: Show all recorded learnings for the current project
user-invocable: true
disable-model-invocation: true
---

# /learnings — Review project learnings

The user wants to see what learnings have been captured for the current project.

## Steps

1. **Find the learnings directory**: Look for `.claude/learnings/` in the current project. If it doesn't exist, tell the user and suggest running `/learn` to initialize it.

2. **Check for a topic filter**: If the user provided an argument (e.g., `/learnings gotchas`), only show that specific topic file. Valid topics: `gotchas`, `architecture`, `patterns`, `debugging`, `performance`.

3. **Display the learnings**:

   **If no filter** — show INDEX.md as overview, then each non-empty topic file.

   **If filtered** — show only that specific topic file.

4. **Format**: Display the learnings in a clean, readable format with topic name as header.
