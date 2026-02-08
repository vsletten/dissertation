---
name: switch-project
description: Navigate to a project directory and load its context
user-invocable: true
disable-model-invocation: true
allowed-tools: Bash, Read, Glob, Grep
---

# Switch Project

Navigate to a project directory and load its context.

## Arguments

$ARGUMENTS — The project name or path to switch to.

## Steps

1. Find the project directory. Search common locations:
   - Current workspace and subdirectories
   - $HOME/projects/ and subdirectories
   - Anywhere under $HOME that matches the name

2. If no match is found, list available project directories and ask the user which one they meant.

3. Once found, `cd` to the project directory.

4. Run the discovery script:
   ```bash
   ~/.claude/scripts/project-discover.sh "\$(pwd)"
   ```

5. If a CLAUDE.md exists, read it and present a concise summary.

6. Report the full path so the user knows where they are.
