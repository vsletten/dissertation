---
name: quick-start
description: Show how to get started with the current project, including setup and first run
user-invocable: true
disable-model-invocation: true
allowed-tools: Bash, Read, Glob, Grep
---

# Quick Start

Generate a quick-start guide for the current project by analyzing its configuration.

## Steps

1. Run the discovery script:
   ```bash
   ~/.claude/scripts/project-discover.sh "$(pwd)"
   ```

2. Read the project's CLAUDE.md and README.md if they exist.

3. Based on the discovered project type, generate a step-by-step quick-start guide covering:
   - Prerequisites
   - Dependency installation
   - Environment setup
   - Running the application or tests

4. Show the commands as a numbered list the user can follow step by step.
