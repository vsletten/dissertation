---
name: project-info
description: Show comprehensive project context including type, commands, dependencies, git state, and services
user-invocable: true
disable-model-invocation: true
allowed-tools: Bash, Read, Glob, Grep
---

# Project Info

Run the project discovery script to analyze the current working directory and display comprehensive project context.

## Steps

1. Run the discovery script:
```bash
~/.claude/scripts/project-discover.sh "$(pwd)"
```

2. Present the results in a clean, readable format with sections for:
   - Project name and type
   - Framework and package manager
   - Available commands (build, test, lint, dev server)
   - Key dependencies
   - Docker services (if any)
   - Git branch and status
   - Whether a CLAUDE.md exists

3. If a CLAUDE.md exists, read it and provide a brief summary of key project guidance.

4. If there are additional subdirectories with their own package.json or pyproject.toml, mention them as sub-projects.
