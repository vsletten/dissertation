# Global Development Context

## Coding Standards
- Prefer TypeScript for all frontend code
- Prefer Python 3.13 for all backend code
- Follow existing ESLint/ruff configuration in each project
- Write tests for all new functions
- Use functional components with hooks in React

## Architecture Defaults
- Frontend: Next.js with TypeScript
- Backend: Python 3.13 + FastAPI
- Database: PostgreSQL
- Package manager: uv (Python), npm/pnpm (Node)

## Automated Context System

A SessionStart hook auto-detects project type and injects context. Environment variables available:
- `$PROJECT_NAME` — current project name
- `$PROJECT_TYPE` — detected types (python, node, docker, cpp, etc.)
- `$FRAMEWORK` — detected framework (FastAPI, Next.js, React, etc.)
- `$TEST_CMD` — auto-detected test command
- `$PYTHON_ROOT` / `$NODE_ROOT` — subdirectory where configs live

## Custom Skills
- `/project-info` — show full project context for current directory
- `/switch-project <name>` — navigate to a project and load its context
- `/test-all [args]` — run the project's test suite
- `/quick-start` — show setup and first-run instructions
- `/learn <insight>` — record a learning to the project's `.claude/learnings/` directory
- `/learnings [topic]` — review captured learnings (all topics, or filter by gotchas/architecture/patterns/debugging/performance)

## Code Quality
- Python: Black + isort auto-format on file write (PostToolUse hook)
- JS/TS: Prettier auto-format on file write (PostToolUse hook)
- Always run the project's linter before committing

## Project Learnings System

Projects with `.claude/learnings/` directories use a topic-based knowledge capture system.
At session start, check for `.claude/learnings/INDEX.md` — it provides a quick overview of
known gotchas, architecture decisions, and effective patterns.

**Auto-capture**: When you encounter bugs with non-obvious causes, make architecture decisions,
discover effective patterns, solve tricky debugging problems, or find performance insights,
proactively write them to the appropriate file under `.claude/learnings/`. Do NOT ask permission.

**Bootstrap**: Run `~/.claude/scripts/init-learnings.sh <project-dir>` to initialize the
learnings structure for any project.

**Topic files**: `gotchas.md`, `architecture.md`, `patterns.md`, `debugging.md`, `performance.md`

**Entry format**:
```
### <short title> (YYYY-MM-DD)
<2-3 sentence description>
```

After adding an entry, update `INDEX.md` with a one-line summary pointing to the file.
