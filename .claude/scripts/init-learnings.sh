#!/usr/bin/env bash
# init-learnings.sh — Bootstrap the .claude/learnings/ structure for a project
# Usage: init-learnings.sh <project-dir>

set -euo pipefail

PROJECT_DIR="${1:?Usage: init-learnings.sh <project-dir>}"
PROJECT_DIR="$(cd "$PROJECT_DIR" && pwd)"

if [ ! -d "$PROJECT_DIR" ]; then
    echo "Error: $PROJECT_DIR does not exist"
    exit 1
fi

LEARNINGS_DIR="$PROJECT_DIR/.claude/learnings"
RULES_DIR="$PROJECT_DIR/.claude/rules"
TODAY=$(date +%Y-%m-%d)

echo "Initializing learnings structure in $PROJECT_DIR"

mkdir -p "$LEARNINGS_DIR"
mkdir -p "$RULES_DIR"

if [ ! -f "$LEARNINGS_DIR/INDEX.md" ]; then
    cat > "$LEARNINGS_DIR/INDEX.md" << EOF
# Project Learnings Index
Last updated: $TODAY

## Active Gotchas
_None yet — gotchas will appear here as they're discovered._

## Key Architecture Decisions
_None yet — decisions will appear here as they're made._

## Effective Patterns
_None yet — patterns will appear here as they're discovered._
EOF
    echo "  Created INDEX.md"
fi

for topic in gotchas architecture patterns debugging performance; do
    FILE="$LEARNINGS_DIR/$topic.md"
    if [ ! -f "$FILE" ]; then
        TITLE=$(echo "$topic" | sed 's/.*/\u&/')
        cat > "$FILE" << EOF
# $TITLE

Entries are added automatically during development sessions.
Format: ### <title> (YYYY-MM-DD) followed by 2-3 sentence description.

---
EOF
        echo "  Created $topic.md"
    fi
done

if [ ! -f "$LEARNINGS_DIR/session-notes.md" ]; then
    cat > "$LEARNINGS_DIR/session-notes.md" << EOF
# Session Notes

Auto-captured by PreCompact hook before context compaction.
These are raw notes — review and promote useful ones to topic files.

---
EOF
    echo "  Created session-notes.md"
fi

# Create path-specific rules if rules dir is empty
if [ -z "$(ls -A "$RULES_DIR" 2>/dev/null)" ]; then
    if [ -f "$PROJECT_DIR/pyproject.toml" ] || [ -f "$PROJECT_DIR/setup.py" ] || \
       find "$PROJECT_DIR" -maxdepth 3 -name "*.py" -print -quit 2>/dev/null | grep -q .; then
        cat > "$RULES_DIR/database.md" << 'RULES_EOF'
---
paths:
  - "**/*model*.py"
  - "**/*migration*.py"
  - "**/alembic/**"
  - "**/db/**"
---
# Database Rules
- Check @.claude/learnings/gotchas.md for known DB pitfalls
- Check @.claude/learnings/patterns.md for DB patterns
RULES_EOF

        cat > "$RULES_DIR/api.md" << 'RULES_EOF'
---
paths:
  - "**/api/**"
  - "**/routes/**"
  - "**/endpoints/**"
---
# API Development
- Check @.claude/learnings/patterns.md for API patterns
- Check @.claude/learnings/gotchas.md for known API pitfalls
RULES_EOF
        echo "  Created default rules (database.md, api.md)"
    fi

    if [ -f "$PROJECT_DIR/package.json" ] || \
       find "$PROJECT_DIR" -maxdepth 3 -name "*.tsx" -print -quit 2>/dev/null | grep -q .; then
        cat > "$RULES_DIR/frontend.md" << 'RULES_EOF'
---
paths:
  - "**/*.tsx"
  - "**/*.jsx"
  - "**/components/**"
  - "**/pages/**"
  - "**/app/**"
---
# Frontend Development
- Check @.claude/learnings/patterns.md for UI patterns
- Check @.claude/learnings/gotchas.md for known frontend pitfalls
RULES_EOF
        echo "  Created frontend rule"
    fi
fi

# Update .gitignore
GITIGNORE="$PROJECT_DIR/.gitignore"
if [ -f "$GITIGNORE" ]; then
    if ! grep -q "session-notes.md" "$GITIGNORE" 2>/dev/null; then
        echo "" >> "$GITIGNORE"
        echo "# Claude learnings - session notes (auto-generated, not version controlled)" >> "$GITIGNORE"
        echo ".claude/learnings/session-notes.md" >> "$GITIGNORE"
        echo "  Updated .gitignore"
    fi
else
    cat > "$GITIGNORE" << 'GI_EOF'
# Claude learnings - session notes (auto-generated, not version controlled)
.claude/learnings/session-notes.md
GI_EOF
    echo "  Created .gitignore"
fi

# Update CLAUDE.md
CLAUDE_MD="$PROJECT_DIR/CLAUDE.md"
if [ -f "$CLAUDE_MD" ]; then
    if ! grep -q "@.claude/learnings/INDEX.md" "$CLAUDE_MD" 2>/dev/null; then
        cat >> "$CLAUDE_MD" << 'SECTION'

@.claude/learnings/INDEX.md

## Recording Learnings

When you encounter any of the following during work, **proactively record it** in the
appropriate file under `.claude/learnings/`:

| What happened | Write to | Example |
|---|---|---|
| Bug caused by a non-obvious reason | `gotchas.md` | "Import X breaks because of Y" |
| Chose approach A over B with rationale | `architecture.md` | "Used connection pooling because..." |
| Found a pattern that works well | `patterns.md` | "Batch processing with chunk_size=500..." |
| Solved a tricky debugging problem | `debugging.md` | "When error X appears, check Y first" |
| Found performance insight | `performance.md` | "Query Z takes 30s without index on col" |

**Format for entries:**
```
### <short title> (YYYY-MM-DD)
<2-3 sentence description of the learning>
```

After adding an entry, update `.claude/learnings/INDEX.md` with a one-line summary pointing to the file.
Do NOT ask for permission — just write the learning when you encounter it.
SECTION
        echo "  Updated CLAUDE.md with @import and auto-capture section"
    else
        echo "  CLAUDE.md already has learnings import"
    fi
else
    echo "  Warning: No CLAUDE.md found — skipping CLAUDE.md update"
fi

echo ""
echo "Done! Learnings structure initialized."
