#!/usr/bin/env bash
# pre-compact-save.sh — PreCompact hook: saves session context before compaction

set -euo pipefail

INPUT=$(cat)

find_learnings_dir() {
    local dir="$PWD"
    while [ "$dir" != "/" ]; do
        if [ -d "$dir/.claude/learnings" ]; then
            echo "$dir/.claude/learnings"
            return 0
        fi
        dir="$(dirname "$dir")"
    done
    return 1
}

LEARNINGS_DIR=$(find_learnings_dir) || exit 0

SESSION_NOTES="$LEARNINGS_DIR/session-notes.md"
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

SUMMARY=$(echo "$INPUT" | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    print(data.get('transcript_summary', 'Context compaction occurred — no summary available'))
except:
    print('Context compaction occurred — no summary available')
" 2>/dev/null)

cat >> "$SESSION_NOTES" << EOF

### Compaction at $TIMESTAMP
$SUMMARY

EOF

exit 0
