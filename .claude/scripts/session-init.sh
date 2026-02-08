#!/bin/bash
# session-init.sh — SessionStart hook that runs project discovery and injects context.
# Called by Claude Code on session startup, resume, and compact events.

set +e  # Never fail — must not block Claude startup

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

input=$(cat)

cwd=$(echo "$input" | jq -r '.cwd // empty' 2>/dev/null)
cwd="${cwd:-$(pwd)}"

discovery=$("$SCRIPT_DIR/project-discover.sh" "$cwd" 2>/dev/null)

if [[ -z "$discovery" ]]; then
    exit 0
fi

project_name=$(echo "$discovery" | grep '^PROJECT_NAME=' | cut -d= -f2-)
project_type=$(echo "$discovery" | grep '^PROJECT_TYPE=' | cut -d= -f2-)
python_root=$(echo "$discovery" | grep '^PYTHON_ROOT=' | cut -d= -f2-)
node_root=$(echo "$discovery" | grep '^NODE_ROOT=' | cut -d= -f2-)
test_cmd=$(echo "$discovery" | grep '^TEST_CMD=' | cut -d= -f2-)
framework=$(echo "$discovery" | grep '^FRAMEWORK=' | cut -d= -f2-)
branch=$(echo "$discovery" | grep '^BRANCH=' | cut -d= -f2-)
dirty=$(echo "$discovery" | grep '^DIRTY_FILES=' | cut -d= -f2-)
has_claude_md=$(echo "$discovery" | grep '^HAS_CLAUDE_MD=' | cut -d= -f2-)

if [[ -n "$CLAUDE_ENV_FILE" ]]; then
    {
        echo "export PROJECT_NAME='$project_name'"
        echo "export PROJECT_TYPE='$project_type'"
        [[ -n "$python_root" ]] && echo "export PYTHON_ROOT='$python_root'"
        [[ -n "$node_root" ]] && echo "export NODE_ROOT='$node_root'"
        [[ -n "$test_cmd" ]] && echo "export TEST_CMD='$test_cmd'"
        [[ -n "$framework" ]] && echo "export FRAMEWORK='$framework'"
        [[ -n "$branch" ]] && echo "export GIT_BRANCH='$branch'"
    } >> "$CLAUDE_ENV_FILE"
fi

{
    echo ""
    echo "--- Project Context Loaded ---"
    echo "Project: $project_name"
    echo "Type: $project_type"
    [[ -n "$framework" ]] && echo "Framework: $framework"
    [[ -n "$python_root" ]] && echo "Python root: $python_root/"
    [[ -n "$node_root" ]] && echo "Node root: $node_root/"
    [[ -n "$branch" ]] && echo "Branch: $branch (${dirty:-0} dirty files)"
    [[ -n "$test_cmd" ]] && echo "Test: $test_cmd"
    [[ -n "$has_claude_md" ]] && echo "CLAUDE.md: present"
    echo "---"
} >&2

exit 0
