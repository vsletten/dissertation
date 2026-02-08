#!/bin/bash
# Claude Code status line - shows cwd, git branch/status, context remaining, model

input=$(cat)

cwd=$(echo "$input" | jq -r '.cwd // empty')
cwd="${cwd/#$HOME/\~}"

branch=""
git_status=""
if git -C "${cwd/#\~/$HOME}" --no-optional-locks rev-parse --git-dir > /dev/null 2>&1; then
    branch=$(git -C "${cwd/#\~/$HOME}" --no-optional-locks rev-parse --abbrev-ref HEAD 2>/dev/null)

    if [ -n "$branch" ]; then
        if [ -z "$(git -C "${cwd/#\~/$HOME}" --no-optional-locks status --porcelain 2>/dev/null)" ]; then
            git_status="\033[32m✓"
        else
            git_status="\033[33m*"
        fi
    fi
fi

remaining=$(echo "$input" | jq -r '.context_window.remaining_percentage // empty')
model=$(echo "$input" | jq -r '.model.display_name // empty')

project_label=""
if [ -n "$PROJECT_NAME" ] && [ -n "$PROJECT_TYPE" ]; then
    project_label="$PROJECT_NAME"
    [ -n "$FRAMEWORK" ] && project_label="$project_label/$FRAMEWORK"
fi

if [ -n "$project_label" ]; then
    printf '\033[01;35m%s\033[00m' "$project_label"
    printf ' \033[2m%s\033[00m' "$cwd"
else
    printf '\033[01;32m%s\033[00m' "$cwd"
fi

if [ -n "$branch" ]; then
    printf ' \033[36m[\033[00m%s %b\033[36m]\033[00m' "$branch" "$git_status"
fi

if [ -n "$remaining" ]; then
    remaining_int=${remaining%.*}
    if [ "$remaining_int" -gt 50 ] 2>/dev/null; then
        color='\033[32m'
    elif [ "$remaining_int" -gt 20 ] 2>/dev/null; then
        color='\033[33m'
    else
        color='\033[31m'
    fi
    printf ' %bCtx: %s%%\033[00m' "$color" "$remaining"
fi

if [ -n "$model" ]; then
    printf ' \033[2m%s\033[00m' "$model"
fi
