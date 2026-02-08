#!/bin/bash
# project-discover.sh — Scans the current directory and outputs structured project context.
# Used by session-init.sh and the /project-info skill.
#
# Output format: structured text to stdout, one section per detected feature.
# Exit code: 0 always (never block Claude startup).

set +e  # Never fail — this must not block Claude startup

DIR="${1:-$(pwd)}"

# ─── Helpers ───────────────────────────────────────────────────────────────────

has_file()  { [[ -f "$DIR/$1" ]]; }
has_dir()   { [[ -d "$DIR/$1" ]]; }
first_n()   { head -n "${1:-10}"; }

# ─── Project name ─────────────────────────────────────────────────────────────

project_name=$(basename "$DIR")
echo "PROJECT_NAME=$project_name"

# ─── Project type detection ───────────────────────────────────────────────────
# Check root AND common subdirectories (backend/, frontend/, src/, v1/, app/)

types=()
python_root=""
node_root=""

# Search paths: root first, then common subdirs
search_paths=("" "backend/" "frontend/" "src/" "v1/" "app/" "server/" "client/" "api/" "web/" "model/" "viewer/" "lib/" "cmd/" "pkg/")

for sp in "${search_paths[@]}"; do
    if has_file "${sp}pyproject.toml" || has_file "${sp}setup.py" || has_file "${sp}setup.cfg" || has_file "${sp}requirements.txt" || has_file "${sp}Pipfile"; then
        if [[ ! " ${types[*]} " == *" python "* ]]; then
            types+=("python")
        fi
        if [[ -z "$python_root" ]]; then
            python_root="${sp%/}"
        fi
    fi

    if has_file "${sp}package.json"; then
        if [[ ! " ${types[*]} " == *" node "* ]]; then
            types+=("node")
        fi
        if [[ -z "$node_root" ]]; then
            node_root="${sp%/}"
        fi
    fi
done

if has_file "tsconfig.json"; then
    types+=("typescript")
fi

if has_file "Cargo.toml"; then
    types+=("rust")
fi

if has_file "go.mod"; then
    types+=("go")
fi

cpp_found=false
for sp in "${search_paths[@]}"; do
    if has_file "${sp}Makefile" || has_file "${sp}CMakeLists.txt"; then
        cpp_found=true
        break
    fi
done
if $cpp_found; then
    types+=("cpp")
fi

if has_file "docker-compose.yml" || has_file "docker-compose.yaml" || has_file "compose.yml" || has_file "compose.yaml"; then
    types+=("docker")
fi

if has_file "Dockerfile"; then
    types+=("dockerfile")
fi

if [[ ${#types[@]} -eq 0 ]]; then
    types+=("unknown")
fi

echo "PROJECT_TYPE=${types[*]}"
[[ -n "$python_root" ]] && echo "PYTHON_ROOT=$python_root"
[[ -n "$node_root" ]] && echo "NODE_ROOT=$node_root"

# ─── Python project details ──────────────────────────────────────────────────

if [[ " ${types[*]} " == *" python "* ]]; then
    echo ""
    echo "=== Python Project ==="

    py_dir="$DIR"
    [[ -n "$python_root" ]] && py_dir="$DIR/$python_root"

    if [[ -f "$py_dir/pyproject.toml" ]]; then
        echo "CONFIG=${python_root:+$python_root/}pyproject.toml"
        if grep -q '\[tool\.uv\]' "$py_dir/pyproject.toml" 2>/dev/null; then
            echo "PKG_MANAGER=uv"
        elif grep -q '\[tool\.poetry\]' "$py_dir/pyproject.toml" 2>/dev/null; then
            echo "PKG_MANAGER=poetry"
        else
            echo "PKG_MANAGER=pip"
        fi
        if grep -q '\[project\.scripts\]' "$py_dir/pyproject.toml" 2>/dev/null; then
            echo "SCRIPTS:"
            sed -n '/\[project\.scripts\]/,/^\[/p' "$py_dir/pyproject.toml" | grep -v '^\[' | grep '=' | first_n 10
        fi
    elif [[ -f "$py_dir/setup.py" ]]; then
        echo "CONFIG=${python_root:+$python_root/}setup.py"
        echo "PKG_MANAGER=pip"
    elif [[ -f "$py_dir/requirements.txt" ]]; then
        echo "CONFIG=${python_root:+$python_root/}requirements.txt"
        echo "PKG_MANAGER=pip"
    fi

    for check_dir in "$DIR" "$py_dir"; do
        for venv_name in venv .venv env; do
            if [[ -d "$check_dir/$venv_name" ]]; then
                rel="${check_dir#$DIR}"
                rel="${rel#/}"
                if [[ -n "$rel" ]]; then
                    echo "VENV=$rel/$venv_name"
                else
                    echo "VENV=$venv_name"
                fi
                break 2
            fi
        done
    done

    if [[ -f "$py_dir/pyproject.toml" ]] && grep -q '\[tool\.pytest' "$py_dir/pyproject.toml" 2>/dev/null; then
        echo "TEST_CMD=cd ${python_root:-.} && pytest"
    elif [[ -d "$py_dir/tests" ]] || [[ -d "$py_dir/test" ]]; then
        echo "TEST_CMD=cd ${python_root:-.} && pytest"
    fi

    lint_tools=()
    if [[ -f "$py_dir/pyproject.toml" ]]; then
        grep -q '\[tool\.ruff\]' "$py_dir/pyproject.toml" 2>/dev/null && lint_tools+=("ruff")
        grep -q '\[tool\.black\]' "$py_dir/pyproject.toml" 2>/dev/null && lint_tools+=("black")
        grep -q '\[tool\.isort\]' "$py_dir/pyproject.toml" 2>/dev/null && lint_tools+=("isort")
        grep -q '\[tool\.mypy\]' "$py_dir/pyproject.toml" 2>/dev/null && lint_tools+=("mypy")
    fi
    [[ -f "$py_dir/.flake8" ]] && lint_tools+=("flake8")
    [[ -f "$py_dir/mypy.ini" ]] && lint_tools+=("mypy")
    [[ -f "$py_dir/.ruff.toml" ]] && lint_tools+=("ruff")
    if [[ ${#lint_tools[@]} -gt 0 ]]; then
        echo "LINT_TOOLS=$(echo "${lint_tools[*]}" | tr ' ' '\n' | sort -u | tr '\n' ' ')"
    fi

    if [[ -f "$py_dir/requirements.txt" ]]; then
        echo "KEY_DEPS:"
        grep -v '^#' "$py_dir/requirements.txt" | grep -v '^$' | sed 's/[>=<].*//' | first_n 15
    elif [[ -f "$py_dir/pyproject.toml" ]]; then
        echo "KEY_DEPS:"
        sed -n '/^dependencies/,/^\]/p' "$py_dir/pyproject.toml" 2>/dev/null | grep -v '^\[' | grep -v '^dependencies' | sed "s/[>=<\"',].*//" | tr -d ' ' | grep -v '^$' | first_n 15
    fi

    for dep_file in "$py_dir/requirements.txt" "$py_dir/pyproject.toml"; do
        if [[ -f "$dep_file" ]]; then
            grep -qi 'fastapi' "$dep_file" 2>/dev/null && echo "FRAMEWORK=FastAPI"
            grep -qi 'django' "$dep_file" 2>/dev/null && echo "FRAMEWORK=Django"
            grep -qi 'flask' "$dep_file" 2>/dev/null && echo "FRAMEWORK=Flask"
            break
        fi
    done
fi

# ─── Node/TypeScript project details ─────────────────────────────────────────

if [[ " ${types[*]} " == *" node "* ]]; then
    echo ""
    echo "=== Node/TypeScript Project ==="

    nd_dir="$DIR"
    [[ -n "$node_root" ]] && nd_dir="$DIR/$node_root"

    if [[ -f "$nd_dir/pnpm-lock.yaml" ]]; then
        echo "PKG_MANAGER=pnpm"
    elif [[ -f "$nd_dir/yarn.lock" ]]; then
        echo "PKG_MANAGER=yarn"
    elif [[ -f "$nd_dir/bun.lockb" ]]; then
        echo "PKG_MANAGER=bun"
    else
        echo "PKG_MANAGER=npm"
    fi

    if [[ -f "$nd_dir/package.json" ]]; then
        scripts=$(python3 -c "
import json, sys
try:
    with open('$nd_dir/package.json') as f:
        pkg = json.load(f)
    scripts = pkg.get('scripts', {})
    if scripts:
        for k, v in list(scripts.items())[:15]:
            print(f'  {k}: {v}')
except: pass
" 2>/dev/null)
        if [[ -n "$scripts" ]]; then
            echo "SCRIPTS:"
            echo "$scripts"
        fi

        deps=$(python3 -c "
import json
try:
    with open('$nd_dir/package.json') as f:
        pkg = json.load(f)
    deps = list(pkg.get('dependencies', {}).keys())[:15]
    if deps:
        for d in deps:
            print(f'  {d}')
except: pass
" 2>/dev/null)
        if [[ -n "$deps" ]]; then
            echo "KEY_DEPS:"
            echo "$deps"
        fi

        python3 -c "
import json
try:
    with open('$nd_dir/package.json') as f:
        pkg = json.load(f)
    all_deps = {**pkg.get('dependencies', {}), **pkg.get('devDependencies', {})}
    if 'next' in all_deps: print('FRAMEWORK=Next.js')
    elif 'react' in all_deps: print('FRAMEWORK=React')
    elif 'vue' in all_deps: print('FRAMEWORK=Vue')
    elif 'svelte' in all_deps: print('FRAMEWORK=Svelte')
    elif '@angular/core' in all_deps: print('FRAMEWORK=Angular')
    elif 'express' in all_deps: print('FRAMEWORK=Express')
except: pass
" 2>/dev/null
    fi

    if [[ -f "$nd_dir/tsconfig.json" ]] || has_file "tsconfig.json"; then
        echo "TYPESCRIPT=true"
    fi
fi

# ─── C++/Make project details ────────────────────────────────────────────────

if [[ " ${types[*]} " == *" cpp "* ]]; then
    echo ""
    echo "=== C++/Make Project ==="

    if has_file "CMakeLists.txt"; then
        echo "BUILD_SYSTEM=cmake"
        echo "BUILD_CMD=cmake --build build"
    elif has_file "Makefile"; then
        echo "BUILD_SYSTEM=make"
        echo "BUILD_CMD=make"
        targets=$(grep -E '^[a-zA-Z_-]+:' "$DIR/Makefile" 2>/dev/null | sed 's/:.*//' | first_n 10)
        if [[ -n "$targets" ]]; then
            echo "MAKE_TARGETS:"
            echo "$targets" | sed 's/^/  /'
        fi
    fi
fi

# ─── Docker details ──────────────────────────────────────────────────────────

if [[ " ${types[*]} " == *" docker "* ]]; then
    echo ""
    echo "=== Docker Services ==="

    compose_file=""
    for f in docker-compose.yml docker-compose.yaml compose.yml compose.yaml; do
        if has_file "$f"; then
            compose_file="$f"
            break
        fi
    done

    if [[ -n "$compose_file" ]]; then
        services=$(grep -E '^\s{2}[a-zA-Z_-]+:' "$DIR/$compose_file" 2>/dev/null | sed 's/:.*//' | tr -d ' ' | first_n 15)
        if [[ -n "$services" ]]; then
            echo "SERVICES:"
            echo "$services" | sed 's/^/  /'
        fi

        if command -v docker &>/dev/null; then
            running=$(docker compose -f "$DIR/$compose_file" ps --format '{{.Name}} {{.Status}}' 2>/dev/null | first_n 10)
            if [[ -n "$running" ]]; then
                echo "RUNNING_CONTAINERS:"
                echo "$running" | sed 's/^/  /'
            fi
        fi
    fi
fi

# ─── Git state ────────────────────────────────────────────────────────────────

if has_dir ".git" || git -C "$DIR" rev-parse --git-dir &>/dev/null; then
    echo ""
    echo "=== Git ==="

    branch=$(git -C "$DIR" --no-optional-locks rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
    echo "BRANCH=$branch"

    dirty=$(git -C "$DIR" --no-optional-locks status --porcelain 2>/dev/null | wc -l | tr -d ' ')
    echo "DIRTY_FILES=$dirty"

    echo "RECENT_COMMITS:"
    git -C "$DIR" --no-optional-locks log --oneline -5 2>/dev/null | sed 's/^/  /' || true

    remote=$(git -C "$DIR" remote get-url origin 2>/dev/null || echo "none")
    echo "REMOTE=$remote"
fi

# ─── Existing CLAUDE.md ──────────────────────────────────────────────────────

if has_file "CLAUDE.md"; then
    echo ""
    echo "=== CLAUDE.md ==="
    echo "HAS_CLAUDE_MD=true"
    lines=$(wc -l < "$DIR/CLAUDE.md")
    echo "CLAUDE_MD_LINES=$lines"
    echo "CLAUDE_MD_SUMMARY:"
    head -20 "$DIR/CLAUDE.md" | sed 's/^/  /'
elif has_file ".claude/CLAUDE.md"; then
    echo ""
    echo "=== CLAUDE.md ==="
    echo "HAS_CLAUDE_MD=true (in .claude/)"
    lines=$(wc -l < "$DIR/.claude/CLAUDE.md")
    echo "CLAUDE_MD_LINES=$lines"
fi

# ─── Directory structure (top level only) ─────────────────────────────────────

echo ""
echo "=== Directory Structure ==="
echo "TOP_LEVEL:"
ls -1 "$DIR" | grep -v node_modules | grep -v __pycache__ | grep -v '.git$' | grep -v venv | grep -v .venv | grep -v dist | grep -v build | grep -v '.egg-info' | first_n 25 | sed 's/^/  /'
