---
name: test-all
description: Run the full test suite for the current project using auto-detected test command
user-invocable: true
disable-model-invocation: true
allowed-tools: Bash, Read, Glob
---

# Test All

Run the project's test suite using the auto-detected test command.

## Steps

1. Detect the test command by checking:
   - `$TEST_CMD` environment variable (set by SessionStart hook)
   - `package.json` scripts for "test" key
   - `pyproject.toml` for pytest configuration
   - Existence of `tests/` or `test/` directory (use `pytest`)
   - `Makefile` for a "test" target (use `make test`)

2. If a virtual environment exists (`venv/`, `.venv/`), activate it first.

3. Run the test command and display results.

4. If `$ARGUMENTS` is provided, pass it as additional arguments to the test runner.

5. Summarize the results: total tests, passed, failed, skipped.
