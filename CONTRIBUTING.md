# Contributing to GluCat

Thank you for contributing to GluCat! To maintain high code quality and consistency, we enforce style formatting, static analysis, and verification checks.

## Development Workflow

Before submitting a Pull Request, please ensure that all checks pass locally.

### Local Verification Gate

You can run all local tests and linters at once using the verification script:
```bash
./verify_all.sh
```
This script will run:
1. Python checks: `ruff` and `pylint` on the `pyclical/` wrapper and demos.
2. C++ style check: `clang-format` (restricted to the `glucat/` directory).
3. C++ static analysis: `cppcheck` (restricted to the `glucat/` directory).
4. Primary regression tests: `make check-local`.

---

## Code Style & Formatting

We use `clang-format` to format C++ files strictly inside the `glucat/` directory.

### git blame Ignore Revisions

A bulk reformatting of the codebase was performed and recorded in `.git-blame-ignore-revs`. To prevent these format-only changes from cluttering your `git blame` output, please configure Git to ignore them:
```bash
git config blame.ignoreRevsFile .git-blame-ignore-revs
```

### Formatting Targets

We provide Makefile targets for style checks and formatting:
- Check for styling violations:
  ```bash
  make lint-check
  ```
- Automatically fix styling violations:
  ```bash
  make lint-fix
  ```
