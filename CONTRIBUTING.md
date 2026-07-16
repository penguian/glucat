# Contributing to GluCat

Thank you for contributing to GluCat! To maintain high code quality and consistency, we enforce style formatting, static analysis, and verification checks.

## Development Workflow

Before submitting a Pull Request, please ensure that all checks pass locally.

### Local Verification Gate

You can run local tests, linters, and verification checks at once using `make check-python` or the Python verification runner:
```bash
make check-python
# or directly:
python3 verify_all.py -q --coverage --examples
```
Available flags for `verify_all.py`:
- `-q`, `--quiet`: Suppress sub-command stdout unless a step fails, emitting 1-line check status indicators.
- `--fast`: Execute ultra-fast (~0.12s) pre-commit checks: Makefile dry-run target validation (`make -n check`), license header scanning (`check_license_headers.py`), and Ruff linting (`ruff check`).
- `--coverage`: Execute C++ template header block and branch coverage tests (`make check-coverage-doctest`). Any extra CLI arguments passed to `verify_all.py` are forwarded directly to `make check-coverage-doctest`.
- `--examples`: Run full Python linter static analysis (`pylint`), PyClical test suite, Jupyter notebook validation (`pyclical/demos/validate_notebooks.py`), and interactive/plotting demo scripts (`pyclical/demos/`).

---

### Continuous Integration (GitHub Actions)

All pull requests and merges to `master` automatically trigger our GitHub Actions CI pipeline, which runs:
- **C++ Build & Test Matrix:** Builds and tests GluCat against GCC 14 and Clang 18 using both Eigen and Armadillo backends across Boost versions 1.84.0 and 1.85.0. Curated compiler warning-to-error options (e.g. `-Werror=return-type`, `-Werror=uninitialized`, `-Werror=format`) are enforced.
- **PyClical Validation Engine:** Executes `python3 verify_all.py --coverage --examples`, exercising all C++ coverage doctests, Python linters (Ruff/Pylint at 10.00/10), notebook validations, and demo script executions.
- **Documentation Build:** Installs doc dependencies (Doxygen, Graphviz, LaTeX) and verifies that documentation PDF and HTML manuals build cleanly.

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

---

## Commit & Push-level Enforcement (Git Hooks)

We provide an automated Git hooks installer script to enable repository-wide verification gates.

### Installation & Activation

To install the pre-commit and pre-push hooks in your local repository, run:
```bash
./admin/setup-git-hooks.sh
```

This installs two executable hooks in your `.git/hooks/` directory:
1. **Pre-commit hook (`.git/hooks/pre-commit`)**: Runs `python3 verify_all.py -q --fast` (~0.12s) on every local `git commit` to verify Makefile targets, license headers, and Ruff styling.
2. **Pre-push hook (`.git/hooks/pre-push`)**: Runs `python3 verify_all.py -q --coverage --examples` (~6.6s) on every `git push` to ensure zero CI regressions before pushing code to GitHub.

### Bypassing Hooks

If you need to bypass the local Git hooks for a specific action (for example, during temporary WIP commits), use the `--no-verify` flag:
```bash
git commit -m "temp commit" --no-verify
git push --no-verify
```

### Licensing & Header Policy

All source files (`.py`, `.h`, `.cpp`, `.hpp`, `.pyx`, `.pxd`) must have a valid license header.
- Core source files must include the **GNU Lesser General Public License (LGPL)** header.
- PyClical tutorials (e.g., `pyclical_tutorial_*.py` files) must include the **Creative Commons BY-SA 3.0** header. Other PyClical source files (including core wrapper modules and general demos) should include the **GNU Lesser General Public License (LGPL)** header.

> [!NOTE]
> **Unique Exemption:** The demo file [pyclical/demos/plotting/plotting_demo_dialog.py](pyclical/demos/plotting/plotting_demo_dialog.py) is a unique exemption. Since it is a derivative work of a TraitsUI/Mayavi template copyrighted by Enthought, Inc., it maintains its original **BSD Style** license header. This is the only BSD-licensed file permitted in the repository.
