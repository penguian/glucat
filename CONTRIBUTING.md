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
1. Python checks: `ruff` and `pylint` on the `pyclical/` wrapper, demos, and benchmarks.
2. C++ style check: `clang-format` (restricted to the `glucat/` directory).
3. C++ static analysis: `cppcheck` (restricted to the `glucat/` directory).
4. Primary regression tests: `make check-local`.
5. Jupyter Notebook validation: Rebuilding and checking generated tutorials and demos against JSON syntax and `nbformat` schema via:
   ```bash
   python3 pyclical/demos/validate_notebooks.py
   ```

### Continuous Integration (GitHub Actions)

All pull requests and merges to `master` automatically trigger our GitHub Actions CI pipeline, which runs:
- **C++ Build & Test Matrix:** Builds and tests GluCat against GCC 14 and Clang 18 using both Eigen and Armadillo backends across Boost versions 1.84.0 and 1.85.0. Curated compiler warning-to-error options (e.g. `-Werror=return-type`, `-Werror=uninitialized`, `-Werror=format`) are enforced.
- **PyClical Validation:** Lints Python files with Ruff/Pylint, validates generated Jupyter notebooks, and runs Python doctests.
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

## Commit-level Enforcement (Pre-commit Hooks)

We use `pre-commit` to automate code quality and compliance checks at the commit level.

### Installation & Activation

To install the hooks in your local Git repository:
1. Install `pre-commit` via pip:
   ```bash
   pip install pre-commit
   ```
2. Activate the pre-commit hooks:
   ```bash
   pre-commit install
   ```

Once installed, these hooks run automatically every time you execute `git commit`. If a hook fails (e.g., due to trailing whitespace, syntax issues, or missing license headers), the commit will be rejected, and any autofixes will be applied to your workspace files.

### Bypassing Hooks

If you need to bypass the pre-commit hooks for a specific commit (for example, during temporary WIP commits), use the `--no-verify` flag:
```bash
git commit -m "temp commit" --no-verify
```

### Licensing & Header Policy

All source files (`.py`, `.h`, `.cpp`, `.hpp`, `.pyx`, `.pxd`) must have a valid license header.
- Core source files must include the **GNU Lesser General Public License (LGPL)** header.
- PyClical tutorials (e.g., `pyclical_tutorial_*.py` files) must include the **Creative Commons BY-SA 3.0** header. Other PyClical source files (including core wrapper modules and general demos) should include the **GNU Lesser General Public License (LGPL)** header.

> [!NOTE]
> **Unique Exemption:** The demo file [pyclical/demos/plotting_demo_dialog.py](pyclical/demos/plotting_demo_dialog.py) is a unique exemption. Since it is a derivative work of a TraitsUI/Mayavi template copyrighted by Enthought, Inc., it maintains its original **BSD Style** license header. This is the only BSD-licensed file permitted in the repository.
