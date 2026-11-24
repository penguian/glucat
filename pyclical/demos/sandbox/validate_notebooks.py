#!/usr/bin/env python3
"""
Build and validate Jupyter and Marimo notebook files for PyClical demos.
"""
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# validate_notebooks.py: Build and validate Jupyter and Marimo notebooks.
#
#    copyright            : (C) 2026 by Paul C. Leopardi
#
#    This library is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library.  If not, see <http://www.gnu.org/licenses/>.
#

import glob
import importlib.util
import json
import os
import subprocess
import sys
import nbformat
from nbformat.validator import validate, NotebookValidationError


def _build_python_env(script_dir):
    """Build environment dict with sandbox, demos, and pyclical on PYTHONPATH."""
    env = os.environ.copy()
    demos_dir = os.path.dirname(script_dir)
    pyclical_dir = os.path.dirname(demos_dir)
    extra_paths = [script_dir, demos_dir, pyclical_dir]
    if "PYTHONPATH" in env:
        env["PYTHONPATH"] = (
            os.pathsep.join(extra_paths) + os.pathsep + env["PYTHONPATH"]
        )
    else:
        env["PYTHONPATH"] = os.pathsep.join(extra_paths)
    return env


def _validate_single_ipynb(nb_path):
    """Validate JSON syntax and nbformat v4 schema for a single .ipynb file."""
    nb_name = os.path.basename(nb_path)
    print(f"Validating {nb_name}...")

    try:
        with open(nb_path, "r", encoding="utf-8") as handle:
            json.load(handle)
    except json.JSONDecodeError as exc:
        print(
            f"  [ERROR] Invalid JSON syntax in {nb_name}: {exc}",
            file=sys.stderr,
        )
        return False

    try:
        with open(nb_path, "r", encoding="utf-8") as handle:
            nb = nbformat.read(handle, as_version=4)
        validate(nb)
        print(f"  {nb_name} is valid.")
        return True
    except NotebookValidationError as exc:
        print(
            f"  [ERROR] Schema validation failed in {nb_name}: {exc}",
            file=sys.stderr,
        )
        return False
    except Exception as exc:  # pylint: disable=broad-exception-caught
        print(
            f"  [ERROR] Unexpected error validating {nb_name}: {exc}",
            file=sys.stderr,
        )
        return False


def _validate_jupyter_notebooks(script_dir, env):
    """Build and validate all generated .ipynb notebooks in script_dir."""
    print("Building notebooks...")
    try:
        subprocess.run(
            [sys.executable, "build_pyclical_notebooks.py"],
            cwd=script_dir,
            env=env,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        print(f"Failed to build notebooks: {exc}", file=sys.stderr)
        return False

    notebook_files = sorted(glob.glob(os.path.join(script_dir, "*.ipynb")))
    if not notebook_files:
        print("No notebook files found to validate.", file=sys.stderr)
        return False

    results = [_validate_single_ipynb(nb_path) for nb_path in notebook_files]
    return all(results)


def _check_marimo_syntax(marimo_files):
    """Compile-check generated *_marimo.py files for Python syntax errors."""
    has_errors = False
    for nb_path in marimo_files:
        nb_name = os.path.basename(nb_path)
        try:
            with open(nb_path, "r", encoding="utf-8") as handle:
                source = handle.read()
            compile(source, nb_path, "exec")
        except (OSError, SyntaxError) as exc:
            print(
                f"  [ERROR] Python syntax check failed in {nb_name}: {exc}",
                file=sys.stderr,
            )
            has_errors = True
    return not has_errors


def _run_marimo_notebook(nb_path):
    """Load and execute a Marimo notebook headlessly and verify session init."""
    base_name = os.path.splitext(os.path.basename(nb_path))[0]
    spec = importlib.util.spec_from_file_location(
        f"_validate_{base_name}", nb_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot create module spec for {nb_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    _, defs = mod.app.run()
    session = defs.get("session")
    if session is None:
        raise RuntimeError("Missing 'session' definition in notebook")
    if session.init_error is not None:
        raise RuntimeError(
            f"PyClical session initialization failed: {session.init_error}"
        )


def _validate_marimo_notebooks(script_dir, env):
    """Build, lint-check, and headlessly execute all *_marimo.py notebooks."""
    print("Building Marimo notebooks...")
    try:
        subprocess.run(
            [sys.executable, "build_marimo_notebooks.py"],
            cwd=script_dir,
            env=env,
            check=True,
        )
    except subprocess.CalledProcessError as exc:
        print(f"Failed to build Marimo notebooks: {exc}", file=sys.stderr)
        return False

    marimo_files = sorted(glob.glob(os.path.join(script_dir, "*_marimo.py")))
    if not marimo_files:
        print("No Marimo notebook files found to validate.", file=sys.stderr)
        return False

    has_errors = not _check_marimo_syntax(marimo_files)

    if importlib.util.find_spec("marimo") is None:
        print(
            "[INFO] 'marimo' is not installed in this Python environment; "
            "skipping marimo check and execution."
        )
        return not has_errors

    print("Checking Marimo notebooks...")
    check_res = subprocess.run(
        [sys.executable, "-m", "marimo", "check", "--strict", *marimo_files],
        cwd=script_dir,
        env=env,
        check=False,
        capture_output=True,
        text=True,
    )
    if check_res.returncode != 0:
        if check_res.stdout:
            print(check_res.stdout, file=sys.stderr)
        if check_res.stderr:
            print(check_res.stderr, file=sys.stderr)
        print("[ERROR] marimo check --strict failed.", file=sys.stderr)
        has_errors = True

    for nb_path in marimo_files:
        nb_name = os.path.basename(nb_path)
        print(f"Validating {nb_name}...")
        try:
            _run_marimo_notebook(nb_path)
            print(f"  {nb_name} is valid.")
        except Exception as exc:  # pylint: disable=broad-exception-caught
            print(
                f"  [ERROR] Execution failed in {nb_name}: {exc}",
                file=sys.stderr,
            )
            has_errors = True

    return not has_errors


def main():
    """Main execution entry point for building and validating notebooks."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    env = _build_python_env(script_dir)
    jupyter_ok = _validate_jupyter_notebooks(script_dir, env)
    marimo_ok = _validate_marimo_notebooks(script_dir, env)
    return 0 if (jupyter_ok and marimo_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
