#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# verify_all.py: Run all verification tests and checks.
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

"""
Verify all GluCat and PyClical test suites, lints, and examples.
"""

import argparse
import os
import shutil
import subprocess
import sys


def find_nbformat_python():
    """Find a Python interpreter executable that has `nbformat` installed."""
    candidates = [sys.executable]
    path_python = shutil.which("python3")
    if path_python and path_python not in candidates:
        candidates.append(path_python)
    for std_path in ("/usr/bin/python3", "/usr/local/bin/python3"):
        if os.path.exists(std_path) and std_path not in candidates:
            candidates.append(std_path)

    for py_bin in candidates:
        try:
            res = subprocess.run(
                [py_bin, "-c", "import nbformat"],
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            if res.returncode == 0:
                return py_bin
        except Exception:
            continue
    return None


def find_cython_python():
    """Find a Python interpreter executable that has `Cython` installed."""
    candidates = [sys.executable]
    path_python = shutil.which("python3")
    if path_python and path_python not in candidates:
        candidates.append(path_python)
    for std_path in ("/usr/bin/python3", "/usr/local/bin/python3"):
        if os.path.exists(std_path) and std_path not in candidates:
            candidates.append(std_path)

    for py_bin in candidates:
        try:
            res = subprocess.run(
                [py_bin, "-c", "import Cython"],
                check=False,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            if res.returncode == 0:
                return py_bin
        except Exception:
            continue
    return sys.executable


def check_makefile_target(cmd, cwd=None):
    """Verify that a requested make target exists using dry-run (`make -n`)."""
    if not cmd or cmd[0] != "make":
        return
    dry_run_cmd = ["make", "-n"] + cmd[1:]
    res = subprocess.run(
        dry_run_cmd,
        cwd=cwd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if res.returncode != 0:
        print(
            f"Error: Target in command '{' '.join(cmd)}' does not exist.",
            file=sys.stderr,
        )
        sys.exit(1)


def run_cmd(cmd, cwd=None, env=None):
    """Execute a subprocess command, exiting on failure."""
    if cmd and cmd[0] == "make":
        check_makefile_target(cmd, cwd=cwd)
    print(f"Running: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, cwd=cwd, env=env, check=True, stdin=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e}", file=sys.stderr)
        sys.exit(e.returncode)


def main():
    """Main verification runner parsing flags and executing checks."""
    os.environ.setdefault("MAKEFLAGS", "-j4")
    python_bin = find_cython_python()
    parser = argparse.ArgumentParser(
        description="GluCat & PyClical verification runner"
    )
    parser.add_argument(
        "--coverage",
        action="store_true",
        help="Run C++ block and branch coverage tests",
    )
    parser.add_argument(
        "--examples",
        action="store_true",
        help="Run Python tests, notebook validation, and examples",
    )
    args, extra_args = parser.parse_known_args()

    if extra_args and not args.coverage:
        parser.error(f"unrecognized arguments: {' '.join(extra_args)}")

    # Get the project root directory
    root_dir = os.path.dirname(os.path.abspath(__file__))

    # Run core verification targets
    print("--- Core build check ---")
    run_cmd(["make", "fast-check"], cwd=root_dir)

    if args.coverage:
        print("--- C++ header coverage check ---")
        coverage_cmd = ["make", "check-coverage-doctest"]
        if extra_args:
            coverage_cmd.extend(extra_args)
        run_cmd(coverage_cmd, cwd=root_dir, env=os.environ.copy())

    if args.examples:
        print("--- License headers check ---")
        run_cmd([python_bin, "check_license_headers.py"], cwd=root_dir)

        print("--- Ruff check ---")
        run_cmd(
            [
                "ruff",
                "check",
                "verify_all.py",
                "check_license_headers.py",
                "pyclical/",
            ],
            cwd=root_dir,
        )
        print("--- Pylint check ---")
        run_cmd(
            ["pylint", "pyclical/", "pyclical/demos/"],
            cwd=root_dir,
        )

        # Run PyClical C++-based tests
        print("--- PyClical test ---")
        run_cmd(
            ["make", "-C", "pyclical", "check", f"PYTHON={python_bin}"],
            cwd=root_dir,
        )

        # Run notebook validation
        print("--- Notebook validation ---")
        nbformat_python = find_nbformat_python()
        if nbformat_python:
            run_cmd(
                [nbformat_python, "pyclical/demos/validate_notebooks.py"],
                cwd=root_dir,
            )
        else:
            print(
                "[WARNING] 'nbformat' is not installed in any available "
                "Python interpreter."
            )
            print("[WARNING] Skipping notebook validation.")

        # Run demos in pyclical/demos
        print("--- Demos check ---")
        demos_dir = os.path.join(root_dir, "pyclical", "demos")
        demo_files = [
            "clifford_demo.py",
            "m_theory_demo.py",
            "sqrt_log_demo.py",
            "versor_demo.py",
        ]

        # Include PYTHONPATH to find PyClical in the built/inplace dir and
        # set GLUCAT_NON_INTERACTIVE
        env = os.environ.copy()
        env["GLUCAT_NON_INTERACTIVE"] = "1"
        pyclical_path = os.path.join(root_dir, "pyclical")
        if "PYTHONPATH" in env:
            env["PYTHONPATH"] = pyclical_path + os.pathsep + env["PYTHONPATH"]
        else:
            env["PYTHONPATH"] = pyclical_path

        for demo in demo_files:
            print(f"Running demo: {demo}")
            run_cmd([python_bin, demo], cwd=demos_dir, env=env)

        # Run plotting demo in non-interactive headless mode if pyvista is available
        try:
            subprocess.run(
                [python_bin, "-c", "import pyvista"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            print("Running plotting demo: plotting_demo_pyvista.py")
            plotting_dir = os.path.join(root_dir, "pyclical", "demos", "plotting")
            run_cmd(
                [python_bin, "plotting_demo_pyvista.py"],
                cwd=plotting_dir,
                env=env,
            )
        except Exception:
            print(
                "PyVista not found in Python environment; skipping "
                "plotting_demo_pyvista.py."
            )

        print("=== Python and examples validation succeeded! ===")


if __name__ == "__main__":
    main()
