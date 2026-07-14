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

import argparse
import os
import shutil
import subprocess
import sys


def run_cmd(cmd, cwd=None, env=None):
    print(f"Running: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, cwd=cwd, env=env, check=True, stdin=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e}", file=sys.stderr)
        sys.exit(e.returncode)


def main():
    parser = argparse.ArgumentParser(description="GluCat & PyClical verification runner")
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

    # Case 1: Coverage
    if args.coverage:
        print("=== Running C++ Block Coverage ===")
        # Check if llvm tools exist
        for tool in ["llvm-profdata", "llvm-cov"]:
            if shutil.which(tool) is None:
                print(f"Error: {tool} not found. Please install LLVM tools.", file=sys.stderr)
                sys.exit(1)
        run_cmd(
            ["bash", "test_coverage/src/run_clang_doctest_coverage.sh"] + extra_args,
            cwd=root_dir,
        )
        return

    # Check if Makefile exists
    if not os.path.isfile(os.path.join(root_dir, "Makefile")):
        print(
            "Error: Makefile not found. Please bootstrap and configure the project first:\n"
            "  make -f admin/Makefile.common bootstrap\n"
            "  ./configure [options]",
            file=sys.stderr,
        )
        sys.exit(1)

    # Case 2: Examples and Python verification only
    if args.examples:
        print("=== Running Python verification and examples ===")

        # Run python lints
        print("--- Ruff Lint check ---")
        run_cmd(
            ["ruff", "check", "pyclical/", "pyclical/demos/", "benchmarks/"],
            cwd=root_dir,
        )
        print("--- Pylint check ---")
        run_cmd(
            ["pylint", "pyclical/", "pyclical/demos/", "benchmarks/"],
            cwd=root_dir,
        )

        # Run PyClical C++-based tests
        print("--- PyClical test ---")
        run_cmd(["make", "-C", "pyclical", "check"], cwd=root_dir)

        # Run notebook validation
        print("--- Notebook validation ---")
        run_cmd(
            [sys.executable, "pyclical/demos/validate_notebooks.py"],
            cwd=root_dir,
        )

        # Run demos in pyclical/demos
        print("--- Demos check ---")
        demos_dir = os.path.join(root_dir, "pyclical", "demos")
        demo_files = [
            "clifford_demo.py",
            "m_theory_demo.py",
            "sqrt_log_demo.py",
            "versor_demo.py",
        ]

        # Include PYTHONPATH to find PyClical in the built/inplace dir and set GLUCAT_NON_INTERACTIVE
        env = os.environ.copy()
        env["GLUCAT_NON_INTERACTIVE"] = "1"
        pyclical_path = os.path.join(root_dir, "pyclical")
        if "PYTHONPATH" in env:
            env["PYTHONPATH"] = pyclical_path + os.pathsep + env["PYTHONPATH"]
        else:
            env["PYTHONPATH"] = pyclical_path

        for demo in demo_files:
            print(f"Running demo: {demo}")
            run_cmd([sys.executable, demo], cwd=demos_dir, env=env)

        # Run plotting demo in non-interactive headless mode
        print("Running plotting demo: plotting_demo_pyvista.py")
        plotting_dir = os.path.join(root_dir, "pyclical", "demos", "plotting")
        run_cmd([sys.executable, "plotting_demo_pyvista.py"], cwd=plotting_dir, env=env)

        print("=== Python and examples validation succeeded! ===")
        return

    # Case 3: Default run (all C++ checks + Python checks)
    print("=== Running all GluCat / PyClical verification tests ===")

    # 1. C++ linting & static analysis
    print("--- C++ linting ---")
    run_cmd(["make", "lint-check"], cwd=root_dir)
    print("--- Cppcheck ---")
    run_cmd(["make", "cppcheck"], cwd=root_dir)

    # 2. Primary C++ regression suite
    print("--- C++ regression tests ---")
    run_cmd(["make", "check-local"], cwd=root_dir)

    # 3. Python lints, tests, and notebook validation
    print("--- Python lints & notebook validation ---")
    run_cmd(["ruff", "check", "pyclical/", "pyclical/demos/", "benchmarks/"], cwd=root_dir)
    run_cmd(["pylint", "pyclical/", "pyclical/demos/", "benchmarks/"], cwd=root_dir)
    run_cmd(["make", "-C", "pyclical", "check"], cwd=root_dir)
    run_cmd([sys.executable, "pyclical/demos/validate_notebooks.py"], cwd=root_dir)

    print("=== All checks passed successfully! ===")


if __name__ == "__main__":
    main()
