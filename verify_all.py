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
    """
Find a Python interpreter executable that has `nbformat` installed.
"""
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
    """
Find a Python interpreter executable that has `Cython` installed.
"""
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


def get_clean_make_env(env=None):
    """
    Sanitize MAKEFLAGS in env to strip jobserver flags that cause
    warnings in sub-make processes.
    """
    base_env = env.copy() if env is not None else os.environ.copy()
    if "MAKEFLAGS" in base_env:
        flags = base_env["MAKEFLAGS"].split()
        cleaned = [f for f in flags if not f.startswith("--jobserver")]
        base_env["MAKEFLAGS"] = " ".join(cleaned)
    return base_env


def check_makefile_target(cmd, cwd=None):
    """
Verify that a requested make target exists using dry-run (`make -n`).
"""
    if not cmd or cmd[0] != "make":
        return
    dry_run_cmd = ["make", "-n"] + cmd[1:]
    res = subprocess.run(
        dry_run_cmd,
        cwd=cwd,
        env=get_clean_make_env(),
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


def run_cmd(cmd, cwd=None, env=None, quiet=False):
    """
Execute a subprocess command, exiting on failure.
"""
    if cmd and cmd[0] == "make":
        check_makefile_target(cmd, cwd=cwd)
        env = get_clean_make_env(env)
    if not quiet:
        print(f"Running: {' '.join(cmd)}")
    try:
        if quiet:
            subprocess.run(
                cmd,
                cwd=cwd,
                env=env,
                check=True,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
        else:
            subprocess.run(
                cmd, cwd=cwd, env=env, check=True, stdin=subprocess.DEVNULL
            )
    except subprocess.CalledProcessError as e:
        if quiet and e.output:
            print(e.output.decode("utf-8", errors="replace"), file=sys.stderr)
        print(f"Command failed: {e}", file=sys.stderr)
        sys.exit(e.returncode)


def run_all_demos(root_dir, python_bin, quiet=False):
    """
Run interactive and plotting python demos in pyclical/demos.
"""
    demos_dir = os.path.join(root_dir, "pyclical", "demos")
    demo_files = [
        "clifford_demo.py",
        "m_theory_demo.py",
        "sqrt_log_demo.py",
        "versor_demo.py",
    ]

    env = os.environ.copy()
    env["GLUCAT_NON_INTERACTIVE"] = "1"
    pyclical_path = os.path.join(root_dir, "pyclical")
    if "PYTHONPATH" in env:
        env["PYTHONPATH"] = pyclical_path + os.pathsep + env["PYTHONPATH"]
    else:
        env["PYTHONPATH"] = pyclical_path

    for demo in demo_files:
        if not quiet:
            print(f"Running demo: {demo}")
        run_cmd([python_bin, demo], cwd=demos_dir, env=env, quiet=quiet)

    try:
        subprocess.run(
            [python_bin, "-c", "import pyvista"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if not quiet:
            print("Running plotting demo: plotting_demo_pyvista.py")
        plotting_dir = os.path.join(root_dir, "pyclical", "demos", "plotting")
        run_cmd(
            [python_bin, "plotting_demo_pyvista.py"],
            cwd=plotting_dir,
            env=env,
            quiet=quiet,
        )
    except Exception:
        if not quiet:
            print(
                "PyVista not found in Python environment; skipping "
                "plotting_demo_pyvista.py."
            )


def main():
    """
Main verification runner parsing flags and executing checks.
"""
    os.environ.setdefault("MAKEFLAGS", "-j4")
    python_bin = find_cython_python()
    parser = argparse.ArgumentParser(
        description="GluCat & PyClical verification runner"
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress verbose command output unless a step fails",
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help="Run ultra-fast checks (core build, license headers, Ruff linter)",
    )
    parser.add_argument(
        "--coverage",
        action="store_true",
        help="Run C++ block and branch coverage tests",
    )
    parser.add_argument(
        "--python",
        "--examples",
        dest="python",
        action="store_true",
        help="Run Python tests, notebook validation, and examples",
    )
    args, extra_args = parser.parse_known_args()

    if extra_args and not args.coverage:
        parser.error(f"unrecognized arguments: {' '.join(extra_args)}")

    root_dir = os.path.dirname(os.path.abspath(__file__))

    def log_step(name):
        if not args.quiet:
            print(f"--- {name} ---")

    def log_success(name):
        if args.quiet:
            print(f"{name}... OK")

    log_step("Build target pre-flight check")
    check_makefile_target(["make", "check"], cwd=root_dir)
    log_success("Build target pre-flight check")

    if args.fast or args.python:
        log_step("License headers check")
        run_cmd(
            [python_bin, "check_license_headers.py"],
            cwd=root_dir,
            quiet=args.quiet,
        )
        log_success("License headers check")

        log_step("Ruff check")
        run_cmd(
            [
                "ruff",
                "check",
                "verify_all.py",
                "check_license_headers.py",
                "pyclical/",
            ],
            cwd=root_dir,
            quiet=args.quiet,
        )
        log_success("Ruff check")

    if args.coverage:
        log_step("C++ header coverage check")
        coverage_cmd = ["make", "check-coverage-doctest"]
        if extra_args:
            coverage_cmd.extend(extra_args)
        run_cmd(coverage_cmd, cwd=root_dir, env=os.environ.copy(), quiet=args.quiet)
        log_success("C++ header coverage check")

    if args.python:
        log_step("Pylint check")
        run_cmd(
            ["pylint", "pyclical/", "pyclical/demos/"],
            cwd=root_dir,
            quiet=args.quiet,
        )
        log_success("Pylint check")

        log_step("PyClical test (pytest)")
        run_cmd(
            [
                python_bin,
                "-m",
                "pytest",
                "pyclical/test_pytest_doctests.py",
                "-v",
            ],
            cwd=root_dir,
            quiet=args.quiet,
        )
        log_success("PyClical test (pytest)")

        log_step("Notebook validation")
        nbformat_python = find_nbformat_python()
        if nbformat_python:
            run_cmd(
                [nbformat_python, "pyclical/demos/validate_notebooks.py"],
                cwd=root_dir,
                quiet=args.quiet,
            )
            log_success("Notebook validation")
        else:
            print(
                "[WARNING] 'nbformat' is not installed in any available "
                "Python interpreter."
            )
            print("[WARNING] Skipping notebook validation.")

        log_step("Demos check")
        run_all_demos(root_dir, python_bin, quiet=args.quiet)
        log_success("Demos check")

    print("=== Python and examples validation succeeded! ===")


if __name__ == "__main__":
    main()
