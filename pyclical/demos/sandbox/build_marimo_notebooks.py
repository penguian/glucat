#!/usr/bin/env python3
"""
Build Marimo reactive Python notebooks for PyClical tutorials and demos.
"""
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# build_marimo_notebooks.py:
#
# Build Marimo reactive notebooks corresponding to PyClical tutorials and demos.
#
#    copyright            : (C) 2026 by Paul C. Leopardi
#
# Licensed under CC BY-SA 3.0 http://creativecommons.org/licenses/by-sa/3.0/

import os
import sys

DEMO_MODULES = (
    "clifford_demo",
    "pyclical_demo",
    "sqrt_log_demo",
    "m_theory_demo",
    "versor_demo",
)


def _configure_paths(script_dir):
    """Ensure sandbox, demos, and pyclical directories are on sys.path."""
    demos_dir = os.path.dirname(script_dir)
    pyclical_dir = os.path.dirname(demos_dir)
    for path in (pyclical_dir, demos_dir, script_dir):
        if path not in sys.path:
            sys.path.insert(0, path)


def build_marimo_notebook(ctx, module_name, title, out_dir):
    """Build a single Marimo reactive notebook from a PyClical module."""
    module = __import__(module_name)
    out_path = os.path.join(out_dir, f"{module_name}_marimo.py")
    original_stdout = sys.stdout
    try:
        with open(out_path, "w", encoding="utf-8") as handle:
            sys.stdout = handle
            ctx.print_notebook_header(title, module_name=module_name)
            module.run(ctx)
            ctx.print_notebook_footer()
    finally:
        sys.stdout = original_stdout
    print(f"Generated {os.path.basename(out_path)}.")


def build_marimo_from_demo(ctx, module_name, out_dir):
    """Build a Marimo reactive notebook for a demo module."""
    build_marimo_notebook(ctx, module_name, module_name, out_dir)


def build_marimo_from_tut(ctx, tut, out_dir):
    """Build a Marimo reactive notebook for a tutorial module."""
    build_marimo_notebook(ctx, tut.module_name, tut.title, out_dir)


def main():
    """Generate all PyClical Marimo demo and tutorial notebooks."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    _configure_paths(script_dir)

    import pyclical_notebook_utils as pn  # pylint: disable=import-outside-toplevel
    import pyclical_tutorials as pt  # pylint: disable=import-outside-toplevel

    ctx = pn.marimo_context(globals())
    for demo_name in DEMO_MODULES:
        build_marimo_from_demo(ctx, demo_name, script_dir)

    for _, tut in sorted(pt.tutorial_dict.items()):
        build_marimo_from_tut(ctx, tut, script_dir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
