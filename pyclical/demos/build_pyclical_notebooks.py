#!/usr/bin/env python3
"""Build IPython notebook files for PyClical tutorials and demos."""
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# build_pyclical_notebooks.py:
#
# This file contains a script that builds IPython notebooks corresponding to
# PyClical tutorials and demos.
#
#    copyright            : (C) 2014 by Paul C. Leopardi
#
# Licensed under CC BY-SA 3.0 http://creativecommons.org/licenses/by-sa/3.0/

import sys
import pyclical_notebook_utils as pn
import pyclical_tutorials as pt


def build_notebook(ctx, module_name, title):
    """Build a single IPython notebook from a PyClical module."""
    module = __import__(module_name)
    original_stdout = sys.stdout
    try:
        with open(module_name + ".ipynb", "w", encoding="utf-8") as f:
            sys.stdout = f
            ctx.print_notebook_header(title)
            module.run(ctx)
            ctx.print_notebook_footer()
    finally:
        sys.stdout = original_stdout


def build_notebook_from_demo(ctx, module_name):
    """Build an IPython notebook for a demo module."""
    build_notebook(ctx, module_name, module_name)


def build_notebook_from_tut(ctx, tut):
    """Build an IPython notebook for a tutorial module."""
    build_notebook(ctx, tut.module_name, tut.title)


ctx = pn.notebook_context(globals())
build_notebook_from_demo(ctx, "clifford_demo")
build_notebook_from_demo(ctx, "pyclical_demo")
build_notebook_from_demo(ctx, "sqrt_log_demo")
build_notebook_from_demo(ctx, "m_theory_demo")
build_notebook_from_demo(ctx, "versor_demo")

for key, tut in sorted(pt.tutorial_dict.items()):
    build_notebook_from_tut(ctx, tut)
