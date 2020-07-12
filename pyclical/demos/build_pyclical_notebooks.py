#!/usr/bin/env python3
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

import pyclical_tutorials as pt
import pyclical_notebook_utils as pn
import sys

def build_notebook(ctx, module_name, title):
    module = __import__(module_name)
    sys.stdout = open(module_name + ".ipynb", "w")
    ctx.print_notebook_header(title)
    module.run(ctx)
    ctx.print_notebook_footer()

def build_notebook_from_demo(ctx, module_name):
    build_notebook(ctx, module_name, module_name)

def build_notebook_from_tut(ctx, tut):
    build_notebook(ctx, tut.module_name, tut.title)

ctx = pn.notebook_context(globals())
build_notebook_from_demo(ctx, "clifford_demo")
build_notebook_from_demo(ctx, "pyclical_demo")
build_notebook_from_demo(ctx, "sqrt_log_demo")

for key, tut in sorted(pt.tutorial_dict.items()):
    build_notebook_from_tut(ctx, tut)
