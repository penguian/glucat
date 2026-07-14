#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# versor_demo.py: Demonstrate versor and versor_exp with PyClical.
#
#    copyright            : (C) 2008-2026 by Paul C. Leopardi
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

from PyClical import *
from pyclical_tutorial_utils import *


def run(ctx=tutorial_context(globals())):
    for name, method in get_object_methods(ctx).items():
        exec("global " + name + ";" + name + "=method")

    print_fill(
        "# versor_demo.py: Demonstrate versor and versor_exp with PyClical."
    )
    print_line()
    print_fill("INITIALIZATION.")
    print_exec("from PyClical import *")
    print_line()

    print_fill("Define a bivector generator A and a multivector X.")
    print_line()
    print_exec("A = clifford('{1,2}') * pi/4")
    print_exec("print(A)")
    print_exec("X = clifford('{1}')")
    print_exec("print(X)")
    print_line()

    print_fill("Compute rotor R = exp(A).")
    print_line()
    print_exec("R = exp(A)")
    print_exec("print(R)")
    print_line()

    print_fill(
        "Compare the explicit sandwich product, operator|, versor(), and versor_exp()."
    )
    print_line()
    print_exec("explicit_sand = R * X * involute(exp(-A))")
    print_exec("print(explicit_sand)")
    print_exec("operator_vert = X | R")
    print_exec("print(operator_vert)")
    print_exec("versor_val = X.versor(R)")
    print_exec("print(versor_val)")
    print_exec("versor_exp_val = X.versor_exp(A)")
    print_exec("print(versor_exp_val)")
    print_exec("versor_exp_pre = X.versor_exp(A, prechecked=True)")
    print_exec("print(versor_exp_pre)")
    print_line()

    print_fill("Check differences (should be 0.0).")
    print_line()
    print_exec("print(abs(operator_vert - explicit_sand))")
    print_exec("print(abs(versor_val - explicit_sand))")
    print_exec("print(abs(versor_exp_val - explicit_sand))")
    print_exec("print(abs(versor_exp_pre - explicit_sand))")

    pause()
    print_line()
    print_fill("You have completed the demonstration file versor_demo.py.")


if __name__ == "__main__":
    try:
        run()
    except (KeyboardInterrupt, Exception):
        print("The demo was interrupted.")
        pass
