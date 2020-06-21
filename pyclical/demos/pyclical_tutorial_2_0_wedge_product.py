# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_2_0_wedge_product.py:
#
# This file contains a tutorial that demonstrates one property of the exterior product.
#
#    copyright            : (C) 2013-2014 by Paul C. Leopardi
#
# Licensed under CC BY-SA 3.0 http://creativecommons.org/licenses/by-sa/3.0/

from pyclical_tutorial_utils import *

def run(ctx):
    for name, method in get_object_methods(ctx).items():
        exec("global "+name+";"+name+"=method")

    print_head("2.0 Exterior product.")
    print_line()
    print_fill("This file contains a tutorial that demonstrates one property of" +
              " the exterior (wedge) product.")
    print_line()
    print_fill("It is recommended that you do the tutorials in order, beginning with" +
              " 0.0 Notation.")
    pause()
    print_line()
    print_exec("from PyClical import *")

    print_line()
    print_fill("Here we are working in R_{4,0}, but the results do not depend on signature.")
    print_line()
    print_exec("frame = istpq(4,0)")

    print_line()
    print_fill("First, we create two random vectors, u and v.")
    print_line()
    print_exec("u = random_clifford(frame)(1); print(u)")
    print_exec("v = random_clifford(frame)(1); print(v)")

    print_line()
    print_fill("We now set the vector w to be a random linear combination of u and v.")
    print_line()
    print_exec("a = random_clifford(ist(0)); print(a)")
    print_exec("b = random_clifford(ist(0)); print(b)")
    print_exec("w = a*u + b*v; print(w)")

    print_line()
    print_fill("The exterior product shows the pairwise independence of u, v and w ...")
    print_line()
    print_exec("print(u^v)")
    print_exec("print(u^w)")
    print_exec("print(v^w)")

    print_line()
    print_fill("... and the linear dependence of the set {u, v, w}.")
    print_line()
    print_exec("print(u^v^w)")

    pause()
    print_line()
    print_fill("You have completed the tutorial file pyclical_tutorial_2_0_wedge_product.py.")

if __name__ == "__main__":
    ctx = tutorial_context(globals())
    try:
        run(ctx)
    except:
        ctx.print_fill("The tutorial was interrupted.")
        pass
