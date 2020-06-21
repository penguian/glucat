# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_1_4_fourth.py:
#
# This file is based on the Plane Geometry section of the Tutorial from
# [LMV] Pertti Lounesto, R. Mikkola, V. Vierros,
# CLICAL Version 2.0 User Manual: Complex Number, Vector Space and
# Clifford Algebra Calculator for MS-DOS Personal Computers,
# Helsinki University of Technology Institute of Mathematics
# Research Reports A248, August 1987.
#
# Portions of [LMV] are reproduced here by permission of Aalto University, 2012.
#
#    copyright            : (C) 2012-2014 by Paul C. Leopardi
#
# Licensed under CC BY-SA 3.0 http://creativecommons.org/licenses/by-sa/3.0/

from pyclical_tutorial_utils import *

def run(ctx):
    for name, method in get_object_methods(ctx).items():
        exec("global "+name+";"+name+"=method")

    print_fill("1.4 The fourth dimension")
    print_line()
    print_fill("This tutorial file contains examples which will introduce you to PyClical" +
              " and the wide range of calculations with Clifford and Grassmann algebras that" +
              " you can use PyClical to perform.")
    print_line()
    print_fill("It is recommended that you do the tutorials in order, beginning with" +
              " 0.0 Notation.")
    print_line()

    pause()
    print_line()
    print_fill("This file is based on the Plane Geometry section of the Tutorial from")
    print_fill("[LMV] Pertti Lounesto, R. Mikkola, V. Vierros,")
    print_fill("CLICAL Version 2.0 User Manual: Complex Number, Vector Space and" +
              " Clifford Algebra Calculator for MS-DOS Personal Computers,")
    print_fill("Helsinki University of Technology Institute of Mathematics" +
              " Research Reports A248, August 1987.")
    print_line()
    print_fill("Example numbers refer to [LMV], e.g. Example 10.2 is the second example" +
              " on p. 10 of [LMV].")
    print_line()
    print_exec("from PyClical import *")

    pause()
    print_line()
    print_fill("Example 15.1. Calculate the length of the vector a == 2{1}+2{2}+4{3}+5{4}.")
    print_line()
    print_exec("a = 2*e(1)+2*e(2)+4*e(3)+5*e(4); print(a)")
    print_exec("length = abs(a); print(length)")

    pause()
    print_line()
    print_fill("Example 15.2. Determine the area of the parallelogram with sides" +
              " a == 2{1}+2{2}+4{3}+5{4} and b == {1}+2{2}+2{3}+7{4}.")
    print_line()
    print_exec("a = 2*e(1)+2*e(2)+4*e(3)+5*e(4); print(a)")
    print_exec("b = e(1)+2*e(2)+2*e(3)+7*e(4); print(b)")
    print_exec("area = abs(a ^ b); print(area)")

    pause()
    print_line()
    print_fill("Example 15.3. Determine the volume of the parallelepiped with sides" +
              " a == 2{1}, b == 3{2}+{4} and c == {1}+4{3}+2{4}.")
    print_line()
    print_exec("a = 2*e(1); print(a)")
    print_exec("b = 3*e(2)+e(4); print(b)")
    print_exec("c = e(1)+4*e(3)+2*e(4); print(c)")

    pause()
    check_exec("set V to be the oriented volume of the parallelepiped with sides a, b and c.",
               "V",
               "a ^ b ^ c")
    check_exec("set vol to be the absolute volume of the parallelepiped with sides a, b and c.",
               "vol",
               "abs(V)")

    pause()
    print_line()
    print_fill("Example 15.4.Rotate the vector x == {1}+{2}+{4} in the four-dimensional Euclidean space" +
              " first in the {1,2} plane by the angle pi/3 to produce the vector y, and then in the" +
              " {3,4} plane by the angle pi/4 to produce the vector z.")
    print_line()
    print_exec("x = e(1)+e(3)+e(4); print(x)")
    print_exec("A = pi/3*e({1,2}); print(A)")
    print_exec("a = exp(-A/2); print(a)")
    print_exec("y = x | a; print(y)")
    print_line()
    print_fill("Alternatively we could have:")
    print_line()
    print_exec("a1 = exp(A/2); print(a1)")
    print_exec("y1 = inv(a1) * x * a1; print(y1)")

    pause()
    check_exec("set B to be the bivector that generates a rotation in the {3,4} plane by the angle pi/4.",
               "B",
               "pi/4*e({3,4})")

    pause()
    print_exec("b = exp(-B/2); print(b)")
    print_exec("z = y | b; print(z)")
    print_line()
    print_fill("Compare with Example 12.3 from Tutorial 1.2 Space geometry and vector algebra.")
    print_line()

    pause()
    print_line()
    print_fill("Example 16.1. How far is the point P == (1,-5,5,6) from the plane B spanned by the vectors" +
              " x == 3{1}-6{2}+7{3}+2{4} and y == 3{1}+2{2}+3{3}-2{4} ?")
    print_line()
    print_exec("R4 = istpq(4,0); print(R4)")
    print_exec("OP = clifford((1,-5, 5, 6), R4); print(OP)")
    print_exec("x = clifford((3,-6, 7, 2), R4); print(x)")
    print_exec("y = clifford((3, 2, 3,-2), R4); print(y)")
    print_exec("B = x ^ y; print(B)")
    print_exec("v = (OP ^ B) / B; print(v)")
    print_exec("dist = abs(v); print(dist)")

    print_line()
    print_fill(" In this example we have used the alternate notation for vectors.")
    print_line()

    pause()
    print_line()
    print_fill("Example 16.2. Find the squares of the multivectors (1 + {1,2} + {3,4} +/- {1,2,3,4})/2.")
    print_line()
    print_exec("a = (1+e({1,2})+e({3,4})+e({1,2,3,4}))/2; print(a)")
    print_exec("a2 = a ** 2; print(a2)")

    pause()
    check_exec("set b to be the multivector (1 + {1,2} + {3,4} - {1,2,3,4})/2.",
               "b",
               "a - e({1,2,3,4})")
    check_exec("set b2 to be the square of b.",
               "b2",
               "b * b")

    pause()
    print_line()
    print_fill("You have completed the tutorial file pyclical_tutorial_1_4_fourth.py.")

if __name__ == "__main__":
    ctx = tutorial_context(globals())
    try:
        run(ctx)
    except:
        ctx.print_fill("The tutorial was interrupted.")
        pass
