# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_1_0_plane.py:
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
    pause      = ctx.pause
    print_line = ctx.print_line
    print_fill = ctx.print_fill
    print_exec = ctx.print_exec
    check_exec = ctx.check_exec
    check_eval = ctx.check_eval

    print_fill("# pyclical_tutorial_1_0_plane.py:")
    print_line()
    print_fill("1.0 Plane Geometry")
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
    print_fill("Example 7.1. Divide the vector r == 5{1}-{2} into components parallel to" +
              " a == {1}-2{2} and b == {1}+{2}, that is, determine the coefficients in" +
              " the decomposition r == alpha*a + beta*b.")
    print_line()
    print_exec("a = 3*e(1)-2*e(2); print(a)")
    print_exec("b = e(1)+e(2); print(b)")
    print_exec("r = 5*e(1)-e(2); print(r)")
    print_exec("alpha = (r ^ b) / (a ^ b); print(alpha)")

    pause()
    check_exec("set beta to be the coefficient of r with respect to b.",
               "beta",
               "(a ^ r) / (a ^ b)")
    pause()
    check_eval("compare with r to check that the values\n" "of alpha and beta are correct.",
               "alpha * a + beta * b",
               "print(abs(r - ({})))")

    pause()
    print_line()
    print_fill(" In this example a ^ b is the oriented area of the parallelogram with" +
              "  sides a and b.")
    print_line()

    pause()
    print_line()
    print_fill("Example 7.2. Compute the angle between the directions a == 3{1}+{2} and" +
              " b == 4{1}-2{2}.")
    print_line()
    print_exec("a = 3*e(1)+e(2); print(a)")
    print_exec("b = 4*e(1)-2*e(2); print(b)")
    print_exec("c = a & b")
    print_exec("d = abs(a * b)")

    pause()
    check_exec("set theta to be the angle between a and b," +
               "computed using the inverse cosine function and the values of c and d.",
               "theta",
               "acos(c / d)")
    pause()
    check_eval("compare with scalar(a & b) to check that the value of theta is correct.",
               "abs(a) * abs(b) * cos(theta)",
               "print(abs(scalar(a & b) - {}))")

    pause()
    print_line()
    print_fill("The value of theta is expressed in radians and corresponds to 45 degrees. " +
              " This can be checked as follows:")
    print_line()
    print_exec("print(theta/pi)")

    pause()
    print_line()
    print_fill("Example 8.1. Find the projection of the vector a == 8{1}-{2} in the" +
              " direction b == 2{1}+{2}, and the component of a perpendicular to b.")
    print_line()
    print_exec("a = 8*e(1)-e(2); print(a)")
    print_exec("b = 2*e(1)+e(2); print(b)")

    pause()
    check_exec("set c to be the projection of a in the direction b.",
               "c",
               "(a & b) / b")
    pause()
    check_exec("set d to be the component of a perpendicular to b.",
               "d",
               "(a ^ b) / b")
    pause()
    check_eval("compare with a to check that the values of c and d are correct.",
               "c + d",
               "print(abs(a - (c + d)))")

    pause()
    print_line()
    print_fill("Example 8.2. Rotate the vector r==2{1}+3{2} by the angle pi/3==60 degrees.")
    print_line()
    print_exec("r = 2*e(1)+3*e(2); print(r)")
    print_exec("i = e({1,2}); print(i)")
    print_exec("z = exp(i * pi / 3); print(z)")

    pause()
    check_exec("set rz to be the vector r rotated by the angle pi/3.",
               "rz",
               "r * z")
    pause()
    check_exec("set zr to be the vector r rotated by the angle -pi/3.",
               "zr",
               "z * r")

    pause()
    print_line()
    print_fill("Note that the rotations r |-> r*z and r |-> z*r have opposite orientations. " +
              " This is because the unit bivector i=={1,2} anticommutes with all vectors in" +
              " the {1,2} plane.")
    print_line()

    pause()
    print_line()
    print_fill("Example 8.3. Combine two successive reflections of the vector r==4{1}-3{2}" +
              " in the lines defined by the vectors a==3{1}-{2} and b==2{1}+{2}.")
    print_line()
    print_exec("r = 4*e(1)-3*e(2); print(r)")
    print_exec("a = 3*e(1)-e(2); print(a)")
    print_exec("b = 2*e(1)+e(2); print(b)")

    pause()
    check_exec("set q to be the vector r reflected in the line defined by the vector a.",
               "q",
               "r | a")
    pause()
    check_exec("set p to be the vector q reflected in the line defined by the vector b.",
               "p",
               "q | b")

    pause()
    print_line()
    print_fill("The notation r | a means a * r / involute(a), and if a is a vector," +
              " this equals -a*r/a == (r*a - 2*(a & r)) / a == r - 2*(a & r)*a / (a & a).")
    print_line()

    pause()
    check_eval("compare with r | a to check the meaning of r | a.",
               "a * r / involute(a)",
               "print(abs( (r | a) - ({}) ))")
    pause()
    check_eval("compare with r | a to check the first equality above.",
               "-a * r / a",
               "print(abs( (r | a) - ({}) ))")
    pause()
    check_eval("compare with r | a to check the second equality above.",
               "(r * a - 2 * (a & r)) / a",
               "print(abs( (r | a) - ({}) ))")
    pause()
    check_eval("compare with r | a to check the third\n" "equality above.",
               "r - 2 * (a & r) * a / (a & a)",
               "print(abs( (r | a) - ({}) ))")

    pause()
    print_line()
    print_fill("You have completed the tutorial file pyclical_tutorial_1_0_plane.py.")

if __name__ == "__main__":
    ctx = tutorial_context(globals())
    try:
        run(ctx)
    except:
        ctx.print_fill("The tutorial was interrupted.")
        pass
