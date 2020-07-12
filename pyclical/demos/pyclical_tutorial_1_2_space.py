# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_1_2_space.py:
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

    print_head("1.2 Space geometry and vector algebra")
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
    print_fill("Example 10.2. Divide the vector r == -{1}+8{2}+{3} into components parallel to a == 2{1}-{2}," +
              " b == 2{1}+3{2}-{3} and c == 2{1}+2{2}+5{3}, that is, determine the coefficients alpha, beta" +
              " and gamma in the decomposition r == alpha*a + beta*b + gamma*c.")
    print_line()

    pause()
    print_exec("a = 2*e(1)-e(2); print(a)")
    print_exec("b = 2*e(1)+3*e(2)-e(3); print(b)")
    print_exec("c = 2*e(1)+2*e(2)+5*e(3); print(c)")
    print_exec("V = a ^ b ^ c; print(V)")
    print_exec("r = -e(1)+8*e(2)+e(3); print(r)")
    print_exec("alpha = (r ^ b ^ c) / V; print(alpha)")

    pause()
    check_exec("set beta to be the coefficient of r with respect to b.",
               "beta",
               "(a ^ r ^ c) / V")
    check_exec("set gamma to be the coefficient of r with respect to c.",
               "gamma",
               "(a ^ b ^ r) / V")
    check_eval("compare with r to check that the values of alpha, beta and gamma are correct.",
               "alpha * a + beta * b + gamma * c",
               "print(abs(r - ({})))")

    pause()
    print_line()
    print_fill("Compare with Example 7.1 from Tutorial 1.0 Plane geometry.")
    print_line()

    pause()
    print_line()
    print_fill("Example 11.1. A plane A is spanned by the vectors x == 4{1}+{3} and y == 3{1}+{2}. " +
              " Compute the projection p of the vector v == -5{1}+7{2} in the plane A and the" +
              " component r of v perpendicular to A.")
    print_line()

    pause()
    print_exec("v = -5*e(1)+7*e(2); print(v)")
    print_exec("x = 4*e(1)+e(3); print(x)")
    print_exec("y = 3*e(1)+e(2); print(y)")
    print_exec("A = x ^ y; print(A)")
    print_exec("p = (v & A) / A; print(p)")

    pause()
    check_exec("set r to be the component of v perpendicular to A.",
               "r",
               "(v ^ A) / A")
    check_eval("compare with v to check that the values of p and r are correct.",
               "p + r",
               "print(abs(v - ({})))")

    pause()
    print_line()
    print_fill("In this example the bivector A == x ^ y represents the oriented plane spanned by the" +
              " vectors x and y. " +
              " More precisely A is the oriented area of the parallelogram with sides x and y.")
    print_line()
    print_fill("Compare with Example 8.1 from Tutorial 1.0 Plane geometry.")
    print_line()

    pause()
    print_line()
    print_fill("Example 11.2. Determine three perpendicular unit vectors, t1, t2 and t3, with t1" +
              " parallel to v1 == 3{1}-{2}, t2 in the plane of v1 and v2 == {1}+2{3}, pointing into" +
              " the half-plane of v2, and t3 pointing into the half-space of v3 == {2}+{3}.")
    print_line()

    pause()
    print_exec("v1 = 3*e(1)-e(2); print(v1)")
    print_exec("v2 = e(1)+2*e(3); print(v2)")
    print_exec("v3 = e(2)+e(3); print(v3)")
    print_exec("A1 = v1; print(A1)")
    print_exec("A2 = A1 ^ v2; print(A2)")
    print_exec("A3 = A2 ^ v3; print(A3)")

    pause()
    print_exec("u1 = v1; print(u1)")
    print_exec("u2 = inv(A1) * A2; print(u2)")
    print_exec("u3 = inv(A2) * A3; print(u3)")
    print_exec("t1 = u1 / abs(u1); print(t1)")
    print_exec("t2 = u2 / abs(u2); print(t2)")

    pause()
    check_exec("set t3 to be the unit vector parallel to u3.",
               "t3",
               "u3 / abs(u3)")

    pause()
    print_line()
    print_fill("Now, check that the vectors t1, t2 and t3 are mutually orthogonal.")
    print_line()
    print_exec("print(t1 * t2 + t2 * t1)")
    print_exec("print(t1 * t3 + t3 * t1)")

    pause()
    check_eval("check that t2 is orthogonal to t3",
               "t2 * t3 + t3 * t2",
               "print({})")

    pause()
    print_line()
    print_fill("Note that the process used to obtain t1, t2 and t3 is essentially" +
              " Gram-Schmidt orthogonalization.")
    print_line()

    pause()
    print_line()
    print_fill("Example 12.1. the force F == 7{1}+4{2}+5{3} is applied to the point P == (4,6,7). " +
              " Determine t, the magnitude of the torque about the origin O.")
    print_line()

    pause()
    print_exec("R3 = istpq(3,0)")
    print_exec("F = clifford((7,4,5),R3); print(F)")
    print_exec("OP = clifford((4,6,7), R3); print(OP)")
    print_exec("t = abs(OP ^ F); print(t)")
    print_line()
    print_fill("Here we have used the alternate input notation for vectors, which uses a tuple or list of" +
              " coordinates, and a basis defined by an index_set, in this case R3 == {1,2,3}.")
    print_line()

    pause()
    print_line()
    print_fill("Example 12.2. Rotate the vector r == 4{1}+2{2}+2{3} about the axis a == 1.5{1}+2{2}" +
              " by the angle alpha == abs(a), to obtain the vector q.")
    print_line()

    pause()
    print_exec("a = 1.5*e(1)+2*e(2); print(a)")
    print_exec("j = e({1,2,3}); print(j)")
    print_exec("s = exp(-a * j / 2); print(s)")
    print_exec("r = e(1)+2*e(2)+2*e(3); print(r)")
    print_exec("q = r | s; print(q)")

    pause()
    print_line()
    print_fill("Note that we have used -a instead of a here, to allow the use of r | s. " +
              " Alternatively we could have used a and used the following method.")
    print_line()
    print_exec("s1 = exp(a * j / 2); print(s1)")
    print_exec("q1 = inv(s1) * r * s1; print(q1)")
    print_exec("print(abs(q-q1))")

    pause()
    print_line()
    print_fill("Because abs(a) == 2.5 < pi, we can represent the rotation uniquely by -a," +
              " except that 2*pi-a yields the same rotation because s and -s yield" +
              " the same rotation.")
    print_line()
    print_exec("s2 = exp((2 * pi - a) * j / 2); print(s2)")

    pause()
    check_exec("set q2 to be the vector r rotated by the rotor s2.",
               "q2",
               "r | s2")
    print_exec("print(abs(q-q2))")

    pause()
    print_line()
    print_fill("Example 12.3. Perform two successive rotations, the first one around the axis OA," +
              " A == (1,-1,1) by the angle 2*pi/3, and the second around the axis OB, B == (0,1,-1)" +
              " by the angle pi. " +
              " What is the axis of the combined rotation?")
    print_line()

    pause()
    print_exec("OA = e(1)-e(2)+e(3); print(OA)")
    print_exec("OB = e(2)-e(3); print(OB)")
    print_exec("a = OA/abs(OA) * 2*pi/3; print(a)")
    print_exec("b = OB/abs(OB) * pi; print(b)")
    print_exec("j = e({1,2,3}); print(j)")
    print_exec("u = exp(j * a / 2); print(u)")
    print_exec("v = exp(j * b / 2); print(v)")
    print_exec("z = log(u * v); print(z)")
    print_exec("c = j * z/abs(z); print(c)")
    print_line()
    print_fill("This means that the axis of the combined rotation is OC, C == (0,0,1).")
    print_line()

    pause()
    print_line()
    print_fill("You have completed the tutorial file pyclical_tutorial_1_2_space.py.")

if __name__ == "__main__":
    ctx = tutorial_context(globals())
    try:
        run(ctx)
    except:
        ctx.print_fill("The tutorial was interrupted.")
        pass
