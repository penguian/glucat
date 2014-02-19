# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_0_2_operations.py:
#
# This file contains a tutorial that explains the unary and binary operations on
# Clifford algebra elements available within PyClical.
#
#    copyright            : (C) 2012-2014 by Paul C. Leopardi
#
# Licensed under CC BY-SA 3.0 http://creativecommons.org/licenses/by-sa/3.0/

from pyclical_tutorial_utils import *

def run(ctx):
    ctx_methods = get_object_methods(ctx)
    for method in ctx_methods:
        exec(method + " = ctx_methods['" + method +"']")

    print_head("0.2 Operations.")
    print_line()
    print_fill("This file contains a tutorial that explains the unary and binary operations" +
              " on Clifford algebra elements available within PyClical.")
    print_line()
    print_fill("It is recommended that you do the tutorials in order, beginning with" +
              " 0.0 Notation.")
    print_line()
    print_exec("from PyClical import *")

    pause()
    print_line()
    print_fill("PyClical is based on two Python classes, clifford, and index_set.")
    print_line()
    print_fill("The clifford class implements Clifford algebras over the double precision" +
              " floating point approximation to the real numbers.")
    print_line()
    print_fill("Tutorial 0.0 Notation showed you how to construct objects of type clifford. " +
              " This tutorial shows you how to use unary and binary operations to" +
              " manipulate these objects.")

    pause()
    print_line()
    print_fill("First we will create some clifford objects to work on.")
    print_line()
    print_exec("w = 1+3*e(1)+5*e({1,3}); print w")
    print_exec("x = 2-5*e(2)+5*e({-1,3}); print x")

    pause()
    print_line()
    print_fill("Examples: Unary operations")
    print_line()
    print_fill("The unary operators + and - are defined for clifford objects.")

    pause()
    print_fill("Try these now.")
    print_line()
    check_exec("set y to be the unary operator + applied to w.",
               "y",
               "+w")
    check_exec("set z to be the unary operator - applied to x, that is, the negative of x.",
               "z",
               "-x")

    pause()
    print_line()
    print_fill("Examples: Comparison operations")
    print_line()
    print_fill("The binary comparison operators == and != are defined for clifford objects," +
              " and are also defined between clifford objects and floating point numbers.")
    print_line()

    pause()
    print_fill("Try these now.")
    print_line()
    check_eval("check that x * 1.0 equals x.",
               "(x * 1.0) == x",
               "print {}")
    check_eval("check that w does not equal x.",
               "w != x",
               "print {}")
    check_eval("check that 2 does not equal x.",
               "2 != x",
               "print {}")

    pause()
    print_line()
    print_fill("Examples: Binary addition and subtraction")
    print_line()
    print_fill("The binary operators + and - are defined for clifford objects and are also" +
              " defined between clifford objects and floating point numbers.")
    print_line()

    pause()
    print_fill("Try these now.")
    print_line()
    check_exec("set y to be the sum of w and x.",
               "y",
               "w + x")
    check_exec("set z to be the difference between w and x, that is, w minus x.",
               "z",
               "w - x")
    check_exec("set y to be the sum of w and 2.5.",
               "y",
               "w + 2.5")
    check_exec("set z to be the difference between 3.5 and x, that is, 3.5 minus x.",
               "z",
               "3.5 - x")

    pause()
    print_line()
    print_fill("Examples: Assignment addition and subtraction")
    print_line()
    print_fill("Objects of type clifford can be modified by using assignment addition and" +
              " subtraction operators. " +
              " These can be used to add or subtract other clifford" +
              " objects or floating point numbers.")
    print_line()
    print_fill("Examples, with w, x, y and z set to known values:")
    print_line()
    print_exec("w = 1+3*e(1)+5*e({1,3}); print w")
    print_exec("x = 2-5*e(2)+5*e({-1,3}); print x")
    print_exec("y = 1+3*e(1)-5*e({1,3}); print y")
    print_exec("z = 2-5*e(2)-5*e({-1,3}); print z")
    print_line()
    print_fill("Using assignment addition to set y to be the sum of y and x:")
    print_line()
    print_exec("y += x; print y")
    print_line()
    print_fill("Using assignment subtraction to set z to be the difference between z and x:")
    print_line()
    print_exec("z -= x; print z")
    print_line()
    print_fill("Using assignment addition to add 2.5 to y:")
    print_line()
    print_exec("y += 2.5; print y")
    print_line()
    print_fill("Using assignment subtraction to subtract 3.5 from z:")
    print_line()
    print_exec("z -= 3.5; print z")

    pause()
    print_line()
    print_fill("Examples: Multiplication operators")
    print_line()
    print_fill("Unlike Python floating point numbers, there are 4 different multiplication" +
              " operators defined for clifford objects, which also work between clifford" +
              " objects and floating point numbers.")
    print_line()
    print_fill("These are: * Clifford or geometric product; also scalar multiplication.")
    print_fill("           ^ outer (exterior, wedge, Grassmann) product.")
    print_fill("           & Hestenes inner product.")
    print_fill("           % Left contraction.")
    print_line()

    pause()
    print_fill("Try these now.")
    print_line()
    check_exec("set y to be the Clifford product of w and x.",
               "y",
               "w * x")
    check_exec("set z to be the wedge product of w and x.",
               "z",
               "w ^ x")
    check_exec("set y to be the Hestenes inner product of w and x.",
               "y",
               "w & x")
    check_exec("set z to be the left contraction between w and x.",
               "z",
               "w % x")

    pause()
    print_line()
    print_fill("Note that the precedence for these operators is the same as in Python and" +
              " may not agree with convention. " +
              " In particular, the precedence of & and ^ is *lower* than that of + and -.")
    print_line()
    print_exec("print x + y ^ z")
    print_exec("print (x + y) ^ z")
    print_exec("print x + (y ^ z)")
    print_exec("print x + y & z")
    print_exec("print (x + y) & z")
    print_exec("print x + (y & z)")

    pause()
    print_line()
    print_fill("Each of these operators has corresponding assignment versions.")
    print_line()
    print_fill("Examples, with w, x, y and z set to known values:")
    print_line()
    print_exec("w = 1+3*e(1)+5*e({1,3}); print w")
    print_exec("x = 2-5*e(2)+5*e({-1,3}); print x")
    print_exec("y = 1+3*e(1)-5*e({1,3}); print y")
    print_exec("z = 2-5*e(2)-5*e({-1,3}); print z")
    print_exec("y *= w; print y")
    print_exec("y *= 2.0; print y")
    print_exec("z ^= w; print z")
    print_exec("z &= x; print z")
    print_exec("y %= x; print y")

    pause()
    print_line()
    print_fill("Examples: Division and transformation operators.")
    print_line()
    print_fill("Recall that the Clifford product * is not commutative. " +
              " The binary operation x / y, where y is a clifford object, means x * inv(y)," +
              " where y*inv(y) == inv(y)*y == 1, if such an inverse exists.")
    print_line()
    print_fill("Try this now.")
    print_line()
    check_exec("set y to be x divided by w.",
               "y",
               "x / w")

    pause()
    print_line()
    print_fill("The binary operation x | y, where x and y are clifford object," +
              " means  y * x / involute(y), where involute(y) is the grade involution" +
              " applied to y. " +
              " This operation yields the 'twisted adjoint' action of y on x.")
    print_line()

    pause()
    print_fill("Try this now.")
    print_line()
    check_exec("set y to be x transformed by the twisted adjoint action of w.",
               "y",
               "x | w")

    pause()
    print_line()
    print_fill("These operators also have corresponding assignment versions.")
    print_line()
    print_fill("Examples, with w, x, y and z set to known values:")
    print_line()
    print_exec("w = 1+e({1,2})+5*e({1,3}); print w")
    print_exec("x = 2-5*e(2)+5*e({-1,3}); print x")
    print_exec("y = 1+3*e(1)-5*e({1,3}); print y")
    print_exec("z = 2*e(1)-5*e(2)-5*e(3); print z")
    print_exec("y /= w; print y")
    print_exec("y /= 2.0; print y")
    print_exec("z |= w; print z")

    pause()
    print_line()
    print_fill("Examples: Exponentiation operators.")
    print_line()
    print_fill("Exponentiation is expressed as usual in Python using the operation x ** y. " +
              " Here either the x or y can be a clifford object.")
    print_line()

    pause()
    print_fill("Try this now.")
    print_line()

    pause()
    check_exec("set y to be the cube of x.",
               "y",
               "x ** 3")
    check_exec("set z to be 3 to the power of x.",
               "z",
               "3 ** x")
    check_exec("set y to be w to the power of x.",
               "y",
               "w ** x")

    pause()
    print_line()
    print_fill("You have completed the tutorial file pyclical_tutorial_0_2_operations.py.")

if __name__ == "__main__":
    ctx = tutorial_context(globals())
    try:
        run(ctx)
    except:
        ctx.print_fill("The tutorial was interrupted.")
        pass
