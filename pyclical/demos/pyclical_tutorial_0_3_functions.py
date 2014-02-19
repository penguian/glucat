# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_0_3_functions.py:
#
# This file contains a tutorial that explains the algebraic functions on
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

    print_head("0.3 Algebraic functions.")
    print_line()
    print_fill("This file contains a tutorial that explains the algebraic functions on" +
              " Clifford algebra elements available within PyClical.")
    print_line()
    print_fill("It is recommended that you do the tutorials in order, beginning with" +
              " 0.0 Notation.")
    print_line()

    pause()
    print_line()
    print_fill("PyClical is based on two Python classes, clifford, and index_set.")
    print_line()
    print_fill("The clifford class implements Clifford algebras over the double precision" +
              " floating point approximation to the real numbers.")
    print_line()
    print_fill("Tutorial 0.0 Notation showed you how to construct objects of type clifford. " +
              " This tutorial shows you how to use some functions defined on these objects.")

    pause()
    print_line()
    print_fill("First we will create some clifford objects to work on.")
    print_line()
    print_exec("w = 1+3*e(1)-2*e(3)+4*e({-2,3}); print w")
    print_exec("x = 2+3*e(1)-5*e(2)+5*e({-1,3}); print x")

    pause()
    print_line()
    print_fill("Examples: Arithmetic functions.")
    print_line()
    print_fill("These can be called either as a member function or as an ordinary function.")
    print_line()
    print_fill("inv(): Geometric multiplicative inverse.")
    print_line()
    print_exec("print x.inv()")
    print_exec("print inv(x)")

    print_line()
    print_fill("pow(m): Power: self to the m.")
    print_line()
    print_exec("print x.pow(0)")
    print_exec("print x.pow(1)")
    print_exec("print x.pow(3)")
    print_exec("print pow(x, 3)")

    print_line()
    print_fill("outer_pow(m): Outer product power.")
    print_line()
    print_fill("The outer product power is only defined for m as a non-negative integer.")
    print_line()
    print_exec("print x.outer_pow(0)")
    print_exec("print x.outer_pow(1)")
    print_exec("print x.outer_pow(3)")
    print_exec("print outer_pow(x, 3)")

    pause()
    print_line()
    print_fill("Try these now.")
    print_line()
    check_exec("set y to be the inverse of w.",
               "y",
               "inv(w)")
    check_exec("set z to be the cube of x.",
               "z",
               "pow(x, 3)")
    check_exec("set y to be w to the power of w.",
               "y",
               "pow(w, w)")
    check_exec("set z to be the outer cube of x.",
               "z",
               "outer_pow(x, 3)")

    pause()
    print_line()
    print_fill("Examples: Parts of clifford objects.")
    print_line()
    print_fill("Most of these can be called as member functions or as ordinary functions.")
    print_line()
    print_fill("One exception is for the grade-vector part of a clifford object. " +
              " The notation for this uses round brackets () immediately after the" +
              " clifford object, as if it was being called as a function. " +
              " The grade must be a non-negative integer. " +
              " This works with constants as well as variables.")
    print_line()
    print_fill("Examples:")
    print_line()
    print_exec("print x(0)")
    print_exec("print x(1)")
    print_exec("print x(2)")
    print_exec("print x(3)")

    check_exec("set y to be the grade 3 part of w.",
               "y",
               "w(3)")

    pause()
    print_line()
    print_fill("Other functions yielding parts of clifford objects are:")
    print_line()
    print_line()
    print_fill("scalar(): Scalar part of clifford object.  Result is a scalar.")
    print_line()
    print_exec("print x.scalar()")
    print_exec("print scalar(x)")

    print_line()
    print_fill("real(): Synonym for scalar().  This is an ordinary function only.")
    print_line()
    print_exec("print scalar(x)")

    print_line()
    print_fill("pure(): Pure part of a clifford object, that is a clifford object minus" +
              " its scalar part.")
    print_line()
    print_exec("print x.pure()")
    print_exec("print pure(x)")

    print_line()
    print_fill("even(): Even part of clifford object, sum of even grade terms.")
    print_line()
    print_exec("print x.even()")
    print_exec("print even(x)")

    print_line()
    print_fill("odd(): Odd part of clifford object, sum of odd grade terms.")
    print_line()
    print_exec("print x.odd()")
    print_exec("print odd(x)")

    print_line()
    print_fill("vector_part(frm): Vector part of clifford object, as a Python list.")
    print_line()
    print_fill("This is only available as a member function. " +
              " It takes an optional parameter frm," +
              " which should be an index set defining a subalgebra large enough to" +
              " contain the clifford object.")
    print_line()
    print_exec("print x.vector_part()")
    print_exec("print x.vector_part(istpq(4,1))")

    pause()
    print_line()
    print_fill("Try these now.")
    print_line()
    check_exec("set y to be the scalar part of w.",
               "y",
               "scalar(w)")
    check_exec("set z to be the pure part of w.",
               "z",
               "pure(w)")
    check_eval("compare with w to check that the scalar and pure parts add to w.",
               "scalar(w) + pure(w)",
               "print abs( w - ({}) )")
    check_exec("set y to be the even part of w.",
               "y",
               "w.even()")
    check_exec("set z to be the odd part of w.",
               "z",
               "w.odd()")
    check_eval("compare with w to check that the even and odd parts add to w.",
               "even(w) + odd(w)",
               "print abs( w - ({}) )")
    check_exec("set y to be the vector part of w with respect to the frame {-2,1,3}",
               "y",
               "w.vector_part(ist({-2,1,3}))")

    pause()
    print_line()
    print_fill("Examples: Involutions of clifford objects.")
    print_line()
    print_fill("Pyclical implements the three Clifford algebra involutions, both as" +
              " member functions, and as ordinary functions:")
    print_line()
    print_fill("involute(): Grade involution, each {i} is replaced by -{i} in each term," +
              " eg. e(1) -> -e(1).")
    print_line()
    print_exec("print e(1).involute()")
    print_exec("print e({1,2}).involute()")
    print_exec("print x.involute()")
    print_exec("print involute(x)")

    pause()
    print_line()
    print_fill("reverse(): Reversion anti-automorphism, eg. e(1)*e(2) -> e(2)*e(1).")
    print_line()
    print_exec("print e(1).reverse()")
    print_exec("print e({1,2}).reverse()")
    print_exec("print x.reverse()")
    print_exec("print reverse(x)")

    print_line()
    print_fill("conj(): Clifford conjugation anti-automorphism.")
    print_line()
    print_exec("print e(1).conj()")
    print_exec("print e({1,2}).conj()")
    print_exec("print x.conj()")
    print_exec("print conj(x)")

    pause()
    print_line()
    print_fill("Try these now.")
    print_line()
    check_exec("set y to be the grade involute of w.",
               "y",
               "involute(w)")
    check_exec("set z to be the reverse of w.",
               "z",
               "reverse(w)")
    check_exec("set y to be the Clifford conjugate of x.",
               "y",
               "conj(x)")

    pause()
    print_line()
    print_fill("Examples: Properties of clifford objects.")
    print_line()
    print_fill("PyClical implements a number of functions yielding properties of clifford objects. " +
              " Four of these give different ideas of the size of an object. " +
              " These are all available as both member functions and ordinary functions.")
    print_line()
    print_fill("norm(): Norm == sum of squares of coordinates.")
    print_line()
    print_exec("print x.norm()")
    print_exec("print norm(x)")

    print_line()
    print_fill("abs(): Absolute value: square root of norm.")
    print_line()
    print_exec("print x.abs()")
    print_exec("print abs(x)")

    pause()
    print_line()
    print_fill("max_abs(): Maximum of absolute values of components of the clifford object.")
    print_line()
    print_exec("print x.max_abs()")
    print_exec("print max_abs(x)")

    print_line()
    print_fill("quad(): Quadratic form == (reverse(x)*x)(0).")
    print_line()
    print_exec("print x.quad()")
    print_exec("print quad(x)")

    pause()
    print_line()
    print_fill("Try these now.")
    print_line()
    check_exec("set y to be the norm of w.",
               "y",
               "norm(w)")
    check_exec("set z to be the absolute value of w.",
               "z",
               "abs(w)")
    check_exec("set y to be the maximum absolute value of the coordinates of x.",
               "y",
               "max_abs(x)")
    check_exec("set z to be the quadratic form of w.",
               "z",
               "quad(w)")

    pause()
    print_line()
    print_fill("Two more properties of a clifford object are available as member functions:")
    print_line()
    print_fill("frame(): Index set defining a subalgebra that contains the clifford object.")
    print_line()
    print_exec("s = w.frame(); print s")
    print_exec("print type(s)")

    print_line()
    print_fill("isnan(): Check if the clifford object contains any IEEE NaN values.")
    print_line()
    print_exec("print x.isnan()")
    print_exec("print (0/clifford(0)).isnan()")

    pause()
    print_line()
    print_fill("Examples: Approximations to clifford objects.")
    print_line()
    print_fill("The following member function simplifies the value of a clifford object. " +
              " This is especially useful when printing.")
    print_line()
    print_fill("truncated(limit): Remove all terms with relative size smaller than limit.")
    print_line()

    print_exec("print clifford('1e8+{1}+1e-8{1,2}').truncated(1.0e-6)")
    print_exec("print clifford('1e4+{1}+1e-4{1,2}').truncated(1.0e-6)")

    pause()
    print_line()
    print_fill("You have completed the tutorial file pyclical_tutorial_0_3_functions.py.")

if __name__ == "__main__":
    ctx = tutorial_context(globals())
    try:
        run(ctx)
    except:
        ctx.print_fill("The tutorial was interrupted.")
        pass
