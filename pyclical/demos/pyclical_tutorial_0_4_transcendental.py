# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_0_4_transcendental.py:
#
# This file contains a tutorial that explains the algebraic functions on
# Clifford algebra elements available within PyClical.
#
#    copyright            : (C) 2012-2014 by Paul C. Leopardi
#
# Licensed under CC BY-SA 3.0 http://creativecommons.org/licenses/by-sa/3.0/

from pyclical_tutorial_utils import *

def run(ctx):
    for name, method in get_object_methods(ctx).items():
        exec("global "+name+";"+name+"=method")

    print_head("0.4 Square root and transcendental functions.")
    print_line()
    print_fill("This file contains a tutorial that describes the square root and " +
              " transcendental functions on Clifford algebra elements available within PyClical.")
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
              " This tutorial shows you how to use the square root and " +
              " transcendental functions defined on these objects.")

    pause()
    print_line()
    print_fill("First we will create a clifford object to work on.")
    print_line()
    print_exec("w = 1+3*e(1)-2*e(3)+4*e({-2,3}); print(w)")
    print_exec("x = 2+3*e(1)-5*e(2)+5*e({-1,3}); print(x)")

    pause()
    print_line()
    print_fill("The complexifier function.")
    print_line()
    print_fill("One key property of Clifford algebras is that they generally contain many" +
              " square roots of -1. " +
              " Some Clifford algebras have a pseudoscalar e(s) that squares to -1 and" +
              " commutes with all algebra elements. " +
              " Such a Clifford algebra is isomorphic to a full matrix algebra over" +
              " the complex field, and in this case the pseudoscalar e(s) corresponds to" +
              " the complex i times the identity matrix.")
    print_line()
    print_fill("The complexifier function takes a clifford object x and returns the pseudoscalar e(s)" +
              " corresponding to the smallest Clifford algebra A containing x," +
              " such that A is isomorphic to a full complex matrix algebra. " +
              " In some cases, there are two such algebras and an arbitrary choice is made.")
    print_line()
    print_exec("print(complexifier(1))")
    print_exec("print(complexifier(e(1)))")
    print_exec("print(complexifier( e({1,2}) ))")
    print_exec("j = complexifier(w); print(j)")
    print_exec("print(j * j)")
    print_exec("print(j * w - w * j)")

    pause()
    check_exec("set i to be the complexifier of x",
               "i",
               "complexifier(x)")
    check_eval("compare with -1 to check that i is a square root of -1",
               "i * i",
               "print(1 + {})")
    check_eval("compare with x * i to check that i commutes with x",
               "i * x",
               "print(x * i - {})")

    pause()
    print_fill("The complexifier function can also take an index set or a string as an argument.")
    print_line()
    print_fill("Examples:")
    print_line()
    print_exec("print(complexifier( ist(1) ))")
    print_exec("print(complexifier( ist({1,2}) ))")
    print_exec("print(complexifier('{1,2}'))")

    pause()
    print_fill("The square root function and many of the transcendental functions" +
              " take a complexifier i as an optional argument, i.e. f(x, i). " +
              " The default is to use complexifier(x) as the complexifier.")
    print_line()
    print_fill("In the case of the square root and logarithm functions, this is because" +
              " the function needs to use a complexifier i in the case where" +
              " the matrix representation of the clifford object x has a negative real eigenvalue. "
              " See [L-2010] for more details.")
    print_line()
    print_fill("The other functions that use a complexifier i," +
              " do so either because they are implemented by an expression using i," +
              " or because the function eventually calls sqrt or log.")

    pause()
    print_line()
    print_fill("sqrt(x, i): The square root function.")
    print_line()
    print_fill("Examples:")

    print_exec("print(sqrt(-1))")
    print_exec("j = sqrt(-1, complexifier(e(1))); print(j); print(j*j)")
    print_exec("j = sqrt(-1,'{1,2,3}'); print(j); print(j*j)")
    print_exec("print(sqrt(2*e(-1)))")

    pause()
    print_fill("For the examples for the transcendental functions," +
              " we will use the following clifford objects.")
    print_exec("i = e(-1)")
    print_exec("j = e({1,2})")
    print_exec("k = e({1,2,3})")

    pause()
    print_fill("exp(x): Exponential function.")
    print_fill("log(x, i): The natural logarithm function.")
    print_line()
    print_fill("Examples:")

    print_exec("print(exp(j * pi/2))")
    print_exec("print((log(j) / pi))")
    print_exec("print((log(j,'{1,2,3}') / pi))")

    pause()
    print_fill("cos(x, i): The cosine function.")
    print_fill("acos(x, i): The inverse cosine function.")
    print_line()
    print_fill("Examples:")

    print_exec("print(cos(acos(j)))")
    print_exec("print(cos(acos(j), k))")

    pause()
    print_fill("cosh(x): The hyperbolic cosine function.")
    print_fill("acosh(x, i): The inverse hyperbolic cosine function.")
    print_line()
    print_fill("Examples:")

    print_exec("print(cosh(j * pi))")
    print_exec("print(acosh(0) / pi)")
    print_exec("print(acosh(0,'{-2,-1,1}') / pi)")
    print_exec("print(cosh(acosh(j)))")
    print_exec("print(cosh(acosh(k)))")

    pause()
    print_fill("sin(x, i): The sine function.")
    print_fill("asin(x, i): The inverse sine function.")
    print_line()
    print_fill("Examples:")

    print_exec("print(asin(sin(i, i), i))")
    print_exec("print(asin(sin(i, i),'{-2,-1,1}'))")
    print_exec("print(asin(sin(k)))")
    print_exec("print(asin(1) / pi)")

    pause()
    print_fill("sinh(x): The hyperbolic sine function.")
    print_fill("asinh(x, i): The inverse hyperbolic sine function.")
    print_line()
    print_fill("Examples:")

    print_exec("print(sinh(j * pi/2))")
    print_exec("print(asinh(j) / pi)")
    print_exec("print(asinh(j, k) / pi)")

    pause()
    print_fill("tan(x, i): The tangent function.")
    print_fill("atan(x, i): The inverse tangent function.")
    print_line()
    print_fill("Examples:")

    print_exec("print(tan(j))")
    print_exec("print(tan(j) * (exp(1) + exp(-1)) / (exp(1) - exp(-1)))")
    print_exec("print(tan(j, k))")
    print_exec("print(tan(atan(e(1))))")
    print_exec("print(tan(atan(e(1), k), k))")

    pause()
    print_fill("tanh(x): The hyperbolic tangent function.")
    print_fill("atanh(x, i): The inverse hyperbolic tangent function.")
    print_line()
    print_fill("Examples:")

    print_exec("print(tanh(j * pi/4))")
    print_exec("print(tanh(atanh(j)))")

    pause()
    print_line()
    print_fill("Try these now.")

    check_exec("set i to be the complexifier of w",
               "i",
               "complexifier(w)")
    check_exec("set y to be the exponential of (i times w)",
               "y",
               "exp(i * w)")
    check_exec("set z to be the cosine of w plus i times the sine of w",
               "z",
               "cos(w) + i * sin(w)")
    print_fill("Checking:")
    print_exec("print(y - z)")

    pause()
    check_exec("set y to be the exponential of x",
               "y",
               "exp(x)")
    check_exec("set z to be the hyperbolic cosine of x plus the hyperbolic sine of x",
               "z",
               "cosh(x) + sinh(x)")
    print_fill("Checking:")
    print_exec("print(y - z)")

    pause()
    print_line()
    print_fill("You have completed the tutorial file pyclical_tutorial_0_4_transcendental.py.")

if __name__ == "__main__":
    ctx = tutorial_context(globals())
    try:
        run(ctx)
    except:
        ctx.print_fill("The tutorial was interrupted.")
        pass
