# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_1_3_lorentz.py:
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
    ctx_methods = get_object_methods(ctx)
    for method in ctx_methods:
        exec(method + " = ctx_methods['" + method +"']")

    print_head("1.3 Electromagnetism and Lorentz transformations")
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
    print_fill("Example 13.1. Using the Clifford algebra of 3 dimensional Euclidean space," +
              " compute the Lorentz invariants, energy density and Poynting vector of the" +
              " electromagnetic field F == E-j*B, with an electric component E == {1}+2{2}+4{3}" +
              " and a magnetic component B == 3{1}+5{2}+7{3}, where j == {1,2,3}.")
    print_line()
    print_exec("E = cl('{1}+2{2}+4{3}'); print E")
    print_exec("B = cl('3{1}+5{2}+7{3}'); print B")
    print_exec("j = e({1,2,3})")
    print_exec("F = E - j*B; print F")
    print_line()
    print_fill("Lorentz invariants")
    print_line()
    print_exec("Lorentz = F * F / 2; print Lorentz")
    print_line()
    print_fill("Check that Lorentz(0) == (E**2 - B**2)/2; Lorentz(3)*j == E & B.")
    print_line()
    print_exec("print Lorentz(0) - (E**2 - B**2)/2")
    print_exec("print Lorentz(3)*j - (E & B)")
    print_line()
    print_fill("Energy density and Poynting vector.")
    print_line()
    print_exec("EP = -involute(F) * F / 2; print EP")
    print_line()
    print_fill("Check that EP(0) == (E**2 + B**2)/2; EP(1) == -(E ^ B)*j.")
    print_line()
    print_exec("print EP(0) - (E**2 + B**2)/2")
    print_exec("print EP(1) + (E ^ B)*j")

    pause()
    print_line()
    print_fill("Example 13.2. An observer O sees another observer O2 moving along the {1} axis" +
              " at half the speed of light. Using the Clifford algebra of 3 dimensional" +
              " Euclidean space, compute how the observer O2 decomposes the electromagnetic" +
              " field F == E-j*B, with an electric component E == {1}+2{2}+4{3} and a magnetic" +
              " component B == 3{1}+5{2}+7{3}, where j == {1,2,3}.")
    print_line()
    print_exec("v = 0.5*e(1); print v")
    print_exec("a = atanh(v); print a")
    print_exec("s = exp(a/2); print s")
    print_exec("u = inv(s); print u")
    print_line()
    print_fill("Use j and F from Example 13.1.")
    print_line()
    print_exec("F2 = s * F * u; print F2")
    print_exec("E2 = F2(1); print E2")
    print_exec("B2 = j * F2(2); print B2")

    pause()
    print_line()
    print_fill("Example 14.1. At 10 o'clock observer O observes the point (1,1,0)." +
              " Using the Clifford algebra of 3 dimensional Euclidean space," +
              " compute coordinates of this event in the frame of another observer O2" +
              " moving along the {1} axis at half the speed of light.")
    print_line()
    print_exec("x = 10+e(1)+e(2); print x")
    print_line()
    print_fill("Use s from Example 13.2.")
    print_line()
    print_exec("x2 = inv(involute(u)) * x * u; print x2")
    print_exec("print involute(x)*x")
    print_exec("print involute(x2)*x2")
    print_exec("print involute(x)*x - involute(x2)*x2")

    pause()
    print_line()
    print_fill("Example 14.2. Using the Clifford algebra defined by the index_set istpq(3,1)," +
              " compute the Lorentz invariants, energy density and Poynting vector of the" +
              " electromagnetic field F == E-j*B, with an electric component E == {-1,1}+2{-1,2}+4{-1,3}" +
              " and a magnetic component B == 3{-1,1}+5{-1,2}+7{-1,3}, where j == {1,2,3}.")
    print_line()
    print_exec("E = cl('{-1,1}+2{-1,2}+4{-1,3}'); print E")
    print_exec("B = cl('3{-1,1}+5{-1,2}+7{-1,3}'); print B")
    print_exec("j = e({-1,1,2,3})")
    print_exec("F = E + j * B; print F")
    print_line()
    print_fill("Lorentz invariants")
    print_line()
    print_exec("Lorentz = F * F / 2; print Lorentz")
    print_line()
    print_fill("Check that Lorentz(0) == (E**2 - B**2)/2; Lorentz(4)/j == E & B.")
    print_line()
    print_exec("print Lorentz(0) - (E**2 - B**2)/2")
    print_exec("print Lorentz(4)/j - (E & B)")
    print_line()
    print_fill("Energy density and Poynting vector.")
    print_line()
    print_exec("EP = -involute(F) * e(-1) * F / 2; print EP")
    print_line()
    print_fill("Check that EP[-1] == (E**2 + B**2)/2.")
    print_line()
    print_exec("print EP[-1] - (E**2 + B**2)/2")
    print_line()
    print_fill("Check that EP - EP[-1]*{-1} == -(({-1}*E) ^ ({-1}*B))*{1,2,3}.")
    print_line()
    print_exec("print EP - EP[-1]*e(-1) + ((e(-1)*E) ^ (e(-1))*B)*e({1,2,3})")

    pause()
    print_line()
    print_fill("You have completed the tutorial file pyclical_tutorial_1_3_lorentz.py.")

if __name__ == "__main__":
    ctx = tutorial_context(globals())
    try:
        run(ctx)
    except:
        ctx.print_fill("The tutorial was interrupted.")
        pass
