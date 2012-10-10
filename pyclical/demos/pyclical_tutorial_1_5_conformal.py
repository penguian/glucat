# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_1_5_conformal.py:
#
# This file contains a tutorial that explains the algebraic functions on
# Clifford algebra elements available within PyClical.
#
#    copyright            : (C) 2012 by Paul C. Leopardi
#
# Licensed under CC BY-SA 3.0 http://creativecommons.org/licenses/by-sa/3.0/

from PyClical import *
from pyclical_tutorial_utils import *

def tutorial():
    ctx = interaction_context(globals())
    print_exec = ctx.print_exec
    check_exec = ctx.check_exec
    check_eval = ctx.check_eval

    print "# pyclical_tutorial_1_5_conformal.py:"
    print ""
    print_fill("1.5 Conformal Geometric Algebra.")
    print ""
    print_fill("This tutorial shows you how to use PyClical to perform the" +
              " operations of Conformal Geometric Algebra.")
    print ""
    print_fill("It is recommended that you do the tutorials in order, beginning with" +
              " 0.0 Notation.")
    print ""

    pause()
    print ""
    print_fill("PyClical includes objects and functions relevant to Conformal" +
              " Geometric Algebra based on Euclidean 3D space." +
              " The conventions used match those of the text:" +
              " C. Doran and A. Lasenby, Geometric algebra for physicists, Cambridge, 2003." +
              " This text is referred to as [DL] in this tutorial.")

    pause()
    print ""
    print_fill("First we create a vector in 3D space.")
    print ""
    print_exec("v = e(1)-2*e(2)+4*e(3); print v")
    print ""

    pause()
    print ""
    print_fill("The function cga3 maps Euclidean 3D vectors from the Clifford" +
              " algebra defined by the index set {1,2,3} into the corresponding" +
              " Conformal Geometric Algebra (CGA) defined by the index set {-1,1,2,3,4}."
              " The mapping is defined by [DL (10.50)].")
    print ""
    print_exec("x = cga3(v); print x")
    print ""

    pause()
    print ""
    print_fill("The CGA image x of the vector v is a null vector.")
    print ""
    print_exec("print x * x")
    print ""

    pause()
    print ""
    print_fill("The function cga3std maps null vectors in CGA into a standard form." +
              " The mapping is defined by [DL (10.52)].")
    print ""
    print_fill("Here we check that the image x of the vector v is already in standard form.")
    print ""
    print_exec("sx = cga3std(x); print sx")
    print_exec("print norm(sx * sx)")

    pause()
    print ""
    print_fill("The function agc3 maps a standard null vector in GGA back to Euclidean 3D space." +
              " The mapping is defined by [DL (10.50)].")

    print ""
    print_exec("v1 = agc3(x); print v1")
    print_exec("print norm(v - v1)")
    print ""

    pause()
    print ""
    print_fill("The CGA includes two distinguished elements, as defined in [DL]:")
    print ""
    print_fill("The element nbar3 corresponds to the origin of Euclidean 3D space:")
    print ""
    print_exec("print nbar3")
    print ""
    print_fill("The element ninf3 corresponds to the point at infinity for Euclidean 3D space:")
    print ""
    print_exec("print ninf3")
    print ""

    pause()
    print ""
    print_fill("Transformations in the spin group corresponding to index set {-1,2,3,4}" +
              " are orthogonal transformations and therefore preserve the quadratic form." +
              " Therefore null vectors are transformed to null vectors.")
    print ""
    
    pause()
    print ""
    print_fill("We first create a bivector as the exterior product of two vectors:")
    print ""
    print_exec("u = cl('{-1}+3{1}+4{2}-3{3}+{4}')")
    print_exec("w = cl('{-1}+2{1}-3{2}-3{3}+{4}')")
    print_exec("b = u ^ w; print b")
    print ""
    print_fill("We then exponentiate to obtain an element of the spin group:")
    print ""
    print_exec("g = exp(b); print g")
    print ""

    pause()
    print ""
    print_fill("We now check that the action of g on x to produce y yields a null vector:")
    print ""
    check_exec("set y to be the result of the action of g on x.",
               "y",
               "x | g")
    print ""
    print_exec("print norm(y * y)")
    print ""

    pause()
    print ""
    check_exec("set sy to be the CGA standard form of y.",
               "sy",
               "cga3std(y)")
    print ""
    print_exec("print norm(sy * sy)")
    print ""

    pause()
    print ""
    print_fill("Here we show that our particular g moves the origin.")
    print ""
    print_exec("print nbar3")
    print_exec("print cga3std(nbar3)")
    print_exec("print agc3( cga3std(nbar3) )")
    print_exec("print nbar3 | g")
    print_exec("print cga3std( nbar3 | g)")
    print_exec("print agc3( cga3std(nbar3 | g) )")
    print ""

    pause()
    print ""
    print_fill("Our particular g does not move the point at infinity.")
    print ""
    print_exec("print ninf3")
    print_exec("print ninf3 | g")
    print ""

    pause()
    print ""
    print_fill("Our g is not an orthogonal transformation in Euclidean 3D space.")
    print ""
    print_exec("print v * v")
    print ""

    pause()
    print ""
    check_exec("set asy to be the Euclidean 3D vector corresponding to sy.",
               "asy",
               "agc3(sy)")
    print ""
    print_exec("print asy * asy")
    print ""
    
    pause()
    print ""
    print_fill("You have completed the tutorial file pyclical_tutorial_1_5_conformal.py.")

if __name__ == "__main__":
    try:
        tutorial()
    except:
        print_fill("The tutorial was interrupted.")
        pass
