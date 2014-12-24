#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# sqrt_log_demo.py: Demonstrate various sqrt and log calculations with PyClical.
#
#    copyright            : (C) 2008-2014 by Paul C. Leopardi
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
    ctx_methods = get_object_methods(ctx)
    for method in ctx_methods:
        exec(method + " = ctx_methods['" + method +"']")

    print_fill("# sqrt_log_demo.py: Demonstrate various sqrt and log calculations with PyClical.")
    print_line()
    print_fill("INITIALIZATION.")
    print_exec("from PyClical import *")
    print_line()
    print_fill("BASIC OPERATIONS.")
    print_line()
    print_exec("e1 = e(1)")
    print_exec("print e1")
    print_line()
    print_exec("sqrt_e1 = sqrt(e1)")
    print_exec("print sqrt_e1")
    print_exec("print sqrt_e1*sqrt_e1")
    print_exec("print abs(sqrt_e1*sqrt_e1 - e1)")
    print_line()
    print_exec("log_e1 = log(e1)")
    print_exec("print log_e1")
    print_exec("print exp(log_e1)")
    print_exec("print abs(exp(log_e1) - e1)")

    pause()
    print_line()
    print_exec("v = clifford('-2{1}+2{2}-3{3}')")
    print_exec("print v")
    print_line()
    print_exec("sqrt_v = sqrt(v)")
    print_exec("print sqrt_v")
    print_exec("print sqrt_v*sqrt_v")
    print_exec("print abs(sqrt_v*sqrt_v - v)")
    print_line()
    print_exec("log_v = log(v)")
    print_exec("print log_v")
    print_exec("print exp(log_v)")
    print_exec("print abs(exp(log_v) - v)")

    pause()
    print_line()
    print_exec("x = clifford('-2{1}+2{2}-3{3}+4{-1,1,2}')")
    print_exec("print x")
    print_line()
    print_exec("sqrt_x = sqrt(x)")
    print_exec("print sqrt_x")
    print_exec("print sqrt_x*sqrt_x")
    print_exec("print abs(sqrt_x*sqrt_x - x)")
    print_line()
    print_exec("log_x = log(x)")
    print_exec("print log_x")
    print_exec("print exp(log_x)")
    print_exec("print abs(exp(log_x) - x)")

    pause()
    print_line()
    print_exec("x = random_clifford(istpq(2,1))")
    print_line()
    print_exec("sqrt_x = sqrt(x)")
    print_exec("print abs(sqrt_x*sqrt_x - x)")
    print_line()
    print_exec("log_x = log(x)")
    print_exec("print abs(exp(log_x) - x)")

    pause()
    print_line()
    print_fill("SQUARE ROOTS AND LOGARITHMS OF ROTORS")
    print_fill("in Conformal Geometric Algebra R_(4,1).")
    print_line()
    print_fill("Reference:")
    print_fill("[D+V]: L. Dorst and R. Valkenburg, 'Square root and logarithm of rotors',")
    print_fill("L. Dorst and J. Lasenby (editors), Guide to geometric algebra in practice,")
    print_fill("Springer, 2011, Chapter 5, pp. 81-104.")
    print_line()
    print_fill("Define a bivector B.")
    print_line()
    print_exec("B = clifford('2{-1,1}+3{-1,2}+4{-1,3}+5{-1,4}+{1,2}+2{1,3}+3{1,4}-4{2,3}-5{2,4}+{3,4}')")
    print_exec("print B")
    print_line()
    print_fill("Exponentiate -B/2 to obtain the rotor R.")
    print_line()
    print_exec("R = exp(-B/2); print R")

    pause()
    print_line()
    print_fill("SQUARE ROOTS OF THE ROTOR R")
    print_line()
    print_fill("Find the square root of R, and check it.")
    print_line()
    print_exec("sqrt_R = sqrt(R)")
    print_exec("print sqrt_R")
    print_exec("print sqrt_R*sqrt_R")
    print_exec("print abs(sqrt_R*sqrt_R - R)")

    pause()
    print_line()
    print_fill("Now use [D+V (5.4)] to try to obtain a square root of R.")
    print_line()
    print_exec("dv_disc = (1+R(0))**2 - (R(4))**2")
    print_exec("print dv_disc")
    print_exec("dv_sqrt_R = (1+R)*(1+R(0)-R(4))/(2*dv_disc)*(1+R(0)+R(4)+sqrt(dv_disc))/sqrt(1+R(0)+sqrt(dv_disc))")
    print_exec("print dv_sqrt_R")
    print_exec("print dv_sqrt_R*dv_sqrt_R")
    print_exec("print abs(dv_sqrt_R*dv_sqrt_R - R)")

    pause()
    print_line()
    print_fill("The [D+V (5.4)] square root of R is even.")
    print_line()
    print_exec("print abs(odd(dv_sqrt_R))")
    print_line()
    print_fill(" The PyClical square root of R is odd.")
    print_line()
    print_exec("print abs(even(sqrt_R))")

    pause()
    print_line()
    print_fill("LOGARITHMS OF THE ROTOR R")
    print_line()
    print_line()
    print_fill("Find the logarithm of R, and check it.")
    print_line()
    print_exec("log_R = log(R)")
    print_exec("print log_R")
    print_exec("print exp(log_R)")
    print_exec("print abs(exp(log_R) - R)")

    pause()
    print_line()
    print_fill("Now use [D+V, Section 5.3] to try to obtain a logarithm of R or -R.")
    print_line()
    print_line()
    print_fill("Obtain the bivector F via the exterior derivative of the action of the rotor R, as per [D+V (5.24)].")
    print_line()
    print_exec("F = 2 * (R(4)-R(0)) * R(2); print F")
    print_line()
    print_fill("Check the split of F into commuting 2-blades, as per [D+V (5.25)].")
    print_line()
    print_exec("dv_norm = lambda F : sqrt(sqrt( (2*scalar(F**2) - F**2) * F**2 ))")
    print_exec("F_m = F*(1 - ((dv_norm(F))**2)/(F**2))/2")
    print_exec("F_p = F*(1 + ((dv_norm(F))**2)/(F**2))/2")
    print_line()
    print_fill("Check that F_m and F_p are commuting 2-blades.")
    print_line()
    print_exec("print abs(F_m-F_m(2))")
    print_exec("print F_m**2")
    print_exec("print abs(F_p-F_p(2))")
    print_exec("print F_p**2")
    print_exec("print abs(F_p*F_m - F_m*F_p)")

    pause()
    print_line()
    print_fill("Set S_{m/p} = F_{m/p}, as per [D+V, remark before (5.26)].")
    print_line()
    print_exec("S_m = F_m")
    print_exec("S_p = F_p")
    print_line()
    print_fill("Obtain C_{m/p} = cosh(B_{-/+}), as per [D+V (5.26)].")
    print_line()
    print_exec("C_m = -scalar((R**2)(2)/S_p)")
    print_exec("C_p = -scalar((R**2)(2)/S_m)")

    pause()
    print_line()
    print_fill("Define atanh2 as per [D+V (5.21)].")
    print_line()
    print_exec("atanh2 = lambda s, c :"
               + " math.asinh(sqrt(scalar(s**2)))/sqrt(scalar(s**2))*s if scalar(s**2) > 0"
               + " else s if scalar(s**2) == 0"
               + " else math.atan2(sqrt(scalar(-s**2)),c)/sqrt(scalar(-s**2))*s if -1 <= scalar(s**2) < 0"
               + " else float('nan')")
    print_line()
    print_fill("Display the values of atanh2 to be used to obtain the logarithm.")
    print_line()
    print_exec("print atanh2(S_m,C_m)")
    print_exec("print atanh2(S_p,C_p)")

    pause()
    print_line()
    print_fill("Obtain dv_log_R = Log(R) as per [D+V (5.27)].")
    print_line()
    print_exec("dv_log_R = -(atanh2(S_m,C_m)+atanh2(S_p,C_p))/2")
    print_exec("print dv_log_R")
    print_line()
    print_fill("Check that exp(dv_log_R) == R or -R as per [D+V, p. 97].")
    print_line()
    print_exec("print exp(dv_log_R)")
    print_exec("print abs(exp(dv_log_R) - R)")
    print_exec("print abs(exp(dv_log_R) + R)")
    print_line()
    print_fill("Use the pseudoscalar i = e({-1,1,2,3,4}) to obtain the log of -R or R.")
    print_line()
    print_exec("i = e({-1,1,2,3,4})")
    print_line()
    print_fill("First, verify that i has the required properties: squares to -1 and commutes with vectors.")
    print_line()
    print_exec("print i*i")
    print_exec("print i*e(-1) - e(-1)*i")
    for k in range(1,5):
        print_exec("print i*e("+str(k)+") - e("+str(k)+")*i")
    print_line()
    print_fill("Now use the properties of i to find an expression for the log of -R or R.")
    print_line()
    print_exec("print abs(exp(dv_log_R + i*pi) - R)")
    print_exec("print abs(exp(dv_log_R + i*pi) + R)")

    pause()
    print_line()
    print_fill("You have completed the demonstration file sqrt_log_demo.py.")

if __name__ == "__main__":
    try:
        run()
    except:
        print("The demo was interrupted.")
        pass
