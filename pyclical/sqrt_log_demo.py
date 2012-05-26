# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# sqrt_log_demo.py: Demonstrate various sqrt and log calculations with PyClical.
#
#    copyright            : (C) 2008-2012 by Paul C. Leopardi
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
import sys

def pause():
    if __name__ != "__main__" or hasattr(sys, 'ps1') or hasattr(sys, 'ipcompleter'):
        if sys.stdin.isatty():
            raw_input("Press ENTER to continue")

def print_exec(str):
  print ">>>", str
  exec str in globals()

print "# sqrt_log_demo.py: Demonstrate various sqrt and log calculations with PyClical."
print ""
print "    BASIC OPERATIONS."
print ""
print_exec("e1 = e(1)")
print_exec("print e1")
print ""
print_exec("sqrt_e1 = sqrt(e1)")
print_exec("print sqrt_e1")
print_exec("print sqrt_e1*sqrt_e1")
print_exec("print abs(sqrt_e1*sqrt_e1 - e1)")
print ""
print_exec("log_e1 = log(e1)")
print_exec("print log_e1")
print_exec("print exp(log_e1)")
print_exec("print abs(exp(log_e1) - e1)")

pause()
print ""
print_exec("v = clifford('-2{1}+2{2}-3{3}')")
print_exec("print v")
print ""
print_exec("sqrt_v = sqrt(v)")
print_exec("print sqrt_v")
print_exec("print sqrt_v*sqrt_v")
print_exec("print abs(sqrt_v*sqrt_v - v)")
print ""
print_exec("log_v = log(v)")
print_exec("print log_v")
print_exec("print exp(log_v)")
print_exec("print abs(exp(log_v) - v)")

pause()
print ""
print_exec("x = clifford('-2{1}+2{2}-3{3}+4{-1,1,2}')")
print_exec("print x")
print ""
print_exec("sqrt_x = sqrt(x)")
print_exec("print sqrt_x")
print_exec("print sqrt_x*sqrt_x")
print_exec("print abs(sqrt_x*sqrt_x - x)")
print ""
print_exec("log_x = log(x)")
print_exec("print log_x")
print_exec("print exp(log_x)")
print_exec("print abs(exp(log_x) - x)")

pause()
print ""
print_exec("x = random_clifford(istpq(2,1))")
print_exec("print x")
print ""
print_exec("sqrt_x = sqrt(x)")
print_exec("print sqrt_x")
print_exec("print sqrt_x*sqrt_x")
print_exec("print abs(sqrt_x*sqrt_x - x)")
print ""
print_exec("log_x = log(x)")
print_exec("print log_x")
print_exec("print exp(log_x)")
print_exec("print abs(exp(log_x) - x)")

pause()
print ""
print "    SQUARE ROOTS AND LOGARITHMS OF ROTORS"
print "    in Conformal Geometric Algebra R_(4,1)."
print ""
print "    Reference:"
print "    [D+V]: L. Dorst and R. Valkenburg, 'Square root and logarithm of rotors',"
print "    L. Dorst and J. Lasenby (editors), Guide to geometric algebra in practice,"
print "    Springer, 2011, Chapter 5, pp. 81-104." 
print ""
print "    Define a bivector B."
print ""
print_exec("B = clifford('2{-1,1}+3{-1,2}+4{-1,3}+5{-1,4}+{1,2}+2{1,3}+3{1,4}-4{2,3}-5{2,4}+{3,4}')")
print_exec("print B")
print ""
print "    Exponentiate -B/2 to obtain the rotor R."
print ""
print_exec("R = exp(-B/2); print R")

pause()
print ""
print "    SQUARE ROOTS OF THE ROTOR R"
print ""
print "    Find the square root of R, and check it."
print ""
print_exec("sqrt_R = sqrt(R)")
print_exec("print sqrt_R")
print_exec("print sqrt_R*sqrt_R")
print_exec("print abs(sqrt_R*sqrt_R - R)")

pause()
print ""
print "    Now use [D+V (5.4)] to try to obtain a square root of R."
print ""
print_exec("dv_disc = (1+R(0))**2 - (R(4))**2")
print_exec("print dv_disc")
print_exec("dv_sqrt_R = (1+R)*(1+R(0)-R(4))/(2*dv_disc)*(1+R(0)+R(4)+sqrt(dv_disc))/sqrt(1+R(0)+sqrt(dv_disc))")
print_exec("print dv_sqrt_R")
print_exec("print dv_sqrt_R*dv_sqrt_R")
print_exec("print abs(dv_sqrt_R*dv_sqrt_R - R)")

pause()
print ""
print "    The [D+V (5.4)] square root of R is even."
print ""
print_exec("print abs(odd(dv_sqrt_R))")
print ""
print "     The PyClical square root of R is odd."
print ""
print_exec("print abs(even(sqrt_R))")

pause()
print ""
print "    LOGARITHMS OF THE ROTOR R"
print ""
print ""
print "    Find the logarithm of R, and check it."
print ""
print_exec("log_R = log(R)")
print_exec("print log_R")
print_exec("print exp(log_R)")
print_exec("print abs(exp(log_R) - R)")

pause()
print ""
print "    Now use [D+V, Section 5.3] to try to obtain a logarithm of R or -R."
print ""
print ""
print "    Obtain the bivector F via the exterior derivative of the action of the rotor R, as per [D+V (5.24)]."
print ""
print_exec("F = 2 * (R(4)-R(0)) * R(2); print F")
print ""
print "    Check the split of F into commuting 2-blades, as per [D+V (5.25)]."
print ""
print_exec("dv_norm = lambda F : sqrt(sqrt( (2*scalar(F**2) - F**2) * F**2 ))")
print_exec("F_m = F*(1 - ((dv_norm(F))**2)/(F**2))/2")
print_exec("F_p = F*(1 + ((dv_norm(F))**2)/(F**2))/2")
print ""
print "    Check that F_m and F_p are commuting 2-blades."
print ""
print_exec("print abs(F_m-F_m(2))")
print_exec("print F_m**2")
print_exec("print abs(F_p-F_p(2))")
print_exec("print F_p**2")
print_exec("print abs(F_p*F_m - F_m*F_p)")

pause()
print ""
print "    Set S_{m/p} = F_{m/p}, as per [D+V, remark before (5.26)]."
print ""
print_exec("S_m = F_m")
print_exec("S_p = F_p")
print ""
print "    Obtain C_{m/p} = cosh(B_{-/+}), as per [D+V (5.26)]."
print ""
print_exec("C_m = -scalar((R**2)(2)/S_p)")
print_exec("C_p = -scalar((R**2)(2)/S_m)")

pause()
print ""
print "    Define atanh2 as per [D+V (5.21)]."
print ""
print_exec("atanh2 = lambda s, c : math.asinh(sqrt(scalar(s**2)))/sqrt(scalar(s**2))*s if scalar(s**2) > 0 else s if scalar(s**2) == 0 else math.atan2(sqrt(scalar(-s**2)),c)/sqrt(scalar(-s**2))*s if -1 <= scalar(s**2) < 0 else float('nan')")
print ""
print "    Display the values of atanh2 to be used to obtain the logarithm."
print ""
print_exec("print atanh2(S_m,C_m)")
print_exec("print atanh2(S_p,C_p)")

pause()
print ""
print "    Obtain dv_log_R = Log(R) as per [D+V (5.27)]."
print ""
print_exec("dv_log_R = -(atanh2(S_m,C_m)+atanh2(S_p,C_p))/2")
print_exec("print dv_log_R")
print ""
print "     Check that exp(dv_log_R) == R or -R as per [D+V, p. 97]."
print ""
print_exec("print exp(dv_log_R)")
print_exec("print abs(exp(dv_log_R) - R)")
print_exec("print abs(exp(dv_log_R) + R)")

pause()
print ""
print "     You have completed the demonstration file sqrt_log_demo.py."
