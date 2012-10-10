# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_demo.py: This file is an almost complete port of the DEMO file from
#                   CLICAL by Pertti Lounesto, R. Mikkola, V. Vierros, 1987-1994.
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
from pyclical_tutorial_utils import *

def demo():
    ctx = interaction_context(globals())
    print_exec = ctx.print_exec

    print "# pyclical_demo.py: This file is an almost complete port of the DEMO file from"
    print "#                   CLICAL by Pertti Lounesto, R. Mikkola, V. Vierros, 1987-1994."
    print ""
    print_fill("This demonstration file contains examples on the following topics:")
    print_fill("  -Rotations in three dimensions R^3")
    print_fill("  -Lorentz transformations of the Minkowski space-time R^(3,1) ")
    print_fill("  -Electromagnetism in the Clifford algebras R_3 and R_(3,1)")
    print_fill("  -Complex Clifford algebra CxR_(3,1)")
    print_fill("  -Cayley algebra of octonions")
    print_fill("  -Conformal transformations and Vahlen matrices")
    print_fill("  -Elementary functions in the exterior algebra")
    print_fill("  -Pure spinors")
    print_fill("  -Conformal covariance of the Maxwell and Dirac equations")

    pause()
    print ""
    print_fill("ROTATIONS IN THREE DIMENSIONS R^3")
    print ""
    print_exec("x = 2*e(1)+3*e(2)+5*e(3); print x")
    print_exec("a = (3*e(1)+2*e(2)+e(3))/10; print a")
    print_exec("s = exp(e({1,2,3})*a/2); print s")

    pause()
    print ""
    print_fill("Rotate the vector x about the axis a by the angle  |a|:")
    print ""
    print_exec("y = s*x/s; print y")
    print ""
    print_fill("Check that the length is preserved by computing the squares of the vectors.")
    print ""
    print_exec("print y*y")
    print_exec("print x*x")

    pause()
    print ""
    print_fill("LORENTZ TRANSFORMATIONS IN THE MINKOWSKI SPACE-TIME R^(3,1)")
    print ""
    print_exec("A = -3*e({-1,1})-e({-1,2})+4*e({1,2})-2*e({-1,3})+3*e({1,3})+e({2,3}); print A")
    print_exec("s = exp(A/2); print s")
    print ""
    print_fill("Lorentz transformation of a space-time point (or event).")
    print ""
    print_exec("x = 2*e(-1)+2*e(1)+3*e(2)+5*e(3); print x")
    print_exec("y = s*x/involute(s); print y")

    pause()
    print ""
    print_fill("Check that the square-norm is preserved.")
    print ""
    print_exec("print y*y")
    print_exec("print x*x")

    pause()
    print ""
    print_fill("Lorentz transformation of an electromagnetic bivector F.")
    print ""
    print_exec("F = -2*e({-1,1})-e({-1,2})+5*e({1,2})-2*e({-1,3})+e({1,3})+7*e({2,3}); print F")
    print_exec("G = s*F/s; print G")
    print ""
    print_fill("Check that the Lorentz invariants are preserved (j == -{-1,1,2,3}).")
    print ""
    print_exec("print G*G/2")
    print_exec("print F*F/2")

    pause()
    print ""
    print_fill("Compute the Poynting vector and the energy density.")
    print ""
    print_exec("print -F*e(-1)*F/2")

    pause()
    print ""
    print_fill("ELECTROMAGNETISM IN THE CLIFFORD ALGEBRAS R_3 AND R_(3,1)")
    print ""
    print_fill("First, begin in the three-dimensional Euclidean space R^3.")
    print ""
    print_exec("E = e(1)+2*e(2)+4*e(3); print E")
    print_exec("B = 3*e(1)+5*e(2)+7*e(3); print B")
    print ""
    print_fill("The role of the imaginary unit is played by the unit volume element.")
    print ""
    print_exec("i = e({1,2,3}); print i")
    print_exec("F = E-i*B; print F")

    pause()
    print ""
    print_fill("Compute the Lorentz invariants.")
    print ""
    print_exec("print (E*E-B*B)/2")
    print_exec("print E&B")
    print ""
    print_fill("The above computations can be combined in one formula.")
    print ""
    print_exec("print F*F/2")

    pause()
    print ""
    print_fill("Compute the energy density and the Poynting vector.")
    print ""
    print_exec("print (E*E+B*B)/2")
    print_exec("print (E^B)/i")
    print ""
    print_fill("The above computations can be combined in one formula.")
    print ""
    print_exec("print -involute(F)*F/2")

    pause()
    print ""
    print_fill("Consider a boost at half the velocity of light" +
              " in the direction of the positive x-axis.")
    print ""
    print_exec("v = 0.5*e(1); print v")
    print_exec("a = atanh(v); print a")
    print_exec("s = exp(a/2); print s")

    pause()
    print ""
    print_fill("Lorentz transformation of the electromagnetic field F.")
    print ""
    print_exec("G = s*F/s; print G")

    pause()
    print ""
    print_fill("Extract out of G its magnetic induction.")
    print ""
    print_exec("print i*G(2)")

    pause()
    print ""
    print_fill("Lorentz transformation of a space-time point.")
    print ""
    print_exec("x = 10+e(1)+e(2); print x")
    print_exec("y = s*x/involute(s); print y")

    pause()
    print ""
    print_fill("Check that the quadratic form is preserved.")
    print ""
    print_exec("print y*conj(y)")
    print_exec("print x*conj(x)")

    pause()
    print ""
    print_fill("Next, perform the same computations in the space-time R^(3,1).")
    print ""

    pause()
    print_exec("E = E*e(-1); print E")
    print_exec("B = B*e(-1); print B")
    print ""
    print_fill("The role of the imaginary unit is played by the unit volume element.")
    print ""
    print_exec("i = -e({-1,1,2,3}); print i")
    print_exec("F = E-i*B; print F")

    pause()
    print ""
    print_fill("Compute the Lorentz invariants.")
    print ""
    print_exec("print F*F/2")

    pause()
    print ""
    print_fill("Compute the Poynting vector and the energy density.")
    print ""
    print_exec("print -F*e(-1)*F/2")

    pause()
    print ""
    print_fill("COMPLEX CLIFFORD ALGEBRA CxR_(3,1)")
    print ""
    print_fill("The role of the imaginary unit is played by a unit bivector of R_(3,3).")
    print ""
    print_exec("i = e({-3,-2}); print i")

    pause()
    print_exec("x = e(-1)+3*e(1)+2*e(2); print x")
    print_exec("y = 5*e(-1)+e(3); print y")
    print_exec("z = x+i*y; print z")

    pause()
    print_exec("A = -2*e({-1,2})-3*e({-1,3})+2*e({1,2}); print A")
    print_exec("B = -3*e({-1,1})-2*e({1,3})+e({2,3}); print B")
    print_exec("C = A+i*B; print C")

    pause()
    print_exec("s = exp(C/2); print s")

    pause()
    print ""
    print_fill("Check that taking logarithm returns C/2.")
    print ""
    print_exec("print log(s)")

    pause()
    print ""
    print_fill("Perform a complex rotation.")
    print ""
    print_exec("w = s*z/s; print w")

    pause()
    print ""
    print_fill("Extract the real and imaginary parts.")
    print ""
    print_exec("u = (w^i)/i; print u")
    print_exec("v = (w-u)/i; print v")

    pause()
    print ""
    print_fill("Check that the invariants of the complex rotation are preserved.")
    print ""
    print_exec("print x*x-y*y")
    print_exec("print u*u-v*v")
    print_exec("print x&y")
    print_exec("print u&v")

    pause()
    print ""
    print_fill("CAYLEY ALGEBRA OF OCTONIONS")
    print ""
    print_fill("Define a Cayley product of two paravectors in R+R^(0,7).")
    print ""
    print_exec("f = (1-e({-7,-6,-4}))*(1-e({-6,-5,-3}))*(1-e({-5,-4,-2}))*(1-e({-4,-3,-1})); print f")
    print_exec("o = lambda x,y:   real(x*y*f)+(x*y*f)(1)")

    pause()
    print_exec("print o(e(-7), 2+3*e(-6))")

    pause()
    print_exec("a = 1+2*e(-7)+4*e(-6)+2*e(-5); print a")
    print_exec("b = 2+e(-7)+2*e(-3); print b")

    pause()
    print_exec("c = o(a,b); print c")

    pause()
    print ""
    print_fill("Check that the absolute value is preserved.")
    print ""
    print_exec("print abs(c)")
    print_exec("print abs(a)*abs(b)")

    pause()
    print ""
    print_fill("Consider the associator of three vectors in R^(0,7).")
    print ""
    print_exec("a = 3*e(-7)+2*e(-4)+e(-3); print a")
    print_exec("b = 4*e(-6)+3*e(-4)+e(-3); print b")
    print_exec("c = 3*e(-5)+e(-4)+2*e(-1); print c")
    print_exec("u = o(o(a,b) ,c) - o(a, o(b,c)); print u")

    pause()
    print ""
    print_fill("Show that the associator is perpendicular to the factors.")
    print ""
    print_exec("V = a^b^c; print V")
    print_exec("print (u&V)/V")

    pause()
    print ""
    print_fill("Check the Moufang identity (a o b) o (c o a) == (a o (b o c)) o a.")
    print ""
    print_exec("print o(o(a,b),o(c,a))")
    print_exec("print o(o(a,o(b,c)),a)")

    pause()
    print ""
    print_fill("CONFORMAL TRANSFORMATIONS AND VAHLEN MATRICES")
    print_fill("in the Minkowski space-time R^(3,1).")
    print ""
    print_fill("We choose the two extra dimensions e(-2), e(4)  so that our Minkowski space-time" +
              " has a basis  { e(-1), e(1), e(2), e(3) }.")
    print ""
    print_exec("eneg = e(-2); print eneg")
    print_exec("epos = e(4); print epos")

    pause()
    print ""
    print_fill("We give entries of a 2x2-matrix in the Lie algebra of the Vahlen group.")
    print ""
    print_exec("A = 3+4*e({1,2})+e({2,3})-2*e({-1,3}); print A")
    print_exec("B = 3*e(-1)+7*e(1)+e(2); print B")
    print_exec("C = e(-1)+5*e(2); print C")
    print_exec("D = -reverse(A); print D")

    pause()
    print ""
    print_fill("We give a matrix basis of the Lie algebra of the Vahlen group.")
    print ""
    print_exec("mat11 = (1+(epos^eneg))/2; print mat11")
    print_exec("mat12 = (epos-eneg)/2; print mat12")
    print_exec("mat21 = (epos+eneg)/2; print mat21")
    print_exec("mat22 = (1-(epos^eneg))/2; print mat22")

    pause()
    print ""
    print_fill("We form the following VAHLEN matrix" +
              " in the Lie algebra of the Vahlen group")
    print ""
    print_exec("VAHLEN = A*mat11+B*mat12+involute(C)*mat21+involute(D)*mat22; print VAHLEN")

    pause()
    print ""
    print_fill("By exponentiation we get the vahlen matrix.")
    print ""
    print_exec("vahlen = exp(VAHLEN/2); print vahlen")

    pause()
    print_exec("Vah11 = (vahlen^epos^eneg)/(epos^eneg); print Vah11")
    print_exec("Vah12 = ((vahlen^eneg)/eneg-Vah11)/epos; print Vah12")
    print_exec("Vah21 = ((vahlen^epos)/epos-Vah11)/eneg; print Vah21")
    print_exec("Vah22 = (vahlen-Vah11-Vah12*epos-Vah21*eneg)/(epos^eneg); print Vah22")

    pause()
    print ""
    print_fill("Here are the entries of the vahlen matrix.")
    print ""
    print_exec("a = Vah11+Vah22; print a")
    print_exec("b = Vah12-Vah21; print b")
    print_exec("c = involute(Vah12+Vah21); print c")
    print_exec("d = involute(Vah11-Vah22); print d")

    pause()
    print ""
    print_fill("Create a Mobius transformation.")
    print ""
    print_exec("g = lambda x:  (a*x+b)/(c*x+d)")

    pause()
    print_exec("x = 4*e(-1)+6*e(1)+e(2)+3*e(3); print x")
    print_exec("y = 5*e(-1)+3*e(1)+4*e(2)+e(3); print y")

    pause()
    print_exec("print g(x)")
    print_exec("print g(y)")

    pause()
    print ""
    print_fill("Check the formula g(x)-g(y) == 1/reverse(c*x+d)*(x-y)/(c*y+d).")
    print ""
    print_exec("print g(x)-g(y)")
    print_exec("print 1/reverse(c*x+d)*(x-y)/(c*y+d)")

    pause()
    print ""
    print_fill("ELEMENTARY FUNCTIONS IN THE EXTERIOR ALGEBRA")
    print_fill("of the Minkowski space-time R^(3,1).")
    print ""
    print_exec("x = 3*e(-1)+3*e(1)+4*e(2)+2*e(3); print x")
    print_exec("B = pi/2*e({1,2})-0.5*e({-1,3}); print B")

    pause()
    print ""
    print_fill("Rotate x by the angle pi/2 and perform a boost in the direction e(3).")
    print ""
    print_exec("u = exp(B/2); print u")
    print_exec("y = u*x/u; print y")

    pause()
    print_exec("print y*y")
    print_exec("print x*x")

    pause()
    print ""
    print_fill("Compare the above ordinary exponential exp(A) = 1+A+AA/2+AAA/6+...")
    print_fill("to the following exterior exponential expo(B) = 1+B+B^B/2+B^B^B/6:")
    print ""
    print_exec("expo = lambda B: 1+B+(B^B)/2+(B^B^B)/6")
    print_exec("v = expo(B); print v")

    pause()
    print_exec("z = v*x/v; print z")
    print_exec("print z*z")
    print_exec("print x*x")

    pause()
    print ""
    print_fill("Check that ordinary logarithm and exterior logarithm return B/2 and B.")
    print ""
    print_exec("print log(u)")

    pause()
    print_exec("logi = lambda x:  x - (x^x)/2 + (x^x^x)/3")
    print_exec("logo = lambda x:  log(real(x)) + logi(x/real(x)-1)")

    pause()
    print_exec("print logo(v)")

    pause()
    print_exec("print u*reverse(u)")
    print_exec("print v*reverse(v)")

    pause()
    print ""
    print_fill("PURE SPINORS")
    print ""
    print_exec("v = (e(-4)+e(1))*(e(-3)+e(2))*(e(-2)+e(3))/8; print v")
    print_exec("f = (1-e({-4,1}))*(1-e({-3,2}))*(1-e({-2,3}))/8; print f")
    print ""
    print_fill("Here v == e({1,2,3})*f as the following computation shows:")
    print ""
    print_exec("print e({1,2,3})*f")

    pause()
    print ""
    print_fill("The following bivectors annul v:")
    print_fill("{-4,1}-{-3,2}, {-3,2}-{-2,3},")
    print_fill("{-3,1}-{1,2}, {-3,1}+{-4,2}, {-4,-1}-{-1,1},")
    print_fill("{-2,1}-{1,3}, {-2,1}+{-4,3}, {-3,-1}-{-1,2},")
    print_fill("{-2,2}-{2,3}, {-2,2}+{-3,3}, {-2,-1}-{-1,3},")
    print_fill("and their exponential stabilizes v.")
    print_fill("As an example, take a linear combination of the above bivectors:")
    print ""
    print_exec("B = 3*(e({-2,-1})-e({-1,3}))-4*(e({-4,1})-e({-3,2}))-7*(e({-3,1})+e({-4,2})); print B")
    print ""
    print_fill("and exponentiate:")
    print ""
    print_exec("s = exp(B/10); print s")

    pause()
    print ""
    print_fill("Check that s*v == v:")
    print ""
    print_exec("print s*v")
    print_exec("print v")

    pause()
    print ""
    print_fill("Take an arbitrary bivector:")
    print ""
    print_exec("F = 2*e({-3,-2})+3*e({-4,-1})-7*e({-3,2})+4*e({1,2})-5*e({-2,3}); print F")
    print ""
    print_fill("and exponentiate:")
    print ""
    print_exec("u = exp(F/10); print u")

    pause()
    print ""
    print_fill("The product u*v is a pure spinor.")
    print ""
    print_exec("w = u*v")

    pause()
    print ""
    print_fill("Check that reverse(w)*Ak*w vanishes for a pure spinor w and a k-vector Ak" +
              " when k == 0,1,2 but not when k == 3  (v and reverse(w)*A3*w are parallel)")
    print ""
    print_exec("A1 = 2*e(1)+3*e(2)+5*e(3)+7*e(-4)+e(-3)+4*e(-2); print A1")
    print_exec("A2 = e({-4,-3})-4*e({-3,1})-3*e({-4,2})+3*e({1,2})-5*e({-4,3})-e({-2,3})+2*e({2,3}); print A2")
    print_exec("A3 = e({-4,-3,-2})+2*e({-4,-3,2})+4*e({-2,1,2})+e({-4,-3,3})+5*e({-3,-2,3})+3*e({-4,2,3})+e({1,2,3}); print A3")

    pause()
    print_exec("print reverse(w)*w")

    pause()
    print_exec("print reverse(w)*A1*w")

    pause()
    print_exec("print reverse(w)*A2*w")

    pause()
    print_exec("print reverse(w)*A3*w")
    print_exec("print v")

    pause()
    print ""
    print_fill("You have completed the demonstration file pyclical_demo.py.")

if __name__ == "__main__":
    try:
        demo()
    except:
        print_fill("The demo was interrupted.")
        pass
