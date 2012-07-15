# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_0_0_notation.py:
#    This file contains a tutorial that explains the notation used by PyClical.
#
#    copyright            : (C) 2012 by Paul C. Leopardi
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

def tutorial():
    ctx = interaction_context(globals())
    print_exec = ctx.print_exec
    check_exec = ctx.check_exec

    print "# pyclical_tutorial_0_0_notation.py:"
    print ""
    print_fill("0.0 Notation.")
    print ""
    print_fill("This file contains a tutorial that explains the notation used by PyClical.")
    print ""

    pause()
    print ""
    print_fill("To use the capabilities of PyClical from within Python, you must either" +
              " import the PyClical extension module or import objects from this module. " +
              " The simplest way to do this is to use the following Python statement:")
    print ""
    print_exec("from PyClical import *")

    pause()
    print ""
    print_fill("PyClical is based on two Python classes, clifford, and index_set.")
    print ""
    print_fill("The clifford class supports all of the operations usually needed in Clifford" +
              " and Grassmann algebras.")
    print ""
    print_fill("The index_set class plays an important supporting role.")
    print ""

    pause()
    print ""
    print_fill("Every clifford object class can be represented as a linear combination of" +
              " basis elements, with each basis element represented as an index set.")
    print ""
    print_fill("Examples:")
    print ""
    print_exec("v = clifford('{1}+4{2}+3{3}'); print v")
    print_exec("w = clifford('1+2{2}+3{1,2}'); print w")

    pause()
    print ""
    print_fill("The index sets are displayed in canonical increasing order of the indices," +
              " corresponding exactly to the multiplication of generators of the Clifford" +
              " algebra in the same canonical increasing order.")
    print ""
    print_fill("Examples:")
    print ""
    print_exec("w = clifford(index_set({-2,3,4})); print w")
    print_exec("print repr(w)")
    print_exec("x = clifford(index_set({-2})) * clifford(index_set({3})) * clifford(index_set({4})); print x")
    print_exec("print repr(x)")

    pause()
    print ""
    print_fill("Generators with a negative index square to -1.")
    print_fill("Generators with a positive index square to  1.")
    print ""
    print_fill("Examples:")
    print ""
    print_exec("print clifford(index_set({-1})) ** 2")
    print_exec("print clifford(index_set({1})) ** 2")

    pause()
    print ""
    print_fill("Because an index set is a set, you can list the elements in any order" +
              " when creating an index set.")
    print ""
    print_fill("Examples:")
    print ""
    print_exec("s = index_set({-2,3,4}); print s")
    print_exec("print repr(s)")
    print_exec("t = index_set({4,3,-2}); print t")
    print_exec("print repr(t)")

    pause()
    print ""
    print_fill("For convenience, PyClical defines some abbreviations.")
    print ""
    print_fill("For example, ist is an abbreviation for index_set,")
    print ""
    print_exec("r = ist(-2); print r")
    print_exec("s = ist({-2,3,4}); print s")
    print ""
    print_fill("and istpq(p,q) is an abbreviation for index_set({-q,...,-1,1,...,p}):")
    print ""
    print_exec("t = istpq(4,1); print t")
    print_exec("u = istpq(0,4); print u")

    pause()
    print ""
    print_fill("The abbreviation e(obj) is often used for clifford(index_set(obj)),")
    print ""
    print_exec("w = e({4,3,-2}); print w")
    print_exec("x = e(-2) * e(3) * e(4); print x")
    print_exec("y = 10 + 2*e(-2) + e({3,4}); print y")
    print_exec("print e(-1)**2")
    print ""
    print_fill("and the abbreviation cl can be used for clifford.")
    print ""
    print_exec("z = cl('10+2{-2}+{3,4}'); print z")

    pause()
    print ""
    print_fill("To obtain the coordinate of a clifford object x with respect to an index set s," +
              " use subscripting notation x[s]. " +
              " You can also use an index, a set of indices, or a string as the subscript,"
              " as long as the subscript defines an index set.")
    print ""
    print_fill("Examples:")
    print ""
    print_exec("y = 10 + 2*e(-2) + e({3,4}); print y")
    print_exec("print y[0]")
    print_exec("print y[index_set('{-2}')]")
    print_exec("print y[ist(-2)]")
    print_exec("print y[-2]")
    print_exec("print y[{3,4}]")
    print_exec("print y['{3,4}']")

    pause()
    print ""
    print_fill("There is a special notation to specify vectors, using a tuple or list of" +
              " coordinates and a basis. " +
              " The basis is specified using an index set.")
    print ""
    print_fill("Examples:")
    print ""
    print_exec("v = clifford((4,1,3), ist({-5,2,7})); print v")
    print_exec("w = clifford([1,5,4,-7], istpq(1,3)); print w")

    pause()
    check_exec("set s to be the index_set {-1,3,4}.",
               "s",
               "ist({-1,3,4})")
    check_exec("set v to be the vector 2{-1}+3{1}.",
               "v",
               "2*e(-1)+3*e(1)")
    check_exec("set x to be the multivector 3{1}+4{2,3}.",
               "x",
               "3*e(1)+4*e({2,3})")
               
    pause()
    print ""
    print_fill("You can also create a random element of a given Clifford algebra. " +
              " The algebra is specified using an index set.")
    print ""
    print_fill("Examples:")
    print ""
    print_exec("v = random_clifford(ist({-5,2,7})); print v")
    print_exec("w = random_clifford(istpq(1,3)); print w")

    pause()
    print ""
    print_fill("You have completed the tutorial file pyclical_tutorial_0_0_notation.py.")

if __name__ == "__main__":
    tutorial()
