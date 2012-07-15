# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_0_1_index_sets.py:
#    This file contains a tutorial that explains how to manipulate index sets
#    within PyClical.
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

    print "# pyclical_tutorial_0_1_index_sets.py:"
    print ""
    print_fill("0.1 Index sets.")
    print ""
    print_fill("This file contains a tutorial that explains how to manipulate index sets" +
              " within PyClical.")
    print ""
    print_fill("It is recommended that you do the tutorials in order, beginning with" +
              " 0.0 Notation.")
    print ""

    pause()
    print ""
    print_fill("PyClical is based on two Python classes, clifford, and index_set.")
    print ""
    print_fill("An index set is essentially a finite set of non-zero integers.")
    print ""

    pause()
    print ""
    print_fill("Examples: Initialization.")
    print ""
    print_fill("Because an index set is a set, you can list the elements in any order" +
              " when creating an index set.")
    print ""
    print_exec("s = index_set({-2,3,4}); print s")
    print_exec("print repr(s)")
    print_exec("t = index_set({4,3,-2}); print t")
    print_exec("print repr(t)")
    print ""
    print_fill("Note that {-2,3,4} is a Python set, but index_set({-2,3,4})" +
              " is a PyClical index set.")
    print ""
    print_exec("print type({-2,3,4})")
    print_exec("print type(index_set({-2,3,4}))")

    pause()
    print ""
    print_fill("You can also create an index set from a string," +
              " either with or without curly braces.")
    print ""
    print_exec("s = index_set('{-2,3,4}'); print s")
    print_exec("print repr(s)")
    print_exec("t = index_set('4,3,-2'); print t")
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
    check_exec("set s to be the index set {-1,4,7}.",
               "s",
               "ist({-1,4,7})")

    pause()
    print ""
    print_fill("Examples: Operations.")
    print ""
    print_fill("The notation to test for set membership uses subscripting:")
    print ""
    print_exec("print index_set(1)[1]")
    print_exec("print index_set(1)[2]")
    print ""
    print_fill("This is also used for the notation to alter the membership of an index set.")
    print ""
    print_exec("s=index_set(1); s[2] = True; print s")
    print_exec("s=index_set({1,2}); s[1] = False; print s")

    pause()
    print ""
    print_fill("Because an index set is essentially a set, the usual unary and binary" +
              " Boolean operations are defined for index sets. " +
              " These are as follows.")
    print ""
    print_fill("Set complement: not is denoted by ~. " +
              " The negation is with respect to the maximum index set. " +
              " On a 32 bit machine, this is istpq(16,16). " +
              " On a 64 bit machine, this is istpq(32,32).")
    print ""
    print_exec("print ~istpq(16,16)")

    pause()
    print ""
    print_fill("Symmetric set difference: exclusive or is denoted by ^.")
    print ""
    print_exec("print ist(1) ^ ist(2)")
    print_exec("print ist({1,2}) ^ ist(2)")
    print_exec("s = ist({1,2}); s ^= ist(2); print s")

    pause()
    print ""
    print_fill("Set intersection: and is denoted by &.")
    print ""
    print_exec("print ist(1) & ist(2)")
    print_exec("print ist({1,2}) & ist(2)")
    print_exec("s = ist({1,2}); s &= ist(2); print s")

    pause()
    print ""
    print_fill("Set union: or is denoted by |.")
    print ""
    print_exec("print ist(1) | ist(2)")
    print_exec("print ist({1,2}) | ist(2)")
    print_exec("s = ist({1,2}); s |= ist(2); print s")

    pause()
    check_exec("set s to be the symmetric set difference between the index sets {1,3} and {3,5}.",
               "s",
               "ist({1,3}) ^ ist({3,5})")
    check_exec("set s to be the intersection of the index sets {1,3} and {3,5}.",
               "s",
               "ist({1,3}) & ist({3,5})")
    check_exec("set s to be the union of the index sets {1,3} and {3,5}.",
               "s",
               "ist({1,3}) | ist({3,5})")

    pause()
    print ""
    print_fill("Examples: Member functions.")
    print ""
    print_fill("A number of member functions give information about an index set.")
    print ""
    print_fill("Here is a list.")
    print ""
    print_fill("count(): Cardinality: Number of indices included in set.")
    print ""
    print_exec("print index_set({-1,1,2}).count()")
    print ""
    print_fill("count_neg(): Number of negative indices included in set.")
    print ""
    print_exec("print index_set({-1,1,2}).count_neg()")
    print ""
    print_fill("count_pos(): Number of positive indices included in set.")
    print ""
    print_exec("print index_set({-1,1,2}).count_pos()")

    pause()
    print ""
    print_fill("min(): Minimum member.")
    print ""
    print_exec("print index_set({-1,1,2}).min()")
    print ""
    print_fill("max(): Maximum member.")
    print ""
    print_exec("print index_set({-1,1,2}).max()")

    pause()
    print ""
    print_fill("sign_of_mult(t): Sign of geometric product of two clifford basis elements.")
    print ""
    print_exec("s = index_set({1,2}); t=index_set(-1); print s.sign_of_mult(t)")
    print ""
    print_fill("Sign of geometric square of a Clifford basis element.")
    print ""
    print_exec("s = index_set({1,2}); print s.sign_of_square()")

    pause()
    print ""
    print_fill("Examples: Ordinary functions.")
    print ""
    print_fill("There are also a few ordinary functions that operate on index sets.")
    print ""
    print_fill("compare(lhs,rhs): 'lexicographic compare' eg. {3,4,5} is less than {3,7,8};" +
             " -1 if a<b, +1 if a>b, 0 if a==b.")
    print ""
    print_exec("print compare(index_set({1,2}),index_set({-1,3}))")
    print_exec("print compare(index_set({-1,4}),index_set({-1,3}))")

    pause()
    print ""
    print_fill("min_neg(obj): Minimum negative index, or 0 if none.")
    print ""
    print_exec("print min_neg(index_set({1,2}))")
    print ""
    print_fill("max_pos(obj): Maximum positive index, or 0 if none.")
    print ""
    print_exec("print max_pos(index_set({1,2}))")

    pause()
    print ""
    print_fill("You have completed the tutorial file pyclical_tutorial_0_1_index_sets.py.")

if __name__ == "__main__":
    tutorial()
