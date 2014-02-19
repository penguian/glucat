# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# pyclical_tutorial_0_1_index_sets.py:
#
# This file contains a tutorial that explains how to manipulate index sets
# within PyClical.
#
#    copyright            : (C) 2012-2014 by Paul C. Leopardi
#
# Licensed under CC BY-SA 3.0 http://creativecommons.org/licenses/by-sa/3.0/

from pyclical_tutorial_utils import *

def run(ctx):
    ctx_methods = get_object_methods(ctx)
    for method in ctx_methods:
        exec(method + " = ctx_methods['" + method +"']")

    print_head("0.1 Index sets.")
    print_line()
    print_fill("This file contains a tutorial that explains how to manipulate index sets" +
              " within PyClical.")
    print_line()
    print_fill("It is recommended that you do the tutorials in order, beginning with" +
              " 0.0 Notation.")
    print_line()
    print_exec("from PyClical import *")

    pause()
    print_line()
    print_fill("PyClical is based on two Python classes, clifford, and index_set.")
    print_line()
    print_fill("An index set is essentially a finite set of non-zero integers.")
    print_line()

    pause()
    print_line()
    print_fill("Examples: Initialization.")
    print_line()
    print_fill("Because an index set is a set, you can list the elements in any order" +
              " when creating an index set.")
    print_line()
    print_exec("s = index_set({-2,3,4}); print s")
    print_exec("print repr(s)")
    print_exec("t = index_set({4,3,-2}); print t")
    print_exec("print repr(t)")
    print_line()
    print_fill("Note that {-2,3,4} is a Python set, but index_set({-2,3,4})" +
              " is a PyClical index set.")
    print_line()
    print_exec("print type({-2,3,4})")
    print_exec("print type(index_set({-2,3,4}))")

    pause()
    print_line()
    print_fill("You can also create an index set from a string," +
              " either with or without curly braces.")
    print_line()
    print_exec("s = index_set('{-2,3,4}'); print s")
    print_exec("print repr(s)")
    print_exec("t = index_set('4,3,-2'); print t")
    print_exec("print repr(t)")

    pause()
    print_line()
    print_fill("For convenience, PyClical defines some abbreviations.")
    print_line()
    print_fill("For example, ist is an abbreviation for index_set,")
    print_line()
    print_exec("r = ist(-2); print r")
    print_exec("s = ist({-2,3,4}); print s")
    print_line()
    print_fill("and istpq(p,q) is an abbreviation for index_set({-q,...,-1,1,...,p}):")
    print_line()
    print_exec("t = istpq(4,1); print t")
    print_exec("u = istpq(0,4); print u")

    pause()
    check_exec("set s to be the index set {-1,4,7}.",
               "s",
               "ist({-1,4,7})")

    pause()
    print_line()
    print_fill("Examples: Operations.")
    print_line()
    print_fill("The notation to test for set membership uses subscripting:")
    print_line()
    print_exec("print index_set(1)[1]")
    print_exec("print index_set(1)[2]")
    print_line()
    print_fill("This is also used for the notation to alter the membership of an index set.")
    print_line()
    print_exec("s=index_set(1); s[2] = True; print s")
    print_exec("s=index_set({1,2}); s[1] = False; print s")

    pause()
    print_line()
    print_fill("You can also use the usual Python 'in' notation for set membership:")
    print_line()
    print_exec("print 1 in index_set(1)")
    print_exec("print 2 in index_set(1)")

    pause()
    print_line()
    print_fill("Because an index set is essentially a set, the usual unary and binary" +
              " Boolean operations are defined for index sets. " +
              " These are as follows.")
    print_line()
    print_fill("Set complement: not is denoted by ~. " +
              " The negation is with respect to the maximum index set. " +
              " On a 32 bit machine, this is istpq(16,16). " +
              " On a 64 bit machine, this is istpq(32,32).")
    print_line()
    print_exec("print ~istpq(16,16)")

    pause()
    print_line()
    print_fill("Symmetric set difference: exclusive or is denoted by ^.")
    print_line()
    print_exec("print ist(1) ^ ist(2)")
    print_exec("print ist({1,2}) ^ ist(2)")
    print_exec("s = ist({1,2}); s ^= ist(2); print s")

    pause()
    print_line()
    print_fill("Set intersection: and is denoted by &.")
    print_line()
    print_exec("print ist(1) & ist(2)")
    print_exec("print ist({1,2}) & ist(2)")
    print_exec("s = ist({1,2}); s &= ist(2); print s")

    pause()
    print_line()
    print_fill("Set union: or is denoted by |.")
    print_line()
    print_exec("print ist(1) | ist(2)")
    print_exec("print ist({1,2}) | ist(2)")
    print_exec("s = ist({1,2}); s |= ist(2); print s")

    pause()
    check_exec("set s to be the symmetric set difference between the index sets {1,3} and {3,5}.",
               "s",
               "ist({1,3}) ^ ist({3,5})")
    pause()
    check_exec("set s to be the intersection of the index sets {1,3} and {3,5}.",
               "s",
               "ist({1,3}) & ist({3,5})")
    pause()
    check_exec("set s to be the union of the index sets {1,3} and {3,5}.",
               "s",
               "ist({1,3}) | ist({3,5})")

    pause()
    print_line()
    print_fill("Examples: Iteration.")
    print_line()
    print_fill("Iteration over the elements of an index set is supported.")
    print_exec("s = ist({1,2}); t = ist({-1,1})")
    print_exec("for i in s: print i, t[i]")

    pause()
    print_line()
    print_fill("Examples: Member functions.")
    print_line()
    print_fill("A number of member functions give information about an index set.")
    print_line()
    print_fill("Here is a list.")
    print_line()
    print_fill("count(): Cardinality: Number of indices included in set.")
    print_line()
    print_exec("print index_set({-1,1,2}).count()")
    print_line()
    print_fill("count_neg(): Number of negative indices included in set.")
    print_line()
    print_exec("print index_set({-1,1,2}).count_neg()")
    print_line()
    print_fill("count_pos(): Number of positive indices included in set.")
    print_line()
    print_exec("print index_set({-1,1,2}).count_pos()")

    pause()
    print_line()
    print_fill("min(): Minimum member.")
    print_line()
    print_exec("print index_set({-1,1,2}).min()")
    print_line()
    print_fill("max(): Maximum member.")
    print_line()
    print_exec("print index_set({-1,1,2}).max()")

    pause()
    print_line()
    print_fill("sign_of_mult(t): Sign of geometric product of two clifford basis elements.")
    print_line()
    print_exec("s = index_set({1,2}); t=index_set(-1); print s.sign_of_mult(t)")
    print_line()
    print_fill("Sign of geometric square of a Clifford basis element.")
    print_line()
    print_exec("s = index_set({1,2}); print s.sign_of_square()")

    pause()
    print_line()
    print_fill("Examples: Ordinary functions.")
    print_line()
    print_fill("There are also a few ordinary functions that operate on index sets.")
    print_line()
    print_fill("compare(lhs,rhs): 'lexicographic compare' eg. {3,4,5} is less than {3,7,8};" +
             " -1 if a<b, +1 if a>b, 0 if a==b.")
    print_line()
    print_exec("print compare(index_set({1,2}),index_set({-1,3}))")
    print_exec("print compare(index_set({-1,4}),index_set({-1,3}))")

    pause()
    print_line()
    print_fill("min_neg(obj): Minimum negative index, or 0 if none.")
    print_line()
    print_exec("print min_neg(index_set({1,2}))")
    print_line()
    print_fill("max_pos(obj): Maximum positive index, or 0 if none.")
    print_line()
    print_exec("print max_pos(index_set({1,2}))")

    pause()
    print_line()
    print_fill("You have completed the tutorial file pyclical_tutorial_0_1_index_sets.py.")

if __name__ == "__main__":
    ctx = tutorial_context(globals())
    try:
        run(ctx)
    except:
        ctx.print_fill("The tutorial was interrupted.")
        pass
