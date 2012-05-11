# -*- coding: utf-8 -*-
# distutils: language = c++
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# PyClical.pyx: Cython definitions visible from Python.
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

import math
import numbers
import collections

from PyClical cimport *

# Forward reference
cdef class index_set 

cdef inline IndexSet toIndexSet(obj):
    """
    Return the C++ IndexSet instance wrapped by index_set(obj).
    """
    return index_set(obj).instance[0]

cdef class index_set:
    """
    Python class index_set wraps C++ class IndexSet.
    """
    cdef IndexSet *instance # Wrapped instance of C++ class IndexSet.

    cdef inline wrap(index_set self, IndexSet other):
        """
        Wrap an instance of the C++ class IndexSet.
        """
        self.instance[0] = other
        return self

    cdef inline IndexSet unwrap(index_set self):
        """
        Return the wrapped C++ IndexSet instance.
        """
        return self.instance[0]

    cpdef copy(index_set self):
        """
        Copy this index_set object.

        >>> s=index_set(1); t=s.copy(); print t
        {1}
        """
        return index_set(self)

    def __cinit__(self, other = 0):
        """
        Construct and object of type index_set.

        >>> print index_set(1)
        {1}
        >>> print index_set("{1,2}")
        {1,2}
        >>> print index_set(index_set("{1,2}"))
        {1,2}
        >>> print index_set({1,2})
        {1,2}
        >>> print index_set({1,2,1})
        {1,2}
        >>> print index_set("{1,2,1}")
        {1,2}
        >>> print index_set("")
        {}
        """
        error_msg_prefix = "Cannot initialize index_set object from"
        if   isinstance(other, index_set):
            self.instance = new IndexSet((<index_set>other).unwrap())
        elif isinstance(other, numbers.Integral):
            self.instance = new IndexSet(<int>other)
        elif isinstance(other, (set, frozenset)):
            try:
                self.instance = new IndexSet()
                for idx in other:
                    self[idx] = True
            except IndexError:
                raise IndexError(error_msg_prefix + " invalid " + repr(other) + ".")
            except (RuntimeError, TypeError):
                raise ValueError(error_msg_prefix + " invalid " + repr(other) + ".")
        elif isinstance(other, str):
            try:
                self.instance = new IndexSet(<char *>other)
            except RuntimeError:
                raise ValueError(error_msg_prefix + " invalid string " + repr(other) + ".")
        else:
            raise TypeError(error_msg_prefix + " " + str(type(other)) + ".")

    def __dealloc__(self):
        """
        Clean up by deallocating the instance of C++ class IndexSet.
        """
        del self.instance

    def __richcmp__(lhs, rhs, int op):
        """
        Compare two objects of class index_set.

        >>> index_set(1) == index_set("{1}")
        True
        >>> index_set("{1}") != index_set("{1}")
        False
        >>> index_set("{1}") != index_set("{2}")
        True
        >>> index_set("{1}") == index_set("{2}")
        False
        >>> index_set("{1}") < index_set("{2}")
        True
        >>> index_set("{1}") <= index_set("{2}")
        True
        >>> index_set("{1}") > index_set("{2}")
        False
        >>> index_set("{1}") >= index_set("{2}")
        False
        """
        if (lhs is None) or (rhs is None):
            eq = bool(lhs is rhs)
            if op == 2: # ==
                return eq
            elif op == 3: # !=
                return not eq
            else:
                if op == 0: # <
                    return False
                elif op == 1: # <=
                    return eq
                elif op == 4: # >
                    return False
                elif op == 5: # >=
                    return eq
                else:
                    return NotImplemented
        else:
            eq = bool( toIndexSet(lhs) == toIndexSet(rhs) )
            if op == 2: # ==
                return eq
            elif op == 3: # !=
                return not eq
            else:
                lt = bool( toIndexSet(lhs) < toIndexSet(rhs) )
                if op == 0: # <
                    return lt
                elif op == 1: # <=
                    return lt or eq
                elif op == 4: # >
                    return not (lt or eq)
                elif op == 5: # >=
                    return not lt
                else:
                    return NotImplemented

    def __getitem__(self, idx):
        """
        Get the value of an index_set object at an index.

        >>> index_set("{1}")[1]
        True
        >>> index_set("{1}")[2]
        False
        >>> index_set("{2}")[-1]
        False
        >>> index_set("{2}")[1]
        False
        >>> index_set("{2}")[2]
        True
        >>> index_set("{2}")[33]
        False
        """
        return self.instance.getitem(idx)

    def __setitem__(self, idx, val):
        """
        Set the value of an index_set object at index idx to value val.

        >>> s=index_set("{1}"); s[2] = True; print s
        {1,2}
        >>> s=index_set("{1,2}"); s[1] = False; print s
        {2}
        """
        self.instance.set(idx, val)
        return

    def __invert__(self):
        """
        Set complement: not.

        >>> print ~index_set("{-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16}")
        {-32,-31,-30,-29,-28,-27,-26,-25,-24,-23,-22,-21,-20,-19,-18,-17,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32}
        """
        return index_set().wrap( self.instance.invert() )

    def __xor__(lhs, rhs):
        """
        Symmetric set difference: exclusive or.

        >>> print index_set("{1}") ^ index_set("{2}")
        {1,2}
        >>> print index_set("{1,2}") ^ index_set("{2}")
        {1}
        """
        return index_set().wrap( toIndexSet(lhs) ^ toIndexSet(rhs) )

    def __ixor__(self, rhs):
        """
        Symmetric set difference: exclusive or.

        >>> x = index_set("{1}"); x ^= index_set("{2}"); print x
        {1,2}
        >>> x = index_set("{1,2}"); x ^= index_set("{2}"); print x
        {1}
        """
        return self.wrap( self.unwrap() ^ toIndexSet(rhs) )

    def __and__(lhs, rhs):
        """
        Set intersection: and.

        >>> print index_set("{1}") & index_set("{2}")
        {}
        >>> print index_set("{1,2}") & index_set("{2}")
        {2}
        """
        return index_set().wrap( toIndexSet(lhs) & toIndexSet(rhs) )

    def __iand__(self, rhs):
        """
        Set intersection: and.

        >>> x = index_set("{1}"); x &= index_set("{2}"); print x
        {}
        >>> x = index_set("{1,2}"); x &= index_set("{2}"); print x
        {2}
        """
        return self.wrap( self.unwrap() & toIndexSet(rhs) )

    def __or__(lhs, rhs):
        """
        Set union: or.

        >>> print index_set("{1}") | index_set("{2}")
        {1,2}
        >>> print index_set("{1,2}") | index_set("{2}")
        {1,2}
        """
        return index_set().wrap( toIndexSet(lhs) | toIndexSet(rhs) )

    def __ior__(self, rhs):
        """
        Set union: or.

        >>> x = index_set("{1}"); x |= index_set("{2}"); print x
        {1,2}
        >>> x = index_set("{1,2}"); x |= index_set("{2}"); print x
        {1,2}
        """
        return self.wrap( self.unwrap() | toIndexSet(rhs) )

    def count(self):
        """
        Cardinality: Number of indices included in set.

        >>> index_set("{-1,1,2}").count()
        3
        """
        return self.instance.count()

    def count_neg(self):
        """
        Number of negative indices included in set.

        >>> index_set("{-1,1,2}").count_neg()
        1
        """
        return self.instance.count_neg()

    def count_pos(self):
        """
        Number of positive indices included in set.

        >>> index_set("{-1,1,2}").count_pos()
        2
        """
        return self.instance.count_pos()

    def min(self):
        """
        Minimum member.

        >>> index_set("{-1,1,2}").min()
        -1
        """
        return self.instance.min()

    def max(self):
        """
        Maximum member.

        >>> index_set("{-1,1,2}").max()
        2
        """
        return self.instance.max()

    def hash_fn(self):
        """
        Hash function.
        """
        return self.instance.hash_fn()

    def sign_of_mult(self, rhs):
        """
        Sign of geometric product of two Clifford basis elements.

        >>> s = index_set("{1,2}"); t=index_set("{-1}"); s.sign_of_mult(t)
        1
        """
        return self.instance.sign_of_mult(toIndexSet(rhs))

    def sign_of_square(self):
        """
        Sign of geometric square of a Clifford basis element.

        >>> s = index_set("{1,2}"); s.sign_of_square()
        -1
        """
        return self.instance.sign_of_square()

    def __repr__(self):
        """
        The “official” string representation of self.

        >>> index_set("{1,2}").__repr__()
        'index_set("{1,2}")'
        >>> repr(index_set("{1,2}"))
        'index_set("{1,2}")'
        """
        return index_set_to_repr( self.unwrap() ).c_str()

    def __str__(self):
        """
        The “informal” string representation of self.

        >>> index_set("{1,2}").__str__()
        '{1,2}'
        >>> str(index_set("{1,2}"))
        '{1,2}'
        """
        return index_set_to_str( self.unwrap() ).c_str()

def index_set_hidden_doctests():
    """
    Tests for functions that Doctest cannot see.
    
    For index_set.__cinit__: Construct index_set.

    >>> print index_set(1)
    {1}
    >>> print index_set("{1,2}")
    {1,2}
    >>> print index_set(index_set("{1,2}"))
    {1,2}
    >>> print index_set({1,2})
    {1,2}
    >>> print index_set({1,2,1})
    {1,2}
    >>> print index_set("{1,2,1}")
    {1,2}
    >>> print index_set("")
    {}
    >>> print index_set("{")
    Traceback (most recent call last):
    ...
    ValueError: Cannot initialize index_set object from invalid string '{'.
    >>> print index_set("{1")
    Traceback (most recent call last):
    ...
    ValueError: Cannot initialize index_set object from invalid string '{1'.
    >>> print index_set("{1,2,100}")
    Traceback (most recent call last):
    ...
    ValueError: Cannot initialize index_set object from invalid string '{1,2,100}'.
    >>> print index_set({1,2,100})
    Traceback (most recent call last):
    ...
    IndexError: Cannot initialize index_set object from invalid set([1, 2, 100]).
    >>> print index_set([1,2])
    Traceback (most recent call last):
    ...
    TypeError: Cannot initialize index_set object from <type 'list'>.

    For index_set.__richcmp__: Compare two objects of class index_set.

    >>> index_set(1) == index_set("{1}")
    True
    >>> index_set("{1}") != index_set("{1}")
    False
    >>> index_set("{1}") != index_set("{2}")
    True
    >>> index_set("{1}") == index_set("{2}")
    False
    >>> index_set("{1}") < index_set("{2}")
    True
    >>> index_set("{1}") <= index_set("{2}")
    True
    >>> index_set("{1}") > index_set("{2}")
    False
    >>> index_set("{1}") >= index_set("{2}")
    False
    >>> None == index_set("{1,2}")
    False
    >>> None != index_set("{1,2}")
    True
    >>> None < index_set("{1,2}")
    False
    >>> None <= index_set("{1,2}")
    False
    >>> None > index_set("{1,2}")
    False
    >>> None >= index_set("{1,2}")
    False
    >>> index_set("{1,2}") == None
    False
    >>> index_set("{1,2}") != None
    True
    >>> index_set("{1,2}") < None
    False
    >>> index_set("{1,2}") <= None
    False
    >>> index_set("{1,2}") > None
    False
    >>> index_set("{1,2}") >= None
    False
    """
    return

cpdef inline compare(lhs,rhs):
    """
    "lexicographic compare" eg. {3,4,5} is less than {3,7,8};
    -1 if a<b, +1 if a>b, 0 if a==b.

    >>> compare(index_set("{1,2}"),index_set("{-1,3}"))
    -1
    >>> compare(index_set("{-1,4}"),index_set("{-1,3}"))
    1
    """
    return glucat.compare( toIndexSet(lhs), toIndexSet(rhs) )

cpdef inline min_neg(obj):
    """
    Minimum negative index, or 0 if none.

    >>> min_neg(index_set("{1,2}"))
    0
    """
    return glucat.min_neg( toIndexSet(obj) )

cpdef inline max_pos(obj):
    """
    Maximum positive index, or 0 if none.

    >>> max_pos(index_set("{1,2}"))
    2
    """
    return glucat.max_pos( toIndexSet(obj) )

cdef inline vector[scalar_t] list_to_vector(lst):
     """
     Create a C++ std:vector[scalar_t] from an iterable Python object.
     """
     cdef vector[scalar_t] v
     for s in lst:
         v.push_back(<scalar_t>s)
     return v

# Forward reference.
cdef class clifford 

cdef inline Clifford toClifford(obj):
    return clifford(obj).instance[0]

cdef class clifford:
    """
    Python class clifford wraps C++ class Clifford.
    """
    cdef Clifford *instance # Wrapped instance of C++ class Clifford.

    cdef inline wrap(clifford self, Clifford other):
        """
        Wrap an instance of the C++ class Clifford.
        """
        self.instance[0] = other
        return self

    cdef inline Clifford unwrap(clifford self):
        """
        Return the wrapped C++ Clifford instance.
        """
        return self.instance[0]

    cpdef copy(clifford self):
        """
        Copy this clifford object.

        >>> x=clifford("1{2}"); y=x.copy(); print y
        {2}
        """
        return clifford(self)

    def __cinit__(self, other = 0, ist = None):
        """
        Construct an object of type clifford.

        >>> print clifford(2)
        2
        >>> print clifford(2L)
        2
        >>> print clifford(2.0)
        2
        >>> print clifford(1.0e-1)
        0.1
        >>> print clifford("2")
        2
        >>> print clifford("2{1,2,3}")
        2{1,2,3}
        >>> print clifford(clifford("2{1,2,3}"))
        2{1,2,3}
        >>> print clifford("-{1}")
        -{1}
        >>> print clifford(2,index_set("{1,2}"))
        2{1,2}
        >>> print clifford([2,3],index_set("{1,2}"))
        2{1}+3{2}
        """
        error_msg_prefix = "Cannot initialize clifford object from"
        if ist is None:
            try:
                if   isinstance(other, clifford):
                    self.instance = new Clifford((<clifford>other).unwrap())
                elif isinstance(other, index_set):
                    self.instance = new Clifford((<index_set>other).unwrap(), <scalar_t>1.0)
                elif isinstance(other, numbers.Real):
                    self.instance = new Clifford(<scalar_t>other)
                elif isinstance(other, str):
                    try:
                        self.instance = new Clifford(<char *>other)
                    except RuntimeError:
                        raise ValueError(error_msg_prefix + " invalid string " + repr(other) + ".")
                else:
                    raise TypeError(error_msg_prefix + " " + str(type(other)) + ".")
            except RuntimeError as runtime_error:
                raise ValueError(error_msg_prefix + " " + str(type(other))
                                                  + " value " + repr(other) + ":"
                                                  + "\n\t" + str(runtime_error))
        elif isinstance(ist, index_set):
            if   isinstance(other, numbers.Real):
                self.instance = new Clifford((<index_set>ist).unwrap(), <scalar_t>other)
            elif isinstance(other, collections.Sequence):
                self.instance = new Clifford(list_to_vector(other), (<index_set>ist).unwrap())
            else:
                raise TypeError(error_msg_prefix + " (" + str(type(other))
                                                 + ", " + repr(ist) + ").")
        else:
            raise TypeError(error_msg_prefix + " (" + str(type(other))
                                             + ", " + str(type(ist)) + ").")

    def __dealloc__(self):
        """
        Clean up by deallocating the instance of C++ class Clifford.
        """
        del self.instance

    def __richcmp__(lhs, rhs, int op):
        """
        Compare objects of type clifford.

        >>> clifford("{1}") == clifford("1{1}")
        True
        >>> clifford("{1}") != clifford("1.0{1}")
        False
        >>> clifford("{1}") != clifford("1.0")
        True
        >>> clifford("{1,2}") == None
        False
        >>> clifford("{1,2}") != None
        True
        >>> None == clifford("{1,2}")
        False
        >>> None != clifford("{1,2}")
        True
        """
        if op == 2: # ==
            if (lhs is None) or (rhs is None):
                return bool(lhs is rhs)
            else:          
                return bool( toClifford(lhs) == toClifford(rhs) )
        elif op == 3: # !=
            if (lhs is None) or (rhs is None):
                return not bool(lhs is rhs)
            else:
                return bool( toClifford(lhs) != toClifford(rhs) )
        else:
            return NotImplemented

    def __getitem__(self, ist):
        """
        Subscripting: map from index set to scalar coordinate.

        >>> clifford("{1}")[index_set("1")]
        1.0
        >>> clifford("{1}")[index_set("{1}")]
        1.0
        >>> clifford("{1}")[index_set("{1,2}")]
        0.0
        >>> clifford("2{1,2}")[index_set("{1,2}")]
        2.0
        """
        return self.instance.getitem(toIndexSet(ist))

    def __neg__(self):
        """
        Unary -.

        >>> print -clifford("{1}")
        -{1}
        """
        return clifford().wrap( self.instance.neg() )

    def __pos__(self):
        """
        Unary +.

        >>> print +clifford("{1}")
        {1}
        """
        return clifford(self)

    def __add__(lhs, rhs):
        """
        Geometric sum.

        >>> print clifford(1) + clifford("{2}")
        1+{2}
        >>> print clifford("{1}") + clifford("{2}")
        {1}+{2}
        """
        return clifford().wrap( toClifford(lhs) + toClifford(rhs) )

    def __iadd__(self, rhs):
        """
        Geometric sum.

        >>> x = clifford(1); x += clifford("{2}"); print x
        1+{2}
        """
        return self.wrap( self.unwrap() + toClifford(rhs) )

    def __sub__(lhs, rhs):
        """
        Geometric difference.

        >>> print clifford(1) - clifford("{2}")
        1-{2}
        >>> print clifford("{1}") - clifford("{2}")
        {1}-{2}
        """
        return clifford().wrap( toClifford(lhs) - toClifford(rhs) )

    def __isub__(self, rhs):
        """
        Geometric difference.

        >>> x = clifford(1); x -= clifford("{2}"); print x
        1-{2}
        """
        return self.wrap( self.unwrap() - toClifford(rhs) )

    def __mul__(lhs, rhs):
        """
        Geometric product.

        >>> print clifford("{1}") * clifford("{2}")
        {1,2}
        >>> print clifford(2) * clifford("{2}")
        2{2}
        >>> print clifford("{1}")*clifford("{1,2}")
        {2}
        """
        return clifford().wrap( toClifford(lhs) * toClifford(rhs) )

    def __imul__(self, rhs):
        """
        Geometric product.

        >>> x = clifford(2); x *= clifford("{2}"); print x
        2{2}
        >>> x = clifford("{1}"); x *= clifford("{2}"); print x
        {1,2}
        >>> x = clifford("{1}"); x *= clifford("{1,2}"); print x
        {2}
        """
        return self.wrap( self.unwrap() * toClifford(rhs) )

    def __mod__(lhs, rhs):
        """
        Contraction.

        >>> print clifford("{1}") % clifford("{2}")
        0
        >>> print clifford(2) % clifford("{2}")
        2{2}
        >>> print clifford("{1}") % clifford("{1}")
        1
        >>> print clifford("{1}") % clifford("{1,2}")
        {2}
        """
        return clifford().wrap( toClifford(lhs) % toClifford(rhs) )

    def __imod__(self, rhs):
        """
        Contraction.

        >>> x = clifford("{1}"); x %= clifford("{2}"); print x
        0
        >>> x = clifford(2); x %= clifford("{2}"); print x
        2{2}
        >>> x = clifford("{1}"); x %= clifford("{1}"); print x
        1
        >>> x = clifford("{1}"); x %= clifford("{1,2}"); print x
        {2}
        """
        return self.wrap( self.unwrap() % toClifford(rhs) )

    def __and__(lhs, rhs):
        """
        Inner product.

        >>> print clifford("{1}") & clifford("{2}")
        0
        >>> print clifford(2) & clifford("{2}")
        0
        >>> print clifford("{1}") & clifford("{1}")
        1
        >>> print clifford("{1}") & clifford("{1,2}")
        {2}
        """
        return clifford().wrap( toClifford(lhs) & toClifford(rhs) )

    def __iand__(self, rhs):
        """
        Inner product.

        >>> x = clifford("{1}"); x &= clifford("{2}"); print x
        0
        >>> x = clifford(2); x &= clifford("{2}"); print x
        0
        >>> x = clifford("{1}"); x &= clifford("{1}"); print x
        1
        >>> x = clifford("{1}"); x &= clifford("{1,2}"); print x
        {2}
        """
        return self.wrap( self.unwrap() & toClifford(rhs) )

    def __xor__(lhs, rhs):
        """
        Outer product.

        >>> print clifford("{1}") ^ clifford("{2}")
        {1,2}
        >>> print clifford(2) ^ clifford("{2}")
        2{2}
        >>> print clifford("{1}") ^ clifford("{1}")
        0
        >>> print clifford("{1}") ^ clifford("{1,2}")
        0
        """
        return clifford().wrap( toClifford(lhs) ^ toClifford(rhs) )

    def __ixor__(self, rhs):
        """
        Outer product.

        >>> x = clifford("{1}"); x ^= clifford("{2}"); print x
        {1,2}
        >>> x = clifford(2); x ^= clifford("{2}"); print x
        2{2}
        >>> x = clifford("{1}"); x ^= clifford("{1}"); print x
        0
        >>> x = clifford("{1}"); x ^= clifford("{1,2}"); print x
        0
        """
        return self.wrap( self.unwrap() ^ toClifford(rhs) )

    def __div__(lhs, rhs):
        """
        Geometric quotient.

        >>> print clifford("{1}") / clifford("{2}")
        {1,2}
        >>> print clifford(2) / clifford("{2}")
        2{2}
        >>> print clifford("{1}") / clifford("{1}")
        1
        >>> print clifford("{1}") / clifford("{1,2}")
        -{2}
        """
        return clifford().wrap( toClifford(lhs) / toClifford(rhs) )

    def __idiv__(self, rhs):
        """
        Geometric quotient.

        >>> x = clifford("{1}"); x /= clifford("{2}"); print x
        {1,2}
        >>> x = clifford(2); x /= clifford("{2}"); print x
        2{2}
        >>> x = clifford("{1}"); x /= clifford("{1}"); print x
        1
        >>> x = clifford("{1}"); x /= clifford("{1,2}"); print x
        -{2}
        """
        return self.wrap( self.unwrap() / toClifford(rhs) )

    def inv(self):
        """
        Geometric multiplicative inverse.

        >>> x = clifford("{1}"); print inv(x)
        {1}
        >>> x = clifford(2); print inv(x)
        0.5
        >>> x = clifford("{1,2}"); print inv(x)
        -{1,2}
        """
        return clifford().wrap( self.instance.inv() )

    def __or__(lhs, rhs):
        """
        Transform left hand side, using right hand side as a transformation.

        >>> x=clifford("{1,2}")*pi/2; y=clifford("{1}"); print y|x
        -{1}
        >>> x=clifford("{1,2}")*pi/2; y=clifford("{1}"); print y|exp(x)
        -{1}
        """
        return clifford().wrap( toClifford(lhs) | toClifford(rhs) )

    def __ior__(self, rhs):
        """
        Transform left hand side, using right hand side as a transformation.

        >>> x=clifford("{1,2}")*pi/2; y=clifford("{1}"); y|=x; print y
        -{1}
        >>> x=clifford("{1,2}")*pi/2; y=clifford("{1}"); y|=exp(x); print y
        -{1}
        """
        return self.wrap( self.unwrap() | toClifford(rhs) )

    def __pow__(self, m, dummy):
        """
        Power: self to the m.

        >>> x=clifford("{1}"); print x ** 2
        1
        >>> x=clifford("2"); print x ** 2
        4
        >>> x=clifford("2+{1}"); print x ** 0
        1
        >>> x=clifford("2+{1}"); print x ** 1
        2+{1}
        >>> x=clifford("2+{1}"); print x ** 2
        5+4{1}
        """
        return pow(self, m)

    def pow(self, m):
        """
        Power: self to the m.

        >>> x=clifford("{1}"); print x.pow(2)
        1
        >>> x=clifford("2"); print x.pow(2)
        4
        >>> x=clifford("2+{1}"); print x.pow(0)
        1
        >>> x=clifford("2+{1}"); print x.pow(1)
        2+{1}
        >>> x=clifford("2+{1}"); print x.pow(2)
        5+4{1}
        >>> print clifford("1+{1}+{1,2}").pow(3)
        1+3{1}+3{1,2}
        """
        return clifford().wrap( self.instance.pow(m) )

    def outer_pow(self, m):
        """
        Outer product power.

        >>> x=clifford("2+{1}"); print x.outer_pow(0)
        1
        >>> x=clifford("2+{1}"); print x.outer_pow(1)
        2+{1}
        >>> x=clifford("2+{1}"); print x.outer_pow(2)
        4+4{1}
        >>> print clifford("1+{1}+{1,2}").outer_pow(3)
        1+3{1}+3{1,2}

        """
        return clifford().wrap( self.instance.outer_pow(m) )

    def __call__(self, grade):
        """
        Pure grade-vector part.

        >>> print clifford("{1}")(1)
        {1}
        >>> print clifford("{1}")(0)
        0
        >>> print clifford("1+{1}+{1,2}")(0)
        1
        >>> print clifford("1+{1}+{1,2}")(1)
        {1}
        >>> print clifford("1+{1}+{1,2}")(2)
        {1,2}
        >>> print clifford("1+{1}+{1,2}")(3)
        0
        """
        return clifford().wrap( self.instance.call(grade) )

    def scalar(self):
        """
        Scalar part.

        >>> clifford("1+{1}+{1,2}").scalar()
        1.0
        >>> clifford("{1,2}").scalar()
        0.0
        """
        return self.instance.scalar()

    def pure(self):
        """
        Pure part.

        >>> print clifford("1+{1}+{1,2}").pure()
        {1}+{1,2}
        >>> print clifford("{1,2}").pure()
        {1,2}
        """
        return clifford().wrap( self.instance.pure() )

    def even(self):
        """
        Even part of multivector, sum of even grade terms.

        >>> print clifford("1+{1}+{1,2}").even()
        1+{1,2}
        """
        return clifford().wrap( self.instance.even() )

    def odd(self):
        """
        Odd part of multivector, sum of odd grade terms.

        >>> print clifford("1+{1}+{1,2}").odd()
        {1}
        """
        return clifford().wrap( self.instance.odd() )

    def vector_part(self):
        """
        Vector part of multivector, as a Python list, with respect to frame.

        >>> print clifford("1+2{1}+3{2}+4{1,2}").vector_part()
        [2.0, 3.0]
        """
        cdef vector[scalar_t] vec = self.instance.vector_part()
        cdef int n = vec.size()
        lst = [0.0]*n
        cdef int i
        for i in xrange(n):
            lst[i] = vec[i]
        return lst

    def involute(self):
        """
        Main involution, each {i} is replaced by -{i} in each term,
        eg. clifford("{1}") -> -clifford("{1}").

        >>> print clifford("{1}").involute()
        -{1}
        >>> print (clifford("{2}")*clifford("{1}")).involute()
        -{1,2}
        >>> print (clifford("{1}")*clifford("{2}")).involute()
        {1,2}
        >>> print clifford("1+{1}+{1,2}").involute()
        1-{1}+{1,2}
        """
        return clifford().wrap( self.instance.involute() )

    def reverse(self):
        """
        Reversion, eg. clifford("{1}")*clifford("{2}") -> clifford("{2}")*clifford("{1}").

        >>> print clifford("{1}").reverse()
        {1}
        >>> print (clifford("{2}")*clifford("{1}")).reverse()
        {1,2}
        >>> print (clifford("{1}")*clifford("{2}")).reverse()
        -{1,2}
        >>> print clifford("1+{1}+{1,2}").reverse()
        1+{1}-{1,2}
        """
        return clifford().wrap( self.instance.reverse() )

    def conj(self):
        """
        Conjugation, reverse o involute == involute o reverse.

        >>> print (clifford("{1}")).conj()
        -{1}
        >>> print (clifford("{2}")*clifford("{1}")).conj()
        {1,2}
        >>> print (clifford("{1}")*clifford("{2}")).conj()
        -{1,2}
        >>> print clifford("1+{1}+{1,2}").conj()
        1-{1}-{1,2}
        """
        return clifford().wrap( self.instance.conj() )

    def quad(self):
        """
        Quadratic form == (rev(x)*x)(0).

        >>> print clifford("1+{1}+{1,2}").quad()
        3.0
        >>> print clifford("1+{-1}+{1,2}+{1,2,3}").quad()
        2.0
        """
        return self.instance.quad()

    def norm(self):
        """
        Norm == sum of squares of coordinates.

        >>> clifford("1+{1}+{1,2}").norm()
        3.0
        >>> clifford("1+{-1}+{1,2}+{1,2,3}").norm()
        4.0
        """
        return self.instance.norm()

    def abs(self):
        """
        Absolute value: square root of norm.

        >>> clifford("1+{-1}+{1,2}+{1,2,3}").abs()
        2.0
        """
        return glucat.abs( self.unwrap() )

    def max_abs(self):
        """
        Maximum of absolute values of components of multivector: multivector infinity norm.

        >>> clifford("1+{-1}+{1,2}+{1,2,3}").max_abs()
        1.0
        >>> clifford("3+2{1}+{1,2}").max_abs()
        3.0
        """
        return self.instance.max_abs()

    def truncated(self, limit):
        """
        Remove all terms of self with relative size smaller than limit.

        >>> clifford("1e8+{1}+1e-8{1,2}").truncated(1.0e-6)
        clifford("100000000")
        >>> clifford("1e4+{1}+1e-4{1,2}").truncated(1.0e-6)
        clifford("10000+{1}")
        """
        return clifford().wrap( self.instance.truncated(limit) )

    def isnan(self):
        """
        Check if a multivector contains any IEEE NaN values.

        >>> clifford().isnan()
        False
        """
        return self.instance.isnan()

    def frame(self):
        """
        Subalgebra generated by all generators of terms of given multivector.

        >>> print clifford("1+3{-1}+2{1,2}+4{-2,7}").frame()
        {-2,-1,1,2,7}
        >>> s=clifford("1+3{-1}+2{1,2}+4{-2,7}").frame(); type(s)
        <type 'PyClical.index_set'>
        """
        return index_set().wrap( self.instance.frame() )

    def __repr__(self):
        """
        The “official” string representation of self.

        >>> clifford("1+3{-1}+2{1,2}+4{-2,7}").__repr__()
        'clifford("1+3{-1}+2{1,2}+4{-2,7}")'
        """
        return clifford_to_repr( self.unwrap() ).c_str()

    def __str__(self):
        """
        The “informal” string representation of self.

        >>> clifford("1+3{-1}+2{1,2}+4{-2,7}").__str__()
        '1+3{-1}+2{1,2}+4{-2,7}'
        """
        return clifford_to_str( self.unwrap() ).c_str()

def clifford_hidden_doctests():
    """
    Tests for functions that Doctest cannot see.

    For clifford.__cinit__: Construct an object of type clifford.

    >>> print clifford(2)
    2
    >>> print clifford(2L)
    2
    >>> print clifford(2.0)
    2
    >>> print clifford(1.0e-1)
    0.1
    >>> print clifford("2")
    2
    >>> print clifford("2{1,2,3}")
    2{1,2,3}
    >>> print clifford(clifford("2{1,2,3}"))
    2{1,2,3}
    >>> print clifford("-{1}")
    -{1}
    >>> print clifford(2,index_set("{1,2}"))
    2{1,2}
    >>> print clifford([2,3],index_set("{1,2}"))
    2{1}+3{2}
    >>> print clifford([1,2])
    Traceback (most recent call last):
      ...
    TypeError: Cannot initialize clifford object from <type 'list'>.
    >>> print clifford(None)
    Traceback (most recent call last):
      ...
    TypeError: Cannot initialize clifford object from <type 'NoneType'>.
    >>> print clifford(None,[1,2])
    Traceback (most recent call last):
      ...
    TypeError: Cannot initialize clifford object from (<type 'NoneType'>, <type 'list'>).
    >>> print clifford([1,2],[1,2])
    Traceback (most recent call last):
      ...
    TypeError: Cannot initialize clifford object from (<type 'list'>, <type 'list'>).
    >>> print clifford("")
    Traceback (most recent call last):
      ...
    ValueError: Cannot initialize clifford object from invalid string ''.
    >>> print clifford("{")
    Traceback (most recent call last):
      ...
    ValueError: Cannot initialize clifford object from invalid string '{'.
    >>> print clifford("{1")
    Traceback (most recent call last):
      ...
    ValueError: Cannot initialize clifford object from invalid string '{1'.
    >>> print clifford("+")
    Traceback (most recent call last):
      ...
    ValueError: Cannot initialize clifford object from invalid string '+'.
    >>> print clifford("-")
    Traceback (most recent call last):
      ...
    ValueError: Cannot initialize clifford object from invalid string '-'.
    >>> print clifford("{1}+")
    Traceback (most recent call last):
      ...
    ValueError: Cannot initialize clifford object from invalid string '{1}+'.

    For clifford.__richcmp__: Compare objects of type clifford.

    >>> clifford("{1}") == clifford("1{1}")
    True
    >>> clifford("{1}") != clifford("1.0{1}")
    False
    >>> clifford("{1}") != clifford("1.0")
    True
    >>> clifford("{1,2}") == None
    False
    >>> clifford("{1,2}") != None
    True
    >>> None == clifford("{1,2}")
    False
    >>> None != clifford("{1,2}")
    True
    """
    return

cpdef inline inv(obj):
    """
    Geometric multiplicative inverse.

    >>> print inv(clifford("{1}"))
    {1}
    >>> print inv(clifford("{-1}"))
    -{-1}
    >>> print inv(clifford("{-2,-1}"))
    -{-2,-1}
    >>> print inv(clifford("{-1}+{1}"))
    nan
    """
    return clifford(obj).inv()

cpdef inline scalar(obj):
    """
    Scalar part.

    >>> scalar(clifford("1+{1}+{1,2}"))
    1.0
    >>> scalar(clifford("{1,2}"))
    0.0
    """
    return clifford(obj).scalar()

cpdef inline real(obj):
    """
    Real part: synonym for scalar part.

    >>> real(clifford("1+{1}+{1,2}"))
    1.0
    >>> real(clifford("{1,2}"))
    0.0
    """
    return clifford(obj).scalar()

cpdef inline imag(obj):
    """
    Imaginary part: deprecated (always 0).

    >>> imag(clifford("1+{1}+{1,2}"))
    0.0
    >>> imag(clifford("{1,2}"))
    0.0
    """
    return 0.0

cpdef inline pure(obj):
    """
    Pure part

    >>> print pure(clifford("1+{1}+{1,2}"))
    {1}+{1,2}
    >>> print pure(clifford("{1,2}"))
    {1,2}
    """
    return clifford(obj).pure()

cpdef inline even(obj):
    """
    Even part of multivector, sum of even grade terms.

    >>> print even(clifford("1+{1}+{1,2}"))
    1+{1,2}
    """
    return clifford(obj).even()

cpdef inline odd(obj):
    """
    Odd part of multivector, sum of odd grade terms.

    >>> print odd(clifford("1+{1}+{1,2}"))
    {1}
    """
    return clifford(obj).odd()

cpdef inline involute(obj):
    """
    Main involution, each {i} is replaced by -{i} in each term, eg. {1}*{2} -> (-{2})*(-{1})

    >>> print involute(clifford("{1}"))
    -{1}
    >>> print involute(clifford("{2}")*clifford("{1}"))
    -{1,2}
    >>> print involute(clifford("{1}")*clifford("{2}"))
    {1,2}
    >>> print involute(clifford("1+{1}+{1,2}"))
    1-{1}+{1,2}
    """
    return clifford(obj).involute()

cpdef inline reverse(obj):
    """
    Reversion, eg. {1}*{2} -> {2}*{1}

    >>> print reverse(clifford("{1}"))
    {1}
    >>> print reverse(clifford("{2}")*clifford("{1}"))
    {1,2}
    >>> print reverse(clifford("{1}")*clifford("{2}"))
    -{1,2}
    >>> print reverse(clifford("1+{1}+{1,2}"))
    1+{1}-{1,2}
    """
    return clifford(obj).reverse()

cpdef inline conj(obj):
    """
    Conjugation, reverse o involute == involute o reverse.

    >>> print conj(clifford("{1}"))
    -{1}
    >>> print conj(clifford("{2}")*clifford("{1}"))
    {1,2}
    >>> print conj(clifford("{1}")*clifford("{2}"))
    -{1,2}
    >>> print conj(clifford("1+{1}+{1,2}"))
    1-{1}-{1,2}
    """
    return clifford(obj).conj()

cpdef inline quad(obj):
    """
    Quadratic form == (rev(x)*x)(0).

    >>> print quad(clifford("1+{1}+{1,2}"))
    3.0
    >>> print quad(clifford("1+{-1}+{1,2}+{1,2,3}"))
    2.0
    """
    return clifford(obj).quad()

cpdef inline norm(obj):
    """
    norm == sum of squares of coordinates.

    >>> norm(clifford("1+{1}+{1,2}"))
    3.0
    >>> norm(clifford("1+{-1}+{1,2}+{1,2,3}"))
    4.0
    """
    return clifford(obj).norm()

cpdef inline abs(obj):
    """
    Absolute value of multivector: multivector 2-norm.

    >>> abs(clifford("1+{-1}+{1,2}+{1,2,3}"))
    2.0
    """
    return glucat.abs(toClifford(obj))

cpdef inline max_abs(obj):
    """
    Maximum absolute value of coordinates multivector: multivector infinity-norm.

    >>> max_abs(clifford("1+{-1}+{1,2}+{1,2,3}"))
    1.0
    >>> max_abs(clifford("3+2{1}+{1,2}"))
    3.0

    """
    return glucat.max_abs(toClifford(obj))

cpdef inline pow(obj, m):
    """
    Integer power of multivector: obj to the m.

    >>> x=clifford("{1}"); print pow(x,2)
    1
    >>> x=clifford("2"); print pow(x,2)
    4
    >>> x=clifford("2+{1}"); print pow(x,0)
    1
    >>> x=clifford("2+{1}"); print pow(x,1)
    2+{1}
    >>> x=clifford("2+{1}"); print pow(x,2)
    5+4{1}
    >>> print pow(clifford("1+{1}+{1,2}"),3)
    1+3{1}+3{1,2}
    """
    try:
        math.pow(obj, m)
    except:    
        return clifford(obj).pow(m)

cpdef inline outer_pow(obj, m):
    """
    Outer product power of multivector.

    >>> print outer_pow(clifford("1+{1}+{1,2}"),3)
    1+3{1}+3{1,2}
    """
    return clifford(obj).outer_pow(m)

cpdef inline complexifier(obj):
    """
    Square root of -1 which commutes with all members of the frame of the given multivector.

    >>> print complexifier(clifford(index_set("{1}")))
    {1,2,3}
    >>> print complexifier(clifford(index_set("{-1}")))
    {-1}
    >>> print complexifier(index_set("{1}"))
    {1,2,3}
    >>> print complexifier(index_set("{-1}"))
    {-1}
    """
    return clifford().wrap( glucat.complexifier(toClifford(obj)) )

cpdef inline sqrt(obj, i = None):
    """
    Square root of multivector with optional complexifier.

    >>> print sqrt(-1)
    {-1}
    >>> print sqrt(clifford("2{-1}"))
    1+{-1}
    >>> j=sqrt(-1,complexifier(index_set("{1}"))); print j; print j*j
    {1,2,3}
    -1
    >>> j=sqrt(-1,"{1,2,3}"); print j; print j*j
    {1,2,3}
    -1
    """
    if not (i is None):
        return clifford().wrap( glucat.sqrt(toClifford(obj), toClifford(i)) )
    else:
        try:
            return math.sqrt(obj)
        except:
            return clifford().wrap( glucat.sqrt(toClifford(obj)) )

cpdef inline exp(obj):
    """
    Exponential of multivector.

    >>> x=clifford("{1,2}")*pi/4; print exp(x)
    0.7071+0.7071{1,2}
    >>> x=clifford("{1,2}")*pi/2; print exp(x)
    {1,2}
    """
    try:
        return math.exp(obj)
    except:
        return clifford().wrap( glucat.exp(toClifford(obj)) )

cpdef inline log(obj,i = None):
    """
    Natural logarithm of multivector with optional complexifier.

    >>> x=clifford("{1,2}"); print (log(x,"{1,2}")*2/pi)
    {1,2}
    >>> x=clifford("{1,2}"); print (log(x,"{1,2,3}")*2/pi)
    {1,2}
    >>> x=clifford("{1,2}"); print (log(x)*2/pi)
    {1,2}
    """
    if not (i is None):
        return clifford().wrap( glucat.log(toClifford(obj), toClifford(i)) )
    else:
        try:
            return math.log(obj)
        except:
            return clifford().wrap( glucat.log(toClifford(obj)) )

cpdef inline cos(obj,i = None):
    """
    Cosine of multivector with optional complexifier.

    >>> x=clifford("{1,2}"); print cos(acos(x),"{1,2,3}")
    {1,2}
    >>> x=clifford("{1,2}"); print cos(acos(x))
    {1,2}
    """
    if not (i is None):
        return clifford().wrap( glucat.cos(toClifford(obj), toClifford(i)) )
    else:
        try:
            return math.cos(obj)
        except:
            return clifford().wrap( glucat.cos(toClifford(obj)) )

cpdef inline acos(obj,i = None):
    """
    Inverse cosine of multivector with optional complexifier.

    >>> x=clifford("{1,2}"); print cos(acos(x),"{1,2,3}")
    {1,2}
    >>> x=clifford("{1,2}"); print cos(acos(x),"{-1,1,2,3,4}")
    {1,2}
    >>> print acos(0) / pi
    0.5
    >>> x=clifford("{1,2}"); print cos(acos(x))
    {1,2}
    """
    if not (i is None):
        return clifford().wrap( glucat.acos(toClifford(obj), toClifford(i)) )
    else:
        try:
            return math.acos(obj)
        except:
            return clifford().wrap( glucat.acos(toClifford(obj)) )

cpdef inline cosh(obj):
    """
    Hyperbolic cosine of multivector.

    >>> x=clifford("{1,2}")*pi; print cosh(x)
    -1
    >>> x=clifford("{1,2,3}"); print cosh(acosh(x))
    {1,2,3}
    >>> x=clifford("{1,2}"); print cosh(acosh(x))
    {1,2}
    """
    try:
        return math.cosh(obj)
    except:
        return clifford().wrap( glucat.cosh(toClifford(obj)) )

cpdef inline acosh(obj,i = None):
    """
    Inverse hyperbolic cosine of multivector with optional complexifier.

    >>> print acosh(0,"{-2,-1,1}")
    1.571{-2,-1,1}
    >>> x=clifford("{1,2,3}"); print cosh(acosh(x,"{-1,1,2,3,4}"))
    {1,2,3}
    >>> print acosh(0)
    1.571{-1}
    >>> x=clifford("{1,2,3}"); print cosh(acosh(x))
    {1,2,3}
    >>> x=clifford("{1,2}"); print cosh(acosh(x))
    {1,2}
    """
    if not (i is None):
        return clifford().wrap( glucat.acosh(toClifford(obj), toClifford(i)) )
    else:
        try:
            return math.acosh(obj)
        except:
            return clifford().wrap( glucat.acosh(toClifford(obj)) )

cpdef inline sin(obj,i = None):
    """
    Sine of multivector with optional complexifier.

    >>> s="{-1}"; x=clifford(s); print asin(sin(x,s),s)
    {-1}
    >>> s="{-1}"; x=clifford(s); print asin(sin(x,s),"{-2,-1,1}")
    {-1}
    >>> x=clifford("{1,2,3}"); print asin(sin(x))
    {1,2,3}
    """
    if not (i is None):
        return clifford().wrap( glucat.sin(toClifford(obj), toClifford(i)) )
    else:
        try:
            return math.sin(obj)
        except:
            return clifford().wrap( glucat.sin(toClifford(obj)) )

cpdef inline asin(obj,i = None):
    """
    Inverse sine of multivector with optional complexifier.

    >>> s="{-1}"; x=clifford(s); print asin(sin(x,s),s)
    {-1}
    >>> s="{-1}"; x=clifford(s); print asin(sin(x,s),"{-2,-1,1}")
    {-1}
    >>> print asin(1) / pi
    0.5
    >>> x=clifford("{1,2,3}"); print asin(sin(x))
    {1,2,3}
    """
    if not (i is None):
        return clifford().wrap( glucat.asin(toClifford(obj), toClifford(i)) )
    else:
        try:
            return math.asin(obj)
        except:
            return clifford().wrap( glucat.asin(toClifford(obj)) )

cpdef inline sinh(obj):
    """
    Hyperbolic sine of multivector.

    >>> x=clifford("{1,2}")*pi/2; print sinh(x)
    {1,2}
    >>> x=clifford("{1,2}")*pi/6; print sinh(x)
    0.5{1,2}
    """
    try:
        return math.sinh(obj)
    except:
        return clifford().wrap( glucat.sinh(toClifford(obj)) )

cpdef inline asinh(obj,i = None):
    """
    Inverse hyperbolic sine of multivector with optional complexifier.

    >>> x=clifford("{1,2}"); print asinh(x,"{1,2,3}")*2/pi
    {1,2}
    >>> x=clifford("{1,2}"); print asinh(x)*2/pi
    {1,2}
    >>> x=clifford("{1,2}") / 2; print asinh(x)*6/pi
    {1,2}
    """
    if not (i is None):
        return clifford().wrap( glucat.asinh(toClifford(obj), toClifford(i)) )
    else:
        try:
            return math.asinh(obj)
        except:
            return clifford().wrap( glucat.asinh(toClifford(obj)) )

cpdef inline tan(obj,i = None):
    """
    Tangent of multivector with optional complexifier.

    >>> x=clifford("{1,2}"); print tan(x,"{1,2,3}")
    0.7616{1,2}
    >>> x=clifford("{1,2}"); print tan(x)
    0.7616{1,2}
    """
    if not (i is None):
        return clifford().wrap( glucat.tan(toClifford(obj), toClifford(i)) )
    else:
        try:
            return math.tan(obj)
        except:
            return clifford().wrap( glucat.tan(toClifford(obj)) )

cpdef inline atan(obj,i = None):
    """
    Inverse tangent of multivector with optional complexifier.

    >>> s=index_set("{1,2,3}"); x=clifford("{1}"); print tan(atan(x,s),s)
    {1}
    >>> x=clifford("{1}"); print tan(atan(x))
    {1}
    """
    if not (i is None):
        return clifford().wrap( glucat.atan(toClifford(obj), toClifford(i)) )
    else:
        try:
            return math.atan(obj)
        except:
            return clifford().wrap( glucat.atan(toClifford(obj)) )

cpdef inline tanh(obj):
    """
    Hyperbolic tangent of multivector.

    >>> x=clifford("{1,2}")*pi/4; print tanh(x)
    {1,2}
    """
    try:
        return math.tanh(obj)
    except:
        return clifford().wrap( glucat.tanh(toClifford(obj)) )

cpdef inline atanh(obj,i = None):
    """
    Inverse hyperbolic tangent of multivector with optional complexifier.

    >>> s=index_set("{1,2,3}"); x=clifford("{1,2}"); print tanh(atanh(x,s))
    {1,2}
    >>> x=clifford("{1,2}"); print tanh(atanh(x))
    {1,2}
    """
    if not (i is None):
        return clifford().wrap( glucat.atanh(toClifford(obj), toClifford(i)) )
    else:
        try:
            return math.atanh(obj)
        except:
            return clifford().wrap( glucat.atanh(toClifford(obj)) )

cpdef inline random_clifford(index_set ist, fill = 1.0):
    """
    Random multivector within a frame.

    >>> print random_clifford(index_set("{-3,-1,2}")).frame()
    {-3,-1,2}
    """
    return clifford().wrap( clifford().instance.random(ist.unwrap(), <scalar_t>fill) )

cpdef inline cga3(obj):
    """
    Convert Euclidean 3D multivector to Conformal Geometric Algebra using Doran and Lasenby definition.

    >>> x=clifford("2{1}+9{2}+{3}"); print cga3(x)
    87{-1}+4{1}+18{2}+2{3}+85{4}
    """
    return clifford().wrap( glucat.cga3(toClifford(obj)) )

cpdef inline cga3std(obj):
    """
    Convert CGA3 null vector to standard conformal null vector using Doran and Lasenby definition.

    >>> x=clifford("2{1}+9{2}+{3}"); print cga3std(cga3(x))
    87{-1}+4{1}+18{2}+2{3}+85{4}
    >>> x=clifford("2{1}+9{2}+{3}"); print cga3std(cga3(x))-cga3(x)
    0
    """
    return clifford().wrap( glucat.cga3std(toClifford(obj)) )

cpdef inline agc3(obj):
    """
    Convert CGA3 null vector to Euclidean 3D vector using Doran and Lasenby definition.

    >>> x=clifford("2{1}+9{2}+{3}"); print agc3(cga3(x))
    2{1}+9{2}+{3}
    >>> x=clifford("2{1}+9{2}+{3}"); print agc3(cga3(x))-x
    0
    """
    return clifford().wrap( glucat.agc3(toClifford(obj)) )
    
# Some abbreviations.
pi = atan(clifford(1.0)) * 4.0

def cl(obj):
    """
    Abbreviation for clifford(obj).

    >>> print cl(2)
    2
    >>> print cl(2L)
    2
    >>> print cl(2.0)
    2
    >>> print cl(5.0e-1)
    0.5
    >>> print cl("2")
    2
    >>> print cl("2{1,2,3}")
    2{1,2,3}
    >>> print cl(cl("2{1,2,3}"))
    2{1,2,3}
    """
    return clifford(obj)

def ist(obj):
    """
    Abbreviation for index_set(obj).

    >>> print ist("{1,2,3}")
    {1,2,3}
    """
    return index_set(obj)

def e(obj):
    """
    Abbreviation for clifford(index_set(obj)).

    >>> print e(1)
    {1}
    >>> print e(-1)
    {-1}
    >>> print e(0)
    1
    """
    return clifford(index_set(obj))

def istpq(p, q):
    """
    Abbreviation for index_set("{-q,...p}").

    >>> print istpq(2,3)
    {-3,-2,-1,1,2}
    """
    return index_set(str(range(-q,p+1)).replace('[','{').replace(']','}'))

ninf3 = e(4) + e(-1) # Null infinity point in 3D Conformal Geometric Algebra.
nbar3 = e(4) - e(-1) # Null bar point in 3D Conformal Geometric Algebra.

# Doctest interface.
def _test():
    import PyClical, doctest
    return doctest.testmod(PyClical)

if __name__ == "__main__":
    _test()
