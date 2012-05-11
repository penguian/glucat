# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# PyClical.pxd: Basic Cython definitions for PyClical
#               corresponding to C++ definitions from PyClical.h.
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

cimport glucat
from glucat cimport IndexSet, String, Clifford, scalar_t, vector

cdef extern from "PyClical.h":
    inline IndexSet operator&(IndexSet Lhs, IndexSet Rhs)
    inline IndexSet operator|(IndexSet Lhs, IndexSet Rhs)
    inline IndexSet operator^(IndexSet Lhs, IndexSet Rhs)

    inline String index_set_to_repr(IndexSet& Ist)
    inline String index_set_to_str(IndexSet& Ist)

    inline Clifford operator+(Clifford Lhs, Clifford Rhs)
    inline Clifford operator-(Clifford Lhs, Clifford Rhs)
    inline Clifford operator*(Clifford Lhs, Clifford Rhs)
    inline Clifford operator&(Clifford Lhs, Clifford Rhs)
    inline Clifford operator%(Clifford Lhs, Clifford Rhs)
    inline Clifford operator^(Clifford Lhs, Clifford Rhs)
    inline Clifford operator/(Clifford Lhs, Clifford Rhs)
    inline Clifford operator|(Clifford Lhs, Clifford Rhs)

    inline String clifford_to_repr(Clifford& Clf)
    inline String clifford_to_str(Clifford& Clf)
