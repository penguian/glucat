# -*- coding: utf-8 -*-
# cython: language_level=3
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
from libcpp.string cimport string

cdef extern from "PyClical.h":
    IndexSet operator&(IndexSet Lhs, IndexSet Rhs)
    IndexSet operator|(IndexSet Lhs, IndexSet Rhs)
    IndexSet operator^(IndexSet Lhs, IndexSet Rhs)

    string index_set_to_repr(IndexSet& Ist)
    string index_set_to_str(IndexSet& Ist)

    Clifford operator+(Clifford Lhs, Clifford Rhs)
    Clifford operator-(Clifford Lhs, Clifford Rhs)
    Clifford operator*(Clifford Lhs, Clifford Rhs)
    Clifford operator&(Clifford Lhs, Clifford Rhs)
    Clifford operator%(Clifford Lhs, Clifford Rhs)
    Clifford operator^(Clifford Lhs, Clifford Rhs)
    Clifford operator/(Clifford Lhs, Clifford Rhs)
    Clifford operator|(Clifford Lhs, Clifford Rhs)

    string clifford_to_repr(Clifford& Clf)
    string clifford_to_str(Clifford& Clf)
