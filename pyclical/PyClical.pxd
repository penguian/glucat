# -*- coding: utf-8 -*-
# cython: language_level=3
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# PyClical.pxd: Basic Cython definitions for PyClical
#               corresponding to C++ definitions from PyClical.h.
#
#    copyright            : (C) 2008-2021 by Paul C. Leopardi
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
    string glucat_package_version

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

    const scalar_t epsilon

cdef extern from *:
    """
    inline IndexSet* new_IndexSet() { return new IndexSet(); }
    inline IndexSet* new_IndexSet_copy(IndexSet other) { return new IndexSet(other); }
    inline IndexSet* new_IndexSet_int(int val) { return new IndexSet(val); }
    inline IndexSet* new_IndexSet_str(const char* val) { return new IndexSet(val); }
    inline void delete_IndexSet(IndexSet* ptr) { delete ptr; }

    inline Clifford* new_Clifford() { return new Clifford(); }
    inline Clifford* new_Clifford_copy(Clifford other) { return new Clifford(other); }
    inline Clifford* new_Clifford_frame_coeff(Clifford other, IndexSet ist) { return new Clifford(other, ist); }
    inline Clifford* new_Clifford_scalar(scalar_t scr) { return new Clifford(scr); }
    inline Clifford* new_Clifford_str(const char* str) { return new Clifford(str); }
    inline Clifford* new_Clifford_frame_scalar(IndexSet ist, scalar_t scr) { return new Clifford(ist, scr); }
    inline Clifford* new_Clifford_vec_frame(std::vector<scalar_t> vec, IndexSet ist) { return new Clifford(vec, ist); }
    inline void delete_Clifford(Clifford* ptr) { delete ptr; }
    """
    IndexSet* new_IndexSet() except +
    IndexSet* new_IndexSet_copy(IndexSet other) except +
    IndexSet* new_IndexSet_int(int val) except +
    IndexSet* new_IndexSet_str(const char* val) except +
    void delete_IndexSet(IndexSet* ptr) except +

    Clifford* new_Clifford() except +
    Clifford* new_Clifford_copy(Clifford other) except +
    Clifford* new_Clifford_frame_coeff(Clifford other, IndexSet ist) except +
    Clifford* new_Clifford_scalar(scalar_t scr) except +
    Clifford* new_Clifford_str(const char* str) except +
    Clifford* new_Clifford_frame_scalar(IndexSet ist, scalar_t scr) except +
    Clifford* new_Clifford_vec_frame(vector[scalar_t] vec, IndexSet ist) except +
    void delete_Clifford(Clifford* ptr) except +
