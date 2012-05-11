# -*- coding: utf-8 -*-
#
# PyClical: Python interface to GluCat:
#           Generic library of universal Clifford algebra templates
#
# glucat.pxd: Basic Cython definitions
#             corresponding to C++ definitions from PyClical.h.
# Kept as a separate module from PyClical.pxd to avoid namespace clashes.
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

from libcpp.vector cimport vector

cdef extern from "PyClical.h":

    cdef cppclass String:
        char* c_str()

    cdef cppclass IndexSet:
        IndexSet ()
        IndexSet (IndexSet Ist) except+
        IndexSet (int idx) except+
        IndexSet (char* str) except+
        inline bint operator==(IndexSet Rhs)
        inline bint operator!=(IndexSet Rhs)
        inline bint operator<(IndexSet Rhs)
        inline IndexSet invert "operator~"()
        inline bint getitem "operator[]"(int idx)
        inline IndexSet set()
        inline IndexSet set(int idx) except+
        inline IndexSet set(int idx, int val) except+
        inline IndexSet reset()
        inline IndexSet reset(int idx) except+
        int count()
        int count_pos()
        int count_neg()
        int min()
        int max()
        int sign_of_mult(IndexSet Rhs)
        int sign_of_square()
        int hash_fn()

    int compare(IndexSet Lhs, IndexSet Rhs)
    int min_neg(IndexSet Ist)
    int max_pos(IndexSet Ist)

    ctypedef double scalar_t

    cdef cppclass Clifford:
        Clifford ()
        Clifford (Clifford Clf) except+
        Clifford (scalar_t scr) except+
        Clifford (char* str) except+
        Clifford (IndexSet ist, scalar_t scr) except+
        Clifford (vector[scalar_t] vec, IndexSet ist) except+
        bint     operator==(Clifford Rhs)
        bint     operator!=(Clifford Rhs)
        Clifford neg "operator-"()
        scalar_t getitem "operator[]"(IndexSet Ist)
        Clifford call "operator()"(int grade)
        scalar_t scalar()
        Clifford pure()
        Clifford even()
        Clifford odd()
        vector[scalar_t] vector_part()
        Clifford involute()
        Clifford reverse()
        Clifford conj()
        Clifford random(IndexSet Ist, scalar_t fill)
        scalar_t norm()
        scalar_t quad()
        IndexSet frame()
        scalar_t max_abs()
        Clifford inv()
        Clifford pow(int m)
        Clifford outer_pow(int m)
        Clifford truncated(scalar_t limit)
        bint     isnan()
        void     write(char* msg)

    scalar_t scalar(Clifford Clf)
    scalar_t real(Clifford Clf)
    scalar_t imag(Clifford Clf)
    Clifford pure(Clifford Clf)
    Clifford even(Clifford Clf)
    Clifford odd(Clifford Clf)
    Clifford involute(Clifford Clf)
    Clifford reverse(Clifford Clf)
    Clifford conj(Clifford Clf)
    scalar_t norm(Clifford Clf)
    scalar_t abs(Clifford Clf)
    scalar_t max_abs(Clifford Clf)
    scalar_t quad(Clifford Clf)
    Clifford inv(Clifford Clf)
    Clifford pow(Clifford Clf,int m)
    Clifford outer_pow(Clifford Clf,int m)

    Clifford complexifier(Clifford Clf)
    Clifford sqrt(Clifford Clf, Clifford I) except+
    Clifford sqrt(Clifford Clf)
    Clifford exp(Clifford Clf)
    Clifford log(Clifford Clf, Clifford I) except+
    Clifford log(Clifford Clf)
    Clifford cos(Clifford Clf, Clifford I) except+
    Clifford cos(Clifford Clf)
    Clifford acos(Clifford Clf, Clifford I) except+
    Clifford acos(Clifford Clf)
    Clifford cosh(Clifford Clf)
    Clifford acosh(Clifford Clf, Clifford I) except+
    Clifford acosh(Clifford Clf)
    Clifford sin(Clifford Clf, Clifford I) except+
    Clifford sin(Clifford Clf)
    Clifford asin(Clifford Clf, Clifford I) except+
    Clifford asin(Clifford Clf)
    Clifford sinh(Clifford Clf)
    Clifford asinh(Clifford Clf, Clifford I) except+
    Clifford asinh(Clifford Clf)
    Clifford tan(Clifford Clf, Clifford I) except+
    Clifford tan(Clifford Clf)
    Clifford atan(Clifford Clf, Clifford I) except+
    Clifford atan(Clifford Clf)
    Clifford tanh(Clifford Clf)
    Clifford atanh(Clifford Clf, Clifford I) except+
    Clifford atanh(Clifford Clf)

cdef extern from "PyClical.h" namespace "cga3":
    Clifford agc3(Clifford Clf)
    Clifford cga3(Clifford Clf)
    Clifford cga3std(Clifford Clf)
