#ifndef _GLUCAT_PORTABILITY_H
#define _GLUCAT_PORTABILITY_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    portability.h : Work around non-standard compilers and libraries
                             -------------------
    begin                : Sun 2001-08-18
    copyright            : (C) 2001-2016 by Paul C. Leopardi
 ***************************************************************************

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library.  If not, see <http://www.gnu.org/licenses/>.

 ***************************************************************************
 This library is based on a prototype written by Arvind Raja and was
 licensed under the LGPL with permission of the author. See Arvind Raja,
 "Object-oriented implementations of Clifford algebras in C++: a prototype",
 in Ablamowicz, Lounesto and Parra (eds.)
 "Clifford algebras with numeric and symbolic computations", Birkhauser, 1996.
 ***************************************************************************
     See also Arvind Raja's original header comments in glucat.h
 ***************************************************************************/

#include <boost/version.hpp>
#include <cmath>

// Workaround for GCC
#if (__GNUC__ == 4 && __GNUC_MINOR__ >= 8)
# define _GLUCAT_GCC_IGNORE_UNUSED_LOCAL_TYPEDEFS
#endif

// Workaround for isnan and isinf
#if __cplusplus > 199711L
# define _GLUCAT_ISNAN(x) (std::isnan(x))
# define _GLUCAT_ISINF(x) (std::isinf(x))
#else
# define _GLUCAT_ISNAN(x) (x != x)
# define _GLUCAT_ISINF(x) (!_GLUCAT_ISNAN(x) && _GLUCAT_ISNAN(x-x))
#endif

// Workaround for abs and sqrt
#if BOOST_VERSION >= 103400
# define UBLAS_ABS  type_abs
# define UBLAS_SQRT type_sqrt
#else
# define UBLAS_ABS  abs
# define UBLAS_SQRT sqrt
#endif

// Use with Cygwin gcc to obtain __WORDSIZE
#if defined(HAVE_BITS_WORDSIZE_H)
# include <bits/wordsize.h>
#endif

#endif // _GLUCAT_PORTABILITY_H
