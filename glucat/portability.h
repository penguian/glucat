#ifndef _GLUCAT_PORTABILITY_H
#define _GLUCAT_PORTABILITY_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    portability.h : Work around non-standard compilers and libraries
                             -------------------
    begin                : Sun 2001-08-18
    copyright            : (C) 2001-2012 by Paul C. Leopardi
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

// Workarounds for ICC and ICPC
#if defined (BOOST_INTEL) || defined (__INTEL_COMPILER) || defined (__ICL) || defined (__ICC)
# pragma warning( disable: 177 ) // variable was declared but never referenced
# pragma warning( disable: 193 ) // zero used for undefined preprocessing identifier
# pragma warning( disable: 279 ) // controlling expression is constant
# pragma warning( disable: 383 ) // value copied to temporary, reference to temporary ...
# pragma warning( disable: 444 ) // destructor for base is not virtual
# pragma warning( disable: 593 ) // variable was set but never used
# pragma warning( disable: 810 ) // conversion from "double" to "int" may lose significant bits
# pragma warning( disable: 858 ) // type qualifier on return type is meaningless
# pragma warning( disable: 869 ) // parameter was never referenced
# pragma warning( disable: 981 ) // operands are evaluated in unspecified order
# pragma warning( disable: 1572 ) // floating-point equality and inequality comparisons ...
# pragma warning( disable: 2259 ) // non-pointer conversion from "double" to "...={float}" may lose significant bits
#endif

// ICPC does not have std::tr1::isnan() or std:tr1::isinf()
#if defined (BOOST_INTEL) || defined (__INTEL_COMPILER) || defined (__ICL) || defined (__ICC)
# define _GLUCAT_ISNAN(x) (x != x)
# define _GLUCAT_ISINF(x) (!_GLUCAT_ISNAN(x) && _GLUCAT_ISNAN(x-x))
#else
# define _GLUCAT_ISNAN(x) (std::isnan(x))
# define _GLUCAT_ISINF(x) (std::isinf(x))
#endif

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
