#ifndef _GLUCAT_LONG_DOUBLE_H
#define _GLUCAT_LONG_DOUBLE_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    long_double.h : Define std functions for long double
                             -------------------
    begin                : 2001-12-18
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
 "Clifford algebras with numeric and symbolic computations, Birkhauser, 1996."
 ***************************************************************************
 See also Arvind Raja's original header comments and references in glucat.h
 ***************************************************************************/

#include "glucat/global.h"
#include "glucat/scalar.h"

namespace glucat
{
#if defined(__USE_GNU)
  static const long double l_pi   = M_PIl;
  static const long double l_ln2 = M_LN2l;
#else
  static const long double l_pi = 3.1415926535897932384626433832795029L;
  static const long double l_ln2 = 0.6931471805599453094172321214581766L;
#endif

  /// Pi for long double
  template<>
  inline
  auto
  numeric_traits<long double>::
  pi() -> long double
  { return l_pi; }

  /// log(2) for long double
  template<>
  inline
  auto
  numeric_traits<long double>::
  ln_2() -> long double
  { return l_ln2; }
}
#endif // _GLUCAT_LONG_DOUBLE_H
