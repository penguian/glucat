#ifndef _GLUCAT_SCALAR_IMP_H
#define _GLUCAT_SCALAR_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    scalar_imp.h : Define functions for scalar_t
                             -------------------
    begin                : 2001-12-20
    copyright            : (C) 2001-2014 by Paul C. Leopardi
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

#include "glucat/portability.h"

#include <boost/numeric/ublas/traits.hpp>

#include <cmath>
#include <limits>

namespace glucat
{
  /// Extra traits which extend numeric limits
  // Reference: [AA], 2.4, p. 30-31

  /// Cast to float
  template< >
  template< typename Other_Scalar_T >
  inline
  float
  numeric_traits<float>::
  to_scalar_t(const Other_Scalar_T& val)
  { return float(numeric_traits<Other_Scalar_T>::to_double(val)); }

  /// Cast to double
  template< >
  template< typename Other_Scalar_T >
  inline
  double
  numeric_traits<double>::
  to_scalar_t(const Other_Scalar_T& val)
  { return numeric_traits<Other_Scalar_T>::to_double(val); }

#if defined(_GLUCAT_USE_QD)
  /// Cast to long double
  template< >
  template< >
  inline
  long double
  numeric_traits<long double>::
  to_scalar_t(const dd_real& val)
  { return static_cast<long double>(val.x[0]) + static_cast<long double>(val.x[1]); }

  /// Cast to long double
  template< >
  template< >
  inline
  long double
  numeric_traits<long double>::
  to_scalar_t(const qd_real& val)
  { return static_cast<long double>(val.x[0]) + static_cast<long double>(val.x[1]); }

  /// Cast to dd_real
  template< >
  template< >
  inline
  dd_real
  numeric_traits<dd_real>::
  to_scalar_t(const long double& val)
  { return dd_real(double(val),double(val - static_cast<long double>(double(val)))); }

  /// Cast to dd_real
  template< >
  template< >
  inline
  dd_real
  numeric_traits<dd_real>::
  to_scalar_t(const qd_real& val)
  { return dd_real(val.x[0],val.x[1]); }

  /// Cast to qd_real
  template< >
  template< >
  inline
  qd_real
  numeric_traits<qd_real>::
  to_scalar_t(const long double& val)
  { return qd_real(double(val),double(val - static_cast<long double>(double(val))),0.0,0.0); }

  /// Cast to qd_real
  template< >
  template< >
  inline
  qd_real
  numeric_traits<qd_real>::
  to_scalar_t(const dd_real& val)
  { return qd_real(val.x[0],val.x[1],0.0,0.0); }
#endif

  /// Cast to promote
  template< typename Scalar_T >
  inline
  typename numeric_traits<Scalar_T>::promoted::type
  to_promote(const Scalar_T& val)
  {
    typedef typename numeric_traits<Scalar_T>::promoted::type promoted_scalar_t;
    return numeric_traits<promoted_scalar_t>::to_scalar_t(val);
  }

  /// Cast to demote
  template< typename Scalar_T >
  inline
  typename numeric_traits<Scalar_T>::demoted::type
  to_demote(const Scalar_T& val)
  {
    typedef typename numeric_traits<Scalar_T>::demoted::type demoted_scalar_t;
    return numeric_traits<demoted_scalar_t>::to_scalar_t(val);
  }
}

#endif // _GLUCAT_SCALAR_IMP_H
