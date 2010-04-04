#ifndef _GLUCAT_QD_H
#define _GLUCAT_QD_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    qd.h : Define functions for dd_real and qd_real as scalar_t
                             -------------------
    begin                : 2010-03-23
    copyright            : (C) 2010 by Paul C. Leopardi
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

#ifdef QD_API
namespace glucat
{
  /// Extra traits which extend numeric limits
  // Reference: [AA], 2.4, p. 30-31

  /// Smart isnan for dd_real
  template<>
  inline
  bool
  numeric_traits<dd_real>::
  isNaN(const dd_real& val)
  { return val.isnan(); }

  /// Smart isnan for qd_real
  template<>
  inline
  bool
  numeric_traits<qd_real>::
  isNaN(const qd_real& val)
  { return val.isnan(); }

  /// Smart isinf for dd_real
  template<>
  inline
  bool
  numeric_traits<dd_real>::
  isInf(const dd_real& val)
  { return val.isinf(); }

  /// Smart isinf for qd_real
  template<>
  inline
  bool
  numeric_traits<qd_real>::
  isInf(const qd_real& val)
  { return val.isinf(); }

  /// Smart isnan or isinf for dd_real
  template<>
  inline
  bool
  numeric_traits<dd_real>::
  isNaN_or_isInf(const dd_real& val)
  { return val.isnan() || val.isinf(); }

  /// Smart isnan or isinf for qd_real
  template<>
  inline
  bool
  numeric_traits<qd_real>::
  isNaN_or_isInf(const qd_real& val)
  { return val.isnan() || val.isinf(); }

  /// to_int for dd_real
  template<>
  inline
  int
  numeric_traits<dd_real>::
  to_int(const dd_real& val)
  { return ::to_int(val); }

  /// to_int for qd_real
  template<>
  inline
  int
  numeric_traits<qd_real>::
  to_int(const qd_real& val)
  { return ::to_int(val); }

  /// Modulo function for dd_real
  template<>
  inline
  const dd_real
  numeric_traits<dd_real>::
  fmod(const dd_real& lhs, const dd_real& rhs)
  { return ::fmod(lhs, rhs); }

  /// Modulo function for qd_real
  template<>
  inline
  const qd_real
  numeric_traits<qd_real>::
  fmod(const qd_real& lhs, const qd_real& rhs)
  { return ::fmod(lhs, rhs); }

  /// pow for dd_real
  template<>
  inline
  const dd_real
  numeric_traits<dd_real>::
  pow(const dd_real& val, int n)
  { 
    dd_real result = 1.0;
    dd_real pow2 = val;
    for (int k = n; k != 0; pow2 *= pow2, k /= 2)
      if (k % 2)
        result *= pow2;
    return result;
  }

  /// pow for qd_real
  template<>
  inline
  const qd_real
  numeric_traits<qd_real>::
  pow(const qd_real& val, int n)
  { 
    qd_real result = 1.0;
    qd_real pow2 = val;
    for (int k = n; k != 0; pow2 *= pow2, k /= 2)
      if (k % 2)
        result *= pow2;
    return result;
  }

  /// Pi for dd_real
  template<> 
  inline
  const dd_real 
  numeric_traits<dd_real>::
  pi()
  { return dd_real::_pi; }

  /// Pi for qd_real
  template<> 
  inline
  const qd_real 
  numeric_traits<qd_real>::
  pi()
  { return qd_real::_pi; }

  /// log(2) for dd_real
  template<> 
  inline
  const dd_real 
  numeric_traits<dd_real>::
  ln_2()
  { return dd_real::_log2; }

  /// log(2) for qd_real
  template<> 
  inline
  const qd_real 
  numeric_traits<qd_real>::
  ln_2()
  { return qd_real::_log2; }

  /// Log for dd_real
  template<>
  inline
  const dd_real
  numeric_traits<dd_real>::
  log(const dd_real& x)
  { return ::log(x); }

  /// Log for qd_real
  template<>
  inline
  const qd_real
  numeric_traits<qd_real>::
  log(const qd_real& x)
  { return ::log(x); }

  /// exp for dd_real
  template<>
  inline
  const dd_real
  numeric_traits<dd_real>::
  exp(const dd_real& val)
  { return ::exp(val); }

  /// exp for qd_real
  template<>
  inline
  const qd_real
  numeric_traits<qd_real>::
  exp(const qd_real& val)
  { return ::exp(val); }

}
#endif

#endif // _GLUCAT_QD_H
