#ifndef _GLUCAT_SCALAR_H
#define _GLUCAT_SCALAR_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    scalar.h : Define functions for scalar_t
                             -------------------
    begin                : 2001-12-20
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

#include "glucat/portability.h"
#include "glucat/global.h"



#include <cmath>
#include <limits>

namespace glucat
{
  /// Extra traits which extend numeric limits
  // Reference: [AA], 2.4, p. 30-31
  template< typename Scalar_T >
  class numeric_traits
  {
  private:
    /// Smart isinf specialised for Scalar_T without infinity
    inline
    static
    auto
    isInf(const Scalar_T& val, bool_to_type<false>) -> bool
    { return false; }

    /// Smart isinf specialised for Scalar_T with infinity
    inline
    static
    auto
    isInf(const Scalar_T& val, bool_to_type<true>) -> bool
    { return _GLUCAT_ISINF(val); }

    /// Smart isnan specialised for Scalar_T without quiet NaN
    inline
    static
    auto
    isNaN(const Scalar_T& val, bool_to_type<false>) -> bool
    { return false; }

    /// Smart isnan specialised for Scalar_T with quiet NaN
    inline
    static
    auto
    isNaN(const Scalar_T& val, bool_to_type<true>) -> bool
    { return _GLUCAT_ISNAN(val); }

  public:
    /// Smart isinf
    inline
    static
    auto
    isInf(const Scalar_T& val) -> bool
    {
      return isInf(val,
             bool_to_type< std::numeric_limits<Scalar_T>::has_infinity >() );
    }

    /// Smart isnan
    inline
    static
    auto
    isNaN(const Scalar_T& val) -> bool
    {
      return isNaN(val,
             bool_to_type< std::numeric_limits<Scalar_T>::has_quiet_NaN >() );
    }

    /// Smart isnan or isinf
    inline
    static
    auto
    isNaN_or_isInf(const Scalar_T& val) -> bool
    {
      return isNaN(val,
             bool_to_type< std::numeric_limits<Scalar_T>::has_quiet_NaN >() )
          || isInf(val,
             bool_to_type< std::numeric_limits<Scalar_T>::has_infinity >() );
    }

    /// Smart NaN
    inline
    static
    auto
    NaN() -> Scalar_T
    {
      return std::numeric_limits<Scalar_T>::has_quiet_NaN
           ? std::numeric_limits<Scalar_T>::quiet_NaN()
           : Scalar_T(std::log(0.0));
    }

    /// Cast to int
    inline
    static
    auto
    to_int(const Scalar_T& val) -> int
    { return static_cast<int>(val); }

    /// Cast to double
    inline
    static
    auto
    to_double(const Scalar_T& val) -> double
    { return static_cast<double>(val); }

    /// Cast to Scalar_T
    template <typename Other_Scalar_T >
    inline
    static
    auto
    to_scalar_t(const Other_Scalar_T& val) -> Scalar_T
    { return static_cast<Scalar_T>(val); }

    /// Promoted type
    struct promoted {using type = double;};

    /// Demoted type
    struct demoted {using type = float;};

    /// Modulo function for scalar
    inline
    static
    auto
    fmod(const Scalar_T& lhs, const Scalar_T& rhs) -> Scalar_T
    { return std::fmod(lhs, rhs); }

    /// Complex conjugate of scalar
    inline
    static
    auto
    conj(const Scalar_T& val) -> Scalar_T
    { return val; }

    /// Real part of scalar
    inline
    static
    auto
    real(const Scalar_T& val) -> Scalar_T
    { return val; }

    /// Imaginary part of scalar
    inline
    static
    auto
    imag(const Scalar_T& val) -> Scalar_T
    { return Scalar_T(0); }

    /// Absolute value of scalar
    inline
    static
    auto
    abs(const Scalar_T& val) -> Scalar_T
    { return std::abs(val); }

    /// Pi
    inline
    static
    auto
    pi() -> Scalar_T
    { return Scalar_T(3.14159265358979323); }

    /// log(2)
    inline
    static
    auto
    ln_2() -> Scalar_T
    { return Scalar_T(0.693147180559945309); }

    /// Integer power
    inline
    static
    auto
    pow(const Scalar_T& val, int n) -> Scalar_T
    { return std::pow(val, n); }

    /// Square root of scalar
    inline
    static
    auto
    sqrt(const Scalar_T& val) -> Scalar_T
    { return std::sqrt(val); }

    /// Exponential
    inline
    static
    auto
    exp(const Scalar_T& val) -> Scalar_T
    { return std::exp(val); }

    /// Logarithm of scalar
    inline
    static
    auto
    log(const Scalar_T& val) -> Scalar_T
    { return std::log(val); }

    /// Log base 2
    inline
    static
    auto
    log2(const Scalar_T& val) -> Scalar_T
    { return log(val)/ln_2(); }

    /// Cosine of scalar
    inline
    static
    auto
    cos(const Scalar_T& val) -> Scalar_T
    { return std::cos(val); }

    /// Inverse cosine of scalar
    inline
    static
    auto
    acos(const Scalar_T& val) -> Scalar_T
    { return std::acos(val); }

    /// Hyperbolic cosine of scalar
    inline
    static
    auto
    cosh(const Scalar_T& val) -> Scalar_T
    { return std::cosh(val); }

    /// Sine of scalar
    inline
    static
    auto
    sin(const Scalar_T& val) -> Scalar_T
    { return std::sin(val); }

    /// Inverse sine of scalar
    inline
    static
    auto
    asin(const Scalar_T& val) -> Scalar_T
    { return std::asin(val); }

    /// Hyperbolic sine of scalar
    inline
    static
    auto
    sinh(const Scalar_T& val) -> Scalar_T
    { return std::sinh(val); }

    /// Tangent of scalar
    inline
    static
    auto
    tan(const Scalar_T& val) -> Scalar_T
    { return std::tan(val); }

    /// Inverse tangent of scalar
    inline
    static
    auto
    atan(const Scalar_T& val) -> Scalar_T
    { return std::atan(val); }

    /// Hyperbolic tangent of scalar
    inline
    static
    auto
    tanh(const Scalar_T& val) -> Scalar_T
    { return std::tanh(val); }

  };

  /// Log base 2 of scalar
  template< typename Scalar_T >
  inline
  auto
  log2(const Scalar_T& x) -> Scalar_T
  { return numeric_traits<Scalar_T>::log2(x); }
}

#endif // _GLUCAT_SCALAR_H
