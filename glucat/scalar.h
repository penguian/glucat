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

#include <boost/numeric/ublas/traits.hpp>

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
    bool
    isInf(const Scalar_T& val)
    {
      return isInf(val,
             bool_to_type< std::numeric_limits<Scalar_T>::has_infinity >() );
    }

    /// Smart isnan
    inline
    static
    bool
    isNaN(const Scalar_T& val)
    {
      return isNaN(val,
             bool_to_type< std::numeric_limits<Scalar_T>::has_quiet_NaN >() );
    }

    /// Smart isnan or isinf
    inline
    static
    bool
    isNaN_or_isInf(const Scalar_T& val)
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
    int
    to_int(const Scalar_T& val)
    { return static_cast<int>(val); }

    /// Cast to double
    inline
    static
    double
    to_double(const Scalar_T& val)
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
    Scalar_T
    fmod(const Scalar_T& lhs, const Scalar_T& rhs)
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
    { return boost::numeric::ublas::type_traits<Scalar_T>::UBLAS_ABS(val); }

    /// Pi
    inline
    static
    Scalar_T
    pi()
    { return Scalar_T(3.14159265358979323); }

    /// log(2)
    inline
    static
    Scalar_T
    ln_2()
    { return Scalar_T(0.693147180559945309); }

    /// Integer power
    inline
    static
    Scalar_T
    pow(const Scalar_T& val, int n)
    { return std::pow(val, n); }

    /// Square root of scalar
    inline
    static
    auto
    sqrt(const Scalar_T& val) -> Scalar_T
    { return boost::numeric::ublas::type_traits<Scalar_T>::UBLAS_SQRT(val); }

    /// Exponential
    inline
    static
    Scalar_T
    exp(const Scalar_T& val)
    { return std::exp(val); }

    /// Logarithm of scalar
    inline
    static
    Scalar_T
    log(const Scalar_T& val)
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
    Scalar_T
    cos(const Scalar_T& val)
    { return std::cos(val); }

    /// Inverse cosine of scalar
    inline
    static
    Scalar_T
    acos(const Scalar_T& val)
    { return std::acos(val); }

    /// Hyperbolic cosine of scalar
    inline
    static
    Scalar_T
    cosh(const Scalar_T& val)
    { return std::cosh(val); }

    /// Sine of scalar
    inline
    static
    Scalar_T
    sin(const Scalar_T& val)
    { return std::sin(val); }

    /// Inverse sine of scalar
    inline
    static
    Scalar_T
    asin(const Scalar_T& val)
    { return std::asin(val); }

    /// Hyperbolic sine of scalar
    inline
    static
    Scalar_T
    sinh(const Scalar_T& val)
    { return std::sinh(val); }

    /// Tangent of scalar
    inline
    static
    Scalar_T
    tan(const Scalar_T& val)
    { return std::tan(val); }

    /// Inverse tangent of scalar
    inline
    static
    Scalar_T
    atan(const Scalar_T& val)
    { return std::atan(val); }

    /// Hyperbolic tangent of scalar
    inline
    static
    Scalar_T
    tanh(const Scalar_T& val)
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
