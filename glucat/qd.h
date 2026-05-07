#ifndef _GLUCAT_QD_H
#define _GLUCAT_QD_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    qd.h : Define functions for dd_real and qd_real as scalar_t
                             -------------------
    begin                : 2010-03-23
    copyright            : (C) 2010-2016 by Paul C. Leopardi
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
#include <type_traits>

#if defined(_GLUCAT_USE_QD)
# include <qd/qd_real.h>
#endif

namespace glucat
{
  // Extra traits which extend numeric limits
  // Reference: [AA], 2.4, p. 30-31

#if defined(_GLUCAT_USE_QD) && defined(QD_API)

  // Macro to apply function _F to type _T
# define _GLUCAT_QD_F(_T, _F) \
  template<>              \
  inline                  \
  auto                    \
  numeric_traits<_T>::    \
  _F(const _T& val) -> _T \
  { return ::_F(val); }

  // Smart isnan for dd_real
  template<>
  inline
  auto
  numeric_traits<dd_real>::
  isNaN(const dd_real& val) -> bool
  { return val.isnan(); }

  // Smart isinf for dd_real
  template<>
  inline
  auto
  numeric_traits<dd_real>::
  isInf(const dd_real& val) -> bool
  { return val.isinf(); }

  // Smart isnan or isinf for dd_real
  template<>
  inline
  auto
  numeric_traits<dd_real>::
  isNaN_or_isInf(const dd_real& val) -> bool
  { return val.isnan() || val.isinf(); }

  // to_int for dd_real
  template<>
  inline
  auto
  numeric_traits<dd_real>::
  to_int(const dd_real& val) -> int
  { return ::to_int(val); }

  // to_double for dd_real
  template<>
  inline
  auto
  numeric_traits<dd_real>::
  to_double(const dd_real& val) -> double
  { return ::to_double(val); }

  // to_scalar_t for dd_real
  template<>
  template<typename Other_Scalar_T>
  inline
  auto
  numeric_traits<dd_real>::
  to_scalar_t(const Other_Scalar_T& val) -> dd_real
  {
    if constexpr (std::is_same_v<Other_Scalar_T, dd_real>) {
      return val;
    } else if constexpr (std::is_same_v<Other_Scalar_T, qd_real>) {
      return ::to_dd_real(val);
    } else if constexpr (std::is_integral_v<Other_Scalar_T>) {
      return dd_real(static_cast<double>(val));
    } else {
      return static_cast<dd_real>(val);
    }
  }

  // Modulo function for dd_real
  template<>
  inline
  auto
  numeric_traits<dd_real>::
  fmod(const dd_real& lhs, const dd_real& rhs) -> dd_real
  { return ::fmod(lhs, rhs); }

  // pow for dd_real
  template<>
  inline
  auto
  numeric_traits<dd_real>::
  pow(const dd_real& val, int n) -> dd_real
  {
    if (val == dd_real(0))
    {
      return
        (n < 0)
        ? NaN()
        : (n == 0)
          ? dd_real(1)
          : dd_real(0);
    }
    auto result = dd_real(1);
    auto power =
      (n < 0)
      ? dd_real(1)/val
      : val;
    for (auto
        k = std::abs(n);
        k != 0;
        k /= 2)
    {
      if (k % 2)
        result *= power;
      power *= power;
    }
    return result;
  }

  // Pi for dd_real
  template<>
  inline
  auto
  numeric_traits<dd_real>::
  pi() -> dd_real
  { return dd_real::_pi; }

  // log(2) for dd_real
  template<>
  inline
  auto
  numeric_traits<dd_real>::
  ln_2() -> dd_real
  { return dd_real::_log2; }

  // Exp of dd_real
  _GLUCAT_QD_F(dd_real, exp)

  // Log of dd_real
  _GLUCAT_QD_F(dd_real, log)

  // Cosine of dd_real
  _GLUCAT_QD_F(dd_real, cos)

  // Inverse cosine of dd_real
  _GLUCAT_QD_F(dd_real, acos)

  // Hyperbolic cosine of dd_real
  _GLUCAT_QD_F(dd_real, cosh)

  // Sine of dd_real
  _GLUCAT_QD_F(dd_real, sin)

  // Inverse sine of dd_real
  _GLUCAT_QD_F(dd_real, asin)

  // Hyperbolic sine of dd_real
  _GLUCAT_QD_F(dd_real, sinh)

  // Tangent of dd_real
  _GLUCAT_QD_F(dd_real, tan)

  // Inverse tangent of dd_real
  _GLUCAT_QD_F(dd_real, atan)

  // Hyperbolic tangent of dd_real
  _GLUCAT_QD_F(dd_real, tanh)

  // Smart isnan for qd_real
  template<>
  inline
  auto
  numeric_traits<qd_real>::
  isNaN(const qd_real& val) -> bool
  { return val.isnan(); }

  // Smart isinf for qd_real
  template<>
  inline
  auto
  numeric_traits<qd_real>::
  isInf(const qd_real& val) -> bool
  { return val.isinf(); }

  // Smart isnan or isinf for qd_real
  template<>
  inline
  auto
  numeric_traits<qd_real>::
  isNaN_or_isInf(const qd_real& val) -> bool
  { return val.isnan() || val.isinf(); }

  // to_int for qd_real
  template<>
  inline
  auto
  numeric_traits<qd_real>::
  to_int(const qd_real& val) -> int
  { return ::to_int(val); }

  // to_double for qd_real
  template<>
  inline
  auto
  numeric_traits<qd_real>::
  to_double(const qd_real& val) -> double
  { return ::to_double(val); }

  // to_scalar_t for qd_real
  template<>
  template<typename Other_Scalar_T>
  inline
  auto
  numeric_traits<qd_real>::
  to_scalar_t(const Other_Scalar_T& val) -> qd_real
  {
    if constexpr (std::is_same_v<Other_Scalar_T, qd_real>) {
      return val;
    } else if constexpr (std::is_integral_v<Other_Scalar_T>) {
      return qd_real(static_cast<double>(val));
    } else {
      return static_cast<qd_real>(val);
    }
  }

  // Modulo function for qd_real
  template<>
  inline
  auto
  numeric_traits<qd_real>::
  fmod(const qd_real& lhs, const qd_real& rhs) -> qd_real
  { return ::fmod(lhs, rhs); }

  // pow for qd_real
  template<>
  inline
  auto
  numeric_traits<qd_real>::
  pow(const qd_real& val, int n) -> qd_real
  {
    if (val == qd_real(0))
    {
      return
        (n < 0)
        ? NaN()
        : (n == 0)
          ? qd_real(1)
          : qd_real(0);
    }
    auto result = qd_real(1);
    auto power =
      (n < 0)
      ? qd_real(1)/val
      : val;
    for (auto
        k = std::abs(n);
        k != 0;
        k /= 2)
    {
      if (k % 2)
        result *= power;
      power *= power;
    }
    return result;
  }

  // Pi for qd_real
  template<>
  inline
  auto
  numeric_traits<qd_real>::
  pi() -> qd_real
  { return qd_real::_pi; }

  // log(2) for qd_real
  template<>
  inline
  auto
  numeric_traits<qd_real>::
  ln_2() -> qd_real
  { return qd_real::_log2; }

  // Exp of qd_real
  _GLUCAT_QD_F(qd_real, exp)

  // Log of qd_real
  _GLUCAT_QD_F(qd_real, log)

  // Cosine of qd_real
  _GLUCAT_QD_F(qd_real, cos)

  // Inverse cosine of qd_real
  _GLUCAT_QD_F(qd_real, acos)

  // Hyperbolic cosine of qd_real
  _GLUCAT_QD_F(qd_real, cosh)

  // Sine of qd_real
  _GLUCAT_QD_F(qd_real, sin)

  // Inverse sine of qd_real
  _GLUCAT_QD_F(qd_real, asin)

  // Hyperbolic sine of qd_real
  _GLUCAT_QD_F(qd_real, sinh)

  // Tangent of qd_real
  _GLUCAT_QD_F(qd_real, tan)

  // Inverse tangent of qd_real
  _GLUCAT_QD_F(qd_real, atan)

  // Hyperbolic tangent of qd_real
  _GLUCAT_QD_F(qd_real, tanh)

#endif // !defined(_GLUCAT_USE_QD) || !defined(QD_API)

} // namespace glucat

#if defined(_GLUCAT_USE_QD) && defined(QD_API)
namespace Eigen {
  template<typename T> struct GenericNumTraits;
  template<typename T> struct NumTraits;

  // Wrapper types for Eigen's internal Real scalar to resolve ambiguity.
  // We use composition to completely hide the base class constructors.
  struct qd_real_eigen {
    ::qd_real val;
    qd_real_eigen() : val(0.0) {}
    qd_real_eigen(const ::qd_real& x) : val(x) {}
    qd_real_eigen(double x) : val(x) {}
    qd_real_eigen(int x) : val(x) {}
    qd_real_eigen(long x) : val(static_cast<double>(x)) {}
    qd_real_eigen(unsigned long x) : val(static_cast<double>(x)) {}
    qd_real_eigen(long long x) : val(static_cast<double>(x)) {}
    qd_real_eigen(unsigned long long x) : val(static_cast<double>(x)) {}

    operator ::qd_real&() { return val; }
    operator const ::qd_real&() const { return val; }

    qd_real_eigen operator-() const { return qd_real_eigen(-val); }
    qd_real_eigen& operator+=(const qd_real_eigen& other) { val += other.val; return *this; }
    qd_real_eigen& operator-=(const qd_real_eigen& other) { val -= other.val; return *this; }
    qd_real_eigen& operator*=(const qd_real_eigen& other) { val *= other.val; return *this; }
    qd_real_eigen& operator/=(const qd_real_eigen& other) { val /= other.val; return *this; }
    friend qd_real_eigen operator*(const qd_real_eigen& a, const qd_real_eigen& b) { return qd_real_eigen(a.val * b.val); }
    friend qd_real_eigen operator+(const qd_real_eigen& a, const qd_real_eigen& b) { return qd_real_eigen(a.val + b.val); }
    friend qd_real_eigen operator-(const qd_real_eigen& a, const qd_real_eigen& b) { return qd_real_eigen(a.val - b.val); }
    friend qd_real_eigen operator/(const qd_real_eigen& a, const qd_real_eigen& b) { return qd_real_eigen(a.val / b.val); }
    friend bool operator<(const qd_real_eigen& a, const qd_real_eigen& b) { return a.val < b.val; }
    friend bool operator==(const qd_real_eigen& a, const qd_real_eigen& b) { return a.val == b.val; }
    friend bool operator!=(const qd_real_eigen& a, const qd_real_eigen& b) { return a.val != b.val; }
    friend bool operator==(const qd_real_eigen& a, int b) { return a.val == b; }
    friend bool operator!=(const qd_real_eigen& a, int b) { return a.val != b; }
  };

  struct dd_real_eigen {
    ::dd_real val;
    dd_real_eigen() : val(0.0) {}
    dd_real_eigen(const ::dd_real& x) : val(x) {}
    dd_real_eigen(double x) : val(x) {}
    dd_real_eigen(int x) : val(x) {}
    dd_real_eigen(long x) : val(static_cast<double>(x)) {}
    dd_real_eigen(unsigned long x) : val(static_cast<double>(x)) {}
    dd_real_eigen(long long x) : val(static_cast<double>(x)) {}
    dd_real_eigen(unsigned long long x) : val(static_cast<double>(x)) {}

    operator ::dd_real&() { return val; }
    operator const ::dd_real&() const { return val; }

    dd_real_eigen operator-() const { return dd_real_eigen(-val); }
    dd_real_eigen& operator+=(const dd_real_eigen& other) { val += other.val; return *this; }
    dd_real_eigen& operator-=(const dd_real_eigen& other) { val -= other.val; return *this; }
    dd_real_eigen& operator*=(const dd_real_eigen& other) { val *= other.val; return *this; }
    dd_real_eigen& operator/=(const dd_real_eigen& other) { val /= other.val; return *this; }
    friend dd_real_eigen operator*(const dd_real_eigen& a, const dd_real_eigen& b) { return dd_real_eigen(a.val * b.val); }
    friend dd_real_eigen operator+(const dd_real_eigen& a, const dd_real_eigen& b) { return dd_real_eigen(a.val + b.val); }
    friend dd_real_eigen operator-(const dd_real_eigen& a, const dd_real_eigen& b) { return dd_real_eigen(a.val - b.val); }
    friend dd_real_eigen operator/(const dd_real_eigen& a, const dd_real_eigen& b) { return dd_real_eigen(a.val / b.val); }
    friend bool operator<(const dd_real_eigen& a, const dd_real_eigen& b) { return a.val < b.val; }
    friend bool operator==(const dd_real_eigen& a, const dd_real_eigen& b) { return a.val == b.val; }
    friend bool operator!=(const dd_real_eigen& a, const dd_real_eigen& b) { return a.val != b.val; }
    friend bool operator==(const dd_real_eigen& a, int b) { return a.val == b; }
    friend bool operator!=(const dd_real_eigen& a, int b) { return a.val != b; }
  };

  template<> struct NumTraits<qd_real> {
    typedef qd_real_eigen Real;
    typedef qd_real NonInteger;
    typedef qd_real Literal;
    enum { IsComplex = 0, IsInteger = 0, IsSigned = 1, RequireInitialization = 1, ReadCost = 4, AddCost = 16, MulCost = 16 };
    static inline Real epsilon() { return qd_real::_eps; }
    static inline Real dummy_precision() { return qd_real::_eps * 1e3; }
    static inline Real highest() { return qd_real::_max; }
    static inline Real lowest() { return -qd_real::_max; }
  };

  template<> struct NumTraits<dd_real> {
    typedef dd_real_eigen Real;
    typedef dd_real NonInteger;
    typedef dd_real Literal;
    enum { IsComplex = 0, IsInteger = 0, IsSigned = 1, RequireInitialization = 1, ReadCost = 2, AddCost = 8, MulCost = 8 };
    static inline Real epsilon() { return dd_real::_eps; }
    static inline Real dummy_precision() { return dd_real::_eps * 1e3; }
    static inline Real highest() { return dd_real::_max; }
    static inline Real lowest() { return -dd_real::_max; }
  };

  // Explicit specializations of Eigen templates to avoid ambiguity
  template<typename T> inline typename NumTraits<T>::Real abs(const T &x);
  template<typename T> inline typename NumTraits<T>::Real real(const T &x);
  template<typename T> inline typename NumTraits<T>::Real imag(const T &x);

  template<> inline qd_real_eigen abs(const ::qd_real &x) { return ::abs(x); }
  template<> inline dd_real_eigen abs(const ::dd_real &x) { return ::abs(x); }
  template<> inline qd_real_eigen real(const ::qd_real &x) { return x; }
  template<> inline dd_real_eigen real(const ::dd_real &x) { return x; }
  template<> inline qd_real_eigen imag(const ::qd_real &) { return 0.0; }
  template<> inline dd_real_eigen imag(const ::dd_real &) { return 0.0; }

  // Math functions for wrappers
  inline qd_real_eigen abs(const qd_real_eigen& x) { return ::abs(x.val); }
  inline dd_real_eigen abs(const dd_real_eigen& x) { return ::abs(x.val); }
  inline qd_real_eigen real(const qd_real_eigen& x) { return x.val; }
  inline dd_real_eigen real(const dd_real_eigen& x) { return x.val; }
  inline qd_real_eigen imag(const qd_real_eigen&) { return 0.0; }
  inline dd_real_eigen imag(const dd_real_eigen&) { return 0.0; }
}
#endif

#if defined(_GLUCAT_USE_QD) && defined(GLUCAT_DOCTEST)
#include <doctest.h>
namespace doctest {
  inline bool operator==(const dd_real& lhs, const Approx& rhs) {
    return operator==(to_double(lhs), rhs);
  }
  inline bool operator==(const qd_real& lhs, const Approx& rhs) {
    return operator==(to_double(lhs), rhs);
  }
  inline bool operator==(const Approx& lhs, const dd_real& rhs) {
    return operator==(lhs, to_double(rhs));
  }
  inline bool operator==(const Approx& lhs, const qd_real& rhs) {
    return operator==(lhs, to_double(rhs));
  }
}
#endif

#endif // _GLUCAT_QD_H
