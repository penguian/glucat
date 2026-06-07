#ifndef _GLUCAT_SCALAR_IMP_H
#define _GLUCAT_SCALAR_IMP_H
/**************************************************************************
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
 ***************************************************************************
 ***************************************************************************/

#include "glucat/scalar.h"
#include "glucat/qd.h"



#include <cmath>
#include <limits>

namespace glucat
{
  /*
   * @brief Extra traits which extend numeric limits
   * @details
   */
  // Reference: [AA], 2.4, p. 30-31

  /*
   * @brief Cast to float
   * @details
   * @param val Value
   * @return Converted value
   */
  template< >
  template< typename Other_Scalar_T >
  inline
  auto
  numeric_traits<float>::
  to_scalar_t(const Other_Scalar_T& val) -> float
  { return static_cast<float>(numeric_traits<Other_Scalar_T>::to_double(val)); }

  /*
   * @brief to_scalar_t
   * @details
   *
   * Usage example:
   * Location: glucat/framed_multi_imp.h:125
   *
   * @code
   *
   * this->insert(term_t(val_term.first, numeric_traits<Scalar_T>::to_scalar_t(val_term.second)));
   * @endcode
   *
   * @tparam Other_Scalar_T
   * @param val Value
   * @return Converted value
   */
  template< >
  template< typename Other_Scalar_T >
  inline
  auto
  numeric_traits<double>::
  to_scalar_t(const Other_Scalar_T& val) -> double
  { return numeric_traits<Other_Scalar_T>::to_double(val); }

#if defined(_GLUCAT_USE_QD)
  /*
   * @brief Cast to long double
   * @details
   * @param val Value
   * @return Converted value
   */
  template< >
  template< >
  inline
  auto
  numeric_traits<long double>::
  to_scalar_t(const dd_real& val) -> long double
  { return static_cast<long double>(val.x[0]) + static_cast<long double>(val.x[1]); }

  /*
   * @brief Cast to long double
   * @details
   * @param val Value
   * @return Converted value
   */
  template< >
  template< >
  inline
  auto
  numeric_traits<long double>::
  to_scalar_t(const qd_real& val) -> long double
  { return static_cast<long double>(val.x[0]) + static_cast<long double>(val.x[1]); }

  /*
   * @brief Cast to dd_real
   * @details
   * @param val Value
   * @return Converted value
   */
  template< >
  template< >
  inline
  auto
  numeric_traits<dd_real>::
  to_scalar_t(const long double& val) -> dd_real
  { return {double(val),double(val - static_cast<long double>(double(val)))}; }

  /*
   * @brief Cast to dd_real
   * @details
   * @param val Value
   * @return Converted value
   */
  template< >
  template< >
  inline
  auto
  numeric_traits<dd_real>::
  to_scalar_t(const qd_real& val) -> dd_real
  { return {val.x[0],val.x[1]}; }

  /*
   * @brief Cast to qd_real
   * @details
   * @param val Value
   * @return Converted value
   */
  template< >
  template< >
  inline
  auto
  numeric_traits<qd_real>::
  to_scalar_t(const long double& val) -> qd_real
  { return {double(val),double(val - static_cast<long double>(double(val))),0.0,0.0}; }

  /*
   * @brief Cast to qd_real
   * @details
   * @param val Value
   * @return Converted value
   */
  template< >
  template< >
  inline
  auto
  numeric_traits<qd_real>::
  to_scalar_t(const dd_real& val) -> qd_real
  { return {val.x[0],val.x[1],0.0,0.0}; }
#endif

#if defined(_GLUCAT_USE_QD) && defined(EIGEN_MAJOR_VERSION)
  // disambiguate for Eigen::Index (typically long or ptrdiff_t)
  // qd_real has constructors for int, double, but not long/long long, causing ambiguity.
  // We can inject a cast or a conversion helper if this file is included before the error site.
  // However, the error is inside Eigen code calling qd_real(Index).
  // The only way to fix that without editing Eigen or QD is if *we* control the trait that Eigen uses
  // OR if we can add a constructor to qd_real (we can't easily, library code).
  // Actually, Eigen uses RealScalar(index) cast.
  // If we specialize Eigen::NumTraits<qd_real>, we might control this.
#endif

  /*
   * @brief Cast to promote
   * @details
   *
   * Usage example:
   * Location: glucat/framed_multi_imp.h:2327
   *
   * @code
   *
   * return to_promote(val);
   * @endcode
   *
   * @tparam Scalar_T
   * @param val Value
   * @return Result
   */
  template< typename Scalar_T >
  inline
  auto
  to_promote(const Scalar_T& val) -> typename numeric_traits<Scalar_T>::promoted::type
  {
    using promoted_scalar_t = typename numeric_traits<Scalar_T>::promoted::type;
    return numeric_traits<promoted_scalar_t>::to_scalar_t(val);
  }

  /*
   * @brief Cast to demote
   * @details
   *
   * Usage example:
   * Location: glucat/framed_multi_imp.h:2310
   *
   * @code
   *
   * return to_demote(val);
   * @endcode
   *
   * @tparam Scalar_T
   * @param val Value
   * @return Result
   */
  template< typename Scalar_T >
  inline
  auto
  to_demote(const Scalar_T& val) -> typename numeric_traits<Scalar_T>::demoted::type
  {
    using demoted_scalar_t = typename numeric_traits<Scalar_T>::demoted::type;
    return numeric_traits<demoted_scalar_t>::to_scalar_t(val);
  }
}

#ifdef GLUCAT_DOCTEST
#include <doctest/doctest.h>
#include <complex>

namespace glucat {

  template <typename T>
  void test_scalar_traits() {
    using traits = numeric_traits<T>;

    T zero(0);
    T one(1);

    SUBCASE("Constants") {
      CHECK(traits::pi() > T(3));
      CHECK(traits::pi() < T(4));
      CHECK(traits::ln_2() > T(0.6));
      CHECK(traits::ln_2() < T(0.7));
    }

    SUBCASE("Classification") {
      CHECK_FALSE(traits::isNaN_or_isInf(zero));
      CHECK_FALSE(traits::isNaN_or_isInf(one));
      CHECK_FALSE(traits::isNaN(zero));
      CHECK_FALSE(traits::isInf(zero));

      T nan_val = traits::NaN();
      CHECK(traits::isNaN(nan_val));
    }

    SUBCASE("Math Functions") {
      T val = T(4);
      CHECK(numeric_traits<T>::to_double(traits::sqrt(val)) == doctest::Approx(2.0));

      T e = T(1);
      CHECK(numeric_traits<T>::to_double(traits::exp(e)) == doctest::Approx(std::exp(1.0)));

      CHECK(numeric_traits<T>::to_double(traits::log(traits::exp(one))) == doctest::Approx(1.0));

      T pi_val = traits::pi();
      CHECK(numeric_traits<T>::to_double(traits::sin(pi_val/T(2))) == doctest::Approx(1.0));
      CHECK(numeric_traits<T>::to_double(traits::cos(pi_val)) == doctest::Approx(-1.0));
    }

    SUBCASE("Promote and Demote") {
      T val(2);
      auto p = to_promote(val);
      static_assert(std::is_same_v<decltype(p), typename traits::promoted::type>);

      auto d = to_demote(val);
      static_assert(std::is_same_v<decltype(d), typename traits::demoted::type>);
    }
  }

} // namespace glucat

TEST_CASE("scalar::traits_and_math") {
  SUBCASE("float")       { glucat::test_scalar_traits<float>(); }
  SUBCASE("double")      { glucat::test_scalar_traits<double>(); }
  SUBCASE("long double") { glucat::test_scalar_traits<long double>(); }
#ifdef _GLUCAT_USE_QD
  SUBCASE("dd_real")     { glucat::test_scalar_traits<dd_real>(); }
  SUBCASE("qd_real")     { glucat::test_scalar_traits<qd_real>(); }
#endif

  SUBCASE("is_complex trait") {
    using namespace glucat;
    CHECK_FALSE(is_complex_v<float>);
    CHECK_FALSE(is_complex_v<double>);
    CHECK_FALSE(is_complex_v<long double>);
#ifdef _GLUCAT_USE_QD
    CHECK_FALSE(is_complex_v<dd_real>);
    CHECK_FALSE(is_complex_v<qd_real>);
#endif
    CHECK(is_complex_v<std::complex<double>>);
  }
}
#endif

#endif // _GLUCAT_SCALAR_IMP_H
