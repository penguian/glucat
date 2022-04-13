#ifndef _GLUCAT_CLIFFORD_ALGEBRA_IMP_H
#define _GLUCAT_CLIFFORD_ALGEBRA_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    clifford_algebra_imp.h : Implement common Clifford algebra functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2021 by Paul C. Leopardi
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

// References for algorithms:
// [AS]:
// Milton Abramowicz and Irene A. Stegun, "Handbook of mathematical functions",
// Dover 1972, first published 1965.
// [CHKL]:
// Sheung Hun Cheng, Nicholas J. Higham, Charles S. Kenney and Alan J. Laub,
// "Approximating the Logarithm of a Matrix to Specified Accuracy", 1999.
// ftp://ftp.ma.man.ac.uk/pub/narep/narep353.ps.gz
// [GL]:
// Gene H. Golub and Charles F. van Loan,
// "Matrix Computations", 3rd ed., Johns Hopkins UP, 1996.
// [GW]:
// C.F. Gerald and P.O. Wheatley, "Applied Numerical Analysis",
// 6th Edition, Addison-Wesley, 1999.
// [H]:
// Nicholas J. Higham
// "The Scaling and Squaring Method for the Matrix Exponential Revisited",
// SIAM Journal on Matrix Analysis and Applications,
// Vol. 26,  Issue 4 (2005), pp. 1179-1193.
// [Z]:
// Doron Zeilberger, "PADE" (Maple code), 2002.
// http://www.math.rutgers.edu/~zeilberg/tokhniot/PADE

#include "glucat/clifford_algebra.h"
#include "glucat/scalar.h"

#include <array>

namespace glucat
{
  template< typename Scalar_T, typename Index_Set_T, typename Multivector_T>
  auto
  clifford_algebra<Scalar_T,Index_Set_T,Multivector_T>::
  classname() -> const std::string
  { return "clifford_algebra"; }

  /// Default for truncation
  template< typename Scalar_T, typename Index_Set_T, typename Multivector_T>
  const
  Scalar_T
  clifford_algebra<Scalar_T,Index_Set_T,Multivector_T>::
  default_truncation = std::numeric_limits<Scalar_T>::epsilon();

  /// Test for inequality of multivectors
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  inline
  auto
  operator!= (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> bool
  { return !(lhs == rhs); }

  /// Test for inequality of multivector and scalar
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  operator!= (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const Scalar_T& scr) -> bool
  { return !(lhs == scr); }

  /// Test for inequality of scalar and multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  operator!= (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI,Tune_P>& rhs) -> bool
  { return !(rhs == scr); }

  /// Quadratic norm tolerance
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  auto
  norm_tol(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> Scalar_T
  {
    using multivector_t = Multivector<Scalar_T,LO,HI,Tune_P>;
    static const auto scalar_eps  = std::numeric_limits<Scalar_T>::epsilon();
    static const auto nbr_different_bits =
      std::numeric_limits<Scalar_T>::digits / Tune_P::denom_different_bits + Tune_P::extra_different_bits;
    static const auto abs_tol = scalar_eps *
      numeric_traits<Scalar_T>::pow(Scalar_T(2), nbr_different_bits);
    using framed_multi_t = typename multivector_t::framed_multi_t;
    const auto nbr_terms = double(framed_multi_t(val).truncated(scalar_eps).nbr_terms());
    return abs_tol * abs_tol * std::max(Scalar_T(nbr_terms), Scalar_T(1));
  }

  /// Quadratic norm of difference of multivectors
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  auto
  norm_of_diff(const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> Scalar_T
  {
    const auto rhs_tol = norm_tol(rhs);
    const auto relative = (norm(rhs) > rhs_tol) && (norm(lhs) > rhs_tol);
    const auto abs_norm_diff = norm(rhs-lhs);
    return (relative)
      ? abs_norm_diff/norm(rhs)
      : abs_norm_diff;
  }

  /// Test for approximate equality of multivectors
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  auto
  approx_equal(const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> bool
  { return norm_of_diff(lhs, rhs) < norm_tol(rhs); }

  /// Geometric sum of multivector and scalar
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  operator+ (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const Scalar_T& scr) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    auto result = lhs;
    return result += scr;
  }

  /// Geometric sum of scalar and multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  operator+ (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    return rhs + scr;
  }

  /// Geometric sum
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  inline
  auto
  operator+ (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    auto result = lhs;
    return result += rhs;
  }

  /// Geometric difference of multivector and scalar
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  operator- (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const Scalar_T& scr) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    auto result = lhs;
    return result -= scr;
  }

  /// Geometric difference of scalar and multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  operator- (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return -rhs + scr; }

  /// Geometric difference
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  inline
  auto
  operator- (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    auto result = lhs;
    return result -= rhs;
  }

  /// Product of multivector and scalar
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  operator* (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const Scalar_T& scr) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    auto result = lhs;
    return result *= scr;
  }

  /// Product of scalar and multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  operator* (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { // Note: this assumes that scalar commutes with multivector.
    // This excludes Clifford algebras over non-commuting rings.
    return rhs * scr;
  }

  /// Geometric product
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  inline
  auto
  operator* (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using multivector_t = Multivector<Scalar_T,LO,HI,Tune_P>;
    return lhs * multivector_t(rhs);
  }

  /// Outer product
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  inline
  auto
  operator^ (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using multivector_t = Multivector<Scalar_T,LO,HI,Tune_P>;
    return lhs ^ multivector_t(rhs);
  }

  /// Inner product
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  inline
  auto
  operator& (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using multivector_t = Multivector<Scalar_T,LO,HI,Tune_P>;
    return lhs & multivector_t(rhs);
  }

  /// Left contraction
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  inline
  auto
  operator% (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using multivector_t = Multivector<Scalar_T,LO,HI,Tune_P>;
    return lhs % multivector_t(rhs);
  }

  /// Hestenes scalar product
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  inline
  auto
  star (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> Scalar_T
  {
    using multivector_t = Multivector<Scalar_T,LO,HI,Tune_P>;
    return star(lhs, multivector_t(rhs));
  }

  /// Quotient of multivector and scalar
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  operator/ (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const Scalar_T& scr) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    auto result = lhs;
    return result /= scr;
  }

  /// Quotient of scalar and multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  operator/ (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    Multivector<Scalar_T,LO,HI,Tune_P> result = scr;
    return result /= rhs;
  }

  /// Geometric quotient
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  inline
  auto
  operator/ (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using multivector_t = Multivector<Scalar_T,LO,HI,Tune_P>;
    return lhs / multivector_t(rhs);
  }

  /// Transformation via twisted adjoint action
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  inline
  auto
  operator| (const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using multivector_t = Multivector<Scalar_T,LO,HI,Tune_P>;
    return lhs | multivector_t(rhs);
  }

  /// Geometric multiplicative inverse
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  inv(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return val.inv(); }

  /// Integer power of multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  pow(const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, int rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using multivector_t = Multivector<Scalar_T,LO,HI,Tune_P>;
    if (lhs == Scalar_T(0))
    {
      using traits_t = numeric_traits<Scalar_T>;
      return
        (rhs < 0)
        ? traits_t::NaN()
        : (rhs == 0)
          ? Scalar_T(1)
          : Scalar_T(0);
    }
    auto result = multivector_t(Scalar_T(1));
    auto power =
      (rhs < 0)
      ? lhs.inv()
      : lhs;
    for (auto
        k = std::abs(rhs);
        k != 0;
        k /= 2)
    {
      if (k % 2)
        result *= power;
      power *= power;
    }
    return result;
  }

  /// Multivector power of multivector
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    template<typename, const index_t, const index_t, typename> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  inline
  auto
  pow(const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, const RHS<Scalar_T,LO,HI,Tune_P>& rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;

    if (lhs == Scalar_T(0))
    {
      const Scalar_T m = rhs.scalar();
      if (rhs == m)
        return
          (m < 0)
            ? traits_t::NaN()
            : (m == 0)
              ? Scalar_T(1)
              : Scalar_T(0);
      else
        return Scalar_T(0);
    }
    return exp(log(lhs) * rhs);
  }

  /// Outer product power of multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  outer_pow(const Multivector<Scalar_T,LO,HI,Tune_P>& lhs, int rhs) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return lhs.outer_pow(rhs); }

  /// Scalar part
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  scalar(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> Scalar_T
  { return val.scalar(); }

  /// Real part: synonym for scalar part
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  real(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> Scalar_T
  { return val.scalar(); }

  /// Imaginary part: deprecated (always 0)
  template
  <
    template<typename, const index_t, const index_t, typename> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P
  >
  inline
  auto
  imag(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> Scalar_T
  { return Scalar_T(0); }

  /// Pure part
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  pure(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return val - val.scalar(); }

  /// Even part
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  even(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return val.even(); }

  /// Odd part
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  odd(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return val.odd(); }

  /// Vector part of multivector, as a vector_t with respect to frame()
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  vector_part(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const std::vector<Scalar_T>
  { return val.vector_part(); }

  /// Main involution, each {i} is replaced by -{i} in each term, eg. {1} ->-{1}
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  involute(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return val.involute(); }

  /// Reversion, eg. {1}*{2} -> {2}*{1}
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  reverse(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return val.reverse(); }

  /// Conjugation, rev o invo == invo o rev
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  conj(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return val.conj(); }

  /// Scalar_T quadratic form == (rev(x)*x)(0)
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  quad(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> Scalar_T
  { return val.quad(); }

  /// Scalar_T norm == sum of norm of coordinates
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  norm(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> Scalar_T
  { return val.norm(); }

  /// Absolute value == sqrt(norm)
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  abs(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> Scalar_T
  { return numeric_traits<Scalar_T>::sqrt(val.norm()); }

  /// Maximum of absolute values of components of multivector: multivector infinity norm
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  max_abs(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> Scalar_T
  { return val.max_abs(); }

  /// Square root of -1 which commutes with all members of the frame of the given multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  complexifier(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using multivector_t = Multivector<Scalar_T,LO,HI,Tune_P>;
    using traits_t = numeric_traits<Scalar_T>;

    auto frm = val.frame();
    using array_t = std::array<index_t, 4>;
    auto incp = array_t{0, 2, 1, 0};
    auto incq = array_t{1, 0, 0, 0};
    auto bott = pos_mod((frm.count_pos() - frm.count_neg()), 4);
    for (auto
        k = index_t(0);
        k != incp[bott];
        k++)
      for (auto
          idx = index_t(1);
          idx != HI+1;
          ++idx)
        if (!frm[idx])
        {
          frm.set(idx);
          break;
        }
    for (auto
        k = index_t(0);
        k != incq[bott];
        k++)
      for (auto
          idx = index_t(-1);
          idx != LO-1;
          --idx)
        if (!frm[idx])
        {
          frm.set(idx);
          break;
        }
    auto new_bott = pos_mod(frm.count_pos() - frm.count_neg(), 4);

    if ((incp[new_bott] == 0) && (incq[new_bott] == 0))
      return multivector_t(frm, Scalar_T(1));
    else
      // Return IEEE NaN or -Inf
      return traits_t::NaN();
  }

  /// Square root of -1 which commutes with all members of the frame of the given multivector
  /// The name "elliptic" is now deprecated: use "complexifier" instead.
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  elliptic(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return complexifier(val); }

  /// Check that i is a valid complexifier for val
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  static
  void
  check_complex(const Multivector<Scalar_T,LO,HI,Tune_P>& val,
                const Multivector<Scalar_T,LO,HI,Tune_P>& i, const bool prechecked = false)
  {
    if (!prechecked)
    {
      using multivector_t = Multivector<Scalar_T,LO,HI,Tune_P>;
      using index_set_t = typename multivector_t::index_set_t;
      using error_t = typename multivector_t::error_t;

      const auto i_frame = i.frame();
      // We need i to be a complexifier whose frame is large enough to represent val
      if (complexifier(i) != i ||
         (val.frame() | i_frame) != i_frame ||
          complexifier(val).frame().count() > i_frame.count())
        throw error_t("check_complex(val, i): i is not a valid complexifier for val");
    }
  }

  /// Square root of multivector with specified complexifier
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  sqrt(const Multivector<Scalar_T,LO,HI,Tune_P>& val, const Multivector<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return sqrt(val, i, prechecked); }

  /// Square root of multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  sqrt(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return sqrt(val, complexifier(val), true); }

  /// Exponential of multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  clifford_exp(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    // Scaling and squaring Pade' approximation of matrix exponential
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [H]

    using traits_t = numeric_traits<Scalar_T>;

    const auto scalar_val = val.scalar();
    const auto scalar_exp = traits_t::exp(scalar_val);
    if (traits_t::isNaN_or_isInf(scalar_exp))
      return traits_t::NaN();
    if (val == scalar_val)
      return scalar_exp;

    using multivector_t = Multivector<Scalar_T,LO,HI,Tune_P>;
    auto A = val - scalar_val;
    const auto pure_scale2 = A.norm();

    if (traits_t::isNaN_or_isInf(pure_scale2))
      return traits_t::NaN();
    if (pure_scale2 == Scalar_T(0))
      return scalar_exp;

    const auto ilog2_scale =
      std::max(0, traits_t::to_int(ceil((log2(pure_scale2) + Scalar_T(A.frame().count()))/Scalar_T(2))) - 3);
    const auto i_scale = traits_t::pow(Scalar_T(2), ilog2_scale);
    if (traits_t::isNaN_or_isInf(i_scale))
      return traits_t::NaN();

    A /= i_scale;
    multivector_t pure_exp;
    {
      using limits_t = std::numeric_limits<Scalar_T>;
      const auto nbr_even_powers = 2*(limits_t::digits / 32) + 4;
      using nbr_t = decltype(nbr_even_powers);

      // Create an array of coefficients
      const auto max_power = 2*nbr_even_powers + 1;
      static std::array<Scalar_T, max_power+1> c;
      if (c[0] != Scalar_T(1))
      {
        c[0] = Scalar_T(1);
        for (auto
            k = decltype(max_power)(0);
            k != max_power;
            ++k)
          c[k+1] = c[k]*(max_power-k) / ((2*max_power-k)*(k+1));
      }

      // Create an array of even powers
      std::array<multivector_t, nbr_even_powers> AA;
      AA[0] = A * A;
      AA[1] = AA[0] * AA[0];
      for (auto
        k = nbr_t(2);
        k != nbr_even_powers;
        ++k)
        AA[k] = AA[k-2] * AA[1];

      // Use compensated summation to calculate U and AV
      auto residual = multivector_t();
      auto U = multivector_t(c[0]);
      for (auto
          k = nbr_t(0);
          k != nbr_even_powers;
          ++k)
      {
        const auto& term = AA[k]*c[2*k + 2] - residual;
        const auto& sum = U + term;
        residual = (sum - U) - term;
        U = sum;
      }
      residual = multivector_t();
      auto AV = multivector_t(c[1]);
      for (auto
          k = nbr_t(0);
          k != nbr_even_powers;
          ++k)
      {
        const auto& term = AA[k]*c[2*k + 3] - residual;
        const auto& sum = AV + term;
        residual = (sum - AV) - term;
        AV = sum;
      }
      AV *= A;
      pure_exp = (U+AV) / (U-AV);
    }
    for (auto
        k = decltype(ilog2_scale)(0);
        k != ilog2_scale;
        ++k)
      pure_exp *= pure_exp;
    return pure_exp * scalar_exp;
  }

  /// Natural logarithm of multivector with specified complexifier
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  log(const Multivector<Scalar_T,LO,HI,Tune_P>& val, const Multivector<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return log(val, i, prechecked); }

  /// Natural logarithm of multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  log(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return log(val, complexifier(val), true); }

  /// Hyperbolic cosine of multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  cosh(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    const auto& s = val.scalar();
    if (val == s)
      return traits_t::cosh(s);
    return (exp(val)+exp(-val)) / Scalar_T(2);
  }

  /// Inverse hyperbolic cosine of multivector with specified complexifier
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  acosh(const Multivector<Scalar_T,LO,HI,Tune_P>& val, const Multivector<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    check_complex(val, i, prechecked);
    if (val.isnan())
      return traits_t::NaN();

    const auto radical = sqrt(val*val - Scalar_T(1), i, true);
    return (norm(val + radical) >= norm(val))
           ?  log(val + radical, i, true)
           : -log(val - radical, i, true);
  }

  /// Inverse hyperbolic cosine of multivector
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  acosh(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return acosh(val, complexifier(val), true); }

  /// Cosine of multivector with specified complexifier
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  cos(const Multivector<Scalar_T,LO,HI,Tune_P>& val, const Multivector<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    const auto& s = val.scalar();
    if (val == s)
      return traits_t::cos(s);

    check_complex(val, i, prechecked);

    static const auto& twopi = Scalar_T(2) * traits_t::pi();
    const auto& z = i *
      (val - s + traits_t::fmod(s, twopi));
    return (exp(z)+exp(-z)) / Scalar_T(2);
  }

  /// Cosine of multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  cos(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return cos(val, complexifier(val), true); }

  /// Inverse cosine of multivector with specified complexifier
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  acos(const Multivector<Scalar_T,LO,HI,Tune_P>& val, const Multivector<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    const auto& s = val.scalar();
    if (val == s && traits_t::abs(s) <= Scalar_T(1))
      return traits_t::acos(s);

    check_complex(val, i, prechecked);
    return i * acosh(val, i, true);
  }

  /// Inverse cosine of multivector
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  acos(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return acos(val, complexifier(val), true); }

  /// Hyperbolic sine of multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  sinh(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    const auto& s = val.scalar();
    if (val == s)
      return traits_t::sinh(s);

    return (exp(val)-exp(-val)) / Scalar_T(2);
  }

  /// Inverse hyperbolic sine of multivector with specified complexifier
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  asinh(const Multivector<Scalar_T,LO,HI,Tune_P>& val, const Multivector<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    check_complex(val, i, prechecked);
    if (val.isnan())
      return traits_t::NaN();

    const auto radical = sqrt(val*val + Scalar_T(1), i, true);
    return (norm(val + radical) >= norm(val))
           ?  log( val + radical, i, true)
           : -log(-val + radical, i, true);
  }

  /// Inverse hyperbolic sine of multivector
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  asinh(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return asinh(val, complexifier(val), true); }

  /// Sine of multivector with specified complexifier
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  sin(const Multivector<Scalar_T,LO,HI,Tune_P>& val, const Multivector<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    const auto& s = val.scalar();
    if (val == s)
      return traits_t::sin(s);

    check_complex(val, i, prechecked);

    static const auto& twopi = Scalar_T(2) * traits_t::pi();
    const auto& z = i *
      (val - s + traits_t::fmod(s, twopi));
    return i * (exp(-z)-exp(z)) / Scalar_T(2);
  }

  /// Sine of multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  sin(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return sin(val, complexifier(val), true); }

  /// Inverse sine of multivector with specified complexifier
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  asin(const Multivector<Scalar_T,LO,HI,Tune_P>& val, const Multivector<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    const auto& s = val.scalar();
    if (val == s && traits_t::abs(s) <= Scalar_T(1))
      return traits_t::asin(s);

    check_complex(val, i, prechecked);
    return -i * asinh(i * val, i, true);
  }

  /// Inverse sine of multivector
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  asin(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return asin(val, complexifier(val), true); }

  /// Hyperbolic tangent of multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  tanh(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    const auto& s = val.scalar();
    if (val == s)
      return traits_t::tanh(s);

    return sinh(val) / cosh(val);
  }

  /// Inverse hyperbolic tangent of multivector with specified complexifier
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  atanh(const Multivector<Scalar_T,LO,HI,Tune_P>& val, const Multivector<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    check_complex(val, i, prechecked);
    return val.isnan()
        ? traits_t::NaN()
        : (norm(val + Scalar_T(1)) > norm(val - Scalar_T(1)))
          ? (log(val + Scalar_T(1), i, true) - log(-val + Scalar_T(1), i, true)) / Scalar_T(2)
          : log((val + Scalar_T(1)) / (-val + Scalar_T(1)), i, true) / Scalar_T(2);
  }

  /// Inverse hyperbolic tangent of multivector
  // Reference: [AS], Section 4.6, p86-89
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  atanh(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return atanh(val, complexifier(val), true); }

  /// Tangent of multivector with specified complexifier
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  tan(const Multivector<Scalar_T,LO,HI,Tune_P>& val, const Multivector<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    const auto& s = val.scalar();
    if (val == s)
      return traits_t::tan(s);

    check_complex(val, i, prechecked);
    return sin(val, i, true) / cos(val, i, true);
  }

  /// Tangent of multivector
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  tan(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return tan(val, complexifier(val), true); }

  /// Inverse tangent of multivector with specified complexifier
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  atan(const Multivector<Scalar_T,LO,HI,Tune_P>& val, const Multivector<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    const auto& s = val.scalar();
    if (val == s)
      return traits_t::atan(s);

    check_complex(val, i, prechecked);
    return -i * atanh(i * val, i, true);
  }

  /// Inverse tangent of multivector
  // Reference: [AS], Section 4.4, p79-83
  template< template<typename, const index_t, const index_t, typename> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  atan(const Multivector<Scalar_T,LO,HI,Tune_P>& val) -> const Multivector<Scalar_T,LO,HI,Tune_P>
  { return atan(val, complexifier(val), true); }

}
#endif  // _GLUCAT_CLIFFORD_ALGEBRA_IMP_H
