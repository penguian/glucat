#ifndef _GLUCAT_CLIFFORD_ALGEBRA_H
#define _GLUCAT_CLIFFORD_ALGEBRA_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    clifford_algebra.h : Declare the operations of a Clifford algebra
                             -------------------
    begin                : Sun 2001-12-09
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

#include <string>
#include <utility>
#include <vector>

namespace glucat
{
  /// clifford_algebra<> declares the operations of a Clifford algebra
  template< typename Scalar_T, typename Index_Set_T, typename Multivector_T>
  class clifford_algebra
  {
  public:
    using scalar_t = Scalar_T;
    using index_set_t = Index_Set_T;
    using multivector_t = Multivector_T;
    using pair_t = std::pair<const index_set_t, Scalar_T>;
    using vector_t = std::vector<Scalar_T>;

    static  auto   classname() -> const std::string;

    virtual ~clifford_algebra() = default;

  // clifford_algebra operations
    /// Test for equality of multivectors
    virtual auto operator==  (const multivector_t& val) const -> bool = 0;
    /// Test for equality of multivector and scalar
    virtual auto operator==  (const Scalar_T& scr) const -> bool = 0;
    /// Geometric sum
    virtual auto operator+=  (const multivector_t& rhs) -> multivector_t& = 0;
    /// Geometric sum of multivector and scalar
    virtual auto operator+=  (const Scalar_T& scr) -> multivector_t& = 0;
    /// Geometric difference
    virtual auto operator-=  (const multivector_t& rhs) -> multivector_t& = 0;
    /// Geometric difference of multivector and scalar
    virtual auto operator-=  (const Scalar_T& scr) -> multivector_t& = 0;
    /// Unary -
    virtual auto operator-   () const -> const multivector_t = 0;
    /// Product of multivector and scalar
    virtual auto operator*=  (const Scalar_T& scr) -> multivector_t& = 0;
    /// Geometric product
    virtual auto operator*=  (const multivector_t& rhs) -> multivector_t& = 0;
    /// Contraction
    virtual auto operator%=  (const multivector_t& rhs) -> multivector_t& = 0;
    /// Inner product
    virtual auto operator&=  (const multivector_t& rhs) -> multivector_t& = 0;
    /// Outer product
    virtual auto operator^=  (const multivector_t& rhs) -> multivector_t& = 0;
    /// Quotient of multivector and scalar
    virtual auto operator/=  (const Scalar_T& scr) -> multivector_t& = 0;
    /// Geometric quotient
    virtual auto operator/=  (const multivector_t& rhs) -> multivector_t& = 0;
    /// Transformation via twisted adjoint action
    virtual auto operator|=  (const multivector_t& rhs) -> multivector_t& = 0;
    /// Geometric multiplicative inverse
    virtual auto inv         () const -> const multivector_t = 0;
    /// *this to the m
    virtual auto pow         (int m) const -> const multivector_t = 0;
    /// Outer product power
    virtual auto outer_pow   (int m) const -> const multivector_t = 0;
    /// Subalgebra generated by all generators of terms of given multivector
    virtual auto frame       () const -> const index_set_t = 0;
    /// Maximum of the grades of each term
    virtual auto grade       () const -> index_t = 0;
    /// Subscripting: map from index set to scalar coordinate
    virtual auto operator[]  (const index_set_t ist) const -> Scalar_T = 0;
    /// Pure grade-vector part
    virtual auto operator()  (index_t grade) const -> const multivector_t = 0;
    /// Scalar part
    virtual auto scalar      () const -> Scalar_T = 0;
    /// Pure part
    virtual auto pure        () const -> const multivector_t = 0;
    /// Even part of multivector, sum of even grade terms
    virtual auto even        () const -> const multivector_t = 0;
    /// Odd part of multivector, sum of odd grade terms
    virtual auto odd         () const -> const multivector_t = 0;
    /// Vector part of multivector, as a vector_t with respect to frame()
    virtual auto vector_part () const -> const vector_t = 0;
    /// Vector part of multivector, as a vector_t with respect to frm
    virtual auto vector_part (const index_set_t frm, const bool prechecked) const -> const vector_t = 0;
    /// Main involution, each {i} is replaced by -{i} in each term, eg. {1} -> -{1}
    virtual auto involute    () const -> const multivector_t = 0;
    /// Reversion, eg. {1}*{2} -> {2}*{1}
    virtual auto reverse     () const -> const multivector_t = 0;
    /// Conjugation, reverse o involute == involute o reverse
    virtual auto conj        () const -> const multivector_t = 0;
    /// Scalar_T quadratic form == (rev(x)*x)(0)
    virtual auto quad        () const -> Scalar_T = 0;
    /// Scalar_T norm == sum of norm of coordinates
    virtual auto norm        () const -> Scalar_T = 0;
    /// Maximum of absolute values of components of multivector: multivector infinity norm
    virtual auto max_abs     () const -> Scalar_T = 0;
    /// Remove all terms with relative size smaller than limit
    virtual auto truncated   (const Scalar_T& limit = Scalar_T(DEFAULT_TRUNCATION)) const -> const multivector_t = 0;
    /// Check if a multivector contains any IEEE NaN values
    virtual auto isnan       () const -> bool = 0;
    /// Write formatted multivector to output
    virtual void write       (const std::string& msg="") const = 0;
    /// Write formatted multivector to file
    virtual void write       (std::ofstream& ofile, const std::string& msg="") const = 0;
  };

#ifndef _GLUCAT_CLIFFORD_ALGEBRA_OPERATIONS
#define _GLUCAT_CLIFFORD_ALGEBRA_OPERATIONS \
    auto operator==  (const multivector_t& val) const  -> bool                override;\
    auto operator==  (const Scalar_T& scr) const       -> bool                override;\
    auto operator+=  (const multivector_t& rhs)        -> multivector_t&      override;\
    auto operator+=  (const Scalar_T& scr)             -> multivector_t&      override;\
    auto operator-=  (const multivector_t& rhs)        -> multivector_t&      override;\
    auto operator-=  (const Scalar_T& scr)             -> multivector_t&      override;\
    auto operator-   () const                          -> const multivector_t override;\
    auto operator*=  (const Scalar_T& scr)             -> multivector_t&      override;\
    auto operator*=  (const multivector_t& rhs)        -> multivector_t&      override;\
    auto operator%=  (const multivector_t& rhs)        -> multivector_t&      override;\
    auto operator&=  (const multivector_t& rhs)        -> multivector_t&      override;\
    auto operator^=  (const multivector_t& rhs)        -> multivector_t&      override;\
    auto operator/=  (const Scalar_T& scr)             -> multivector_t&      override;\
    auto operator/=  (const multivector_t& rhs)        -> multivector_t&      override;\
    auto operator|=  (const multivector_t& rhs)        -> multivector_t&      override;\
    auto inv         () const                          -> const multivector_t override;\
    auto pow         (int m) const                     -> const multivector_t override;\
    auto outer_pow   (int m) const                     -> const multivector_t override;\
    auto frame       () const                          -> const index_set_t   override;\
    auto grade       () const                          -> index_t             override;\
    auto operator[]  (const index_set_t ist) const     -> Scalar_T            override;\
    auto operator()  (index_t grade) const             -> const multivector_t override;\
    auto scalar      () const                          -> Scalar_T            override;\
    auto pure        () const                          -> const multivector_t override;\
    auto even        () const                          -> const multivector_t override;\
    auto odd         () const                          -> const multivector_t override;\
    auto vector_part () const                          -> const vector_t      override;\
    auto vector_part (const index_set_t frm, const bool prechecked = false) const      \
                                                       -> const vector_t      override;\
    auto involute    () const                          -> const multivector_t override;\
    auto reverse     () const                          -> const multivector_t override;\
    auto conj        () const                          -> const multivector_t override;\
    auto quad        () const                          -> Scalar_T            override;\
    auto norm        () const                          -> Scalar_T            override;\
    auto max_abs     () const                          -> Scalar_T            override;\
    auto truncated   (const Scalar_T& limit = Scalar_T(DEFAULT_TRUNCATION)) const      \
                                                       -> const multivector_t override;\
    auto isnan       () const                          -> bool                override;\
    void write       (const std::string& msg="") const                        override;\
    void write       (std::ofstream& ofile, const std::string& msg="") const  override;
#endif // _GLUCAT_CLIFFORD_ALGEBRA_OPERATIONS

  /// Test for inequality of multivectors
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator!= (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs) -> bool;

  /// Test for inequality of multivector and scalar
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator!= (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr) -> bool;

  /// Test for inequality of scalar and multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator!= (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs) -> bool;

  /// Geometric sum of multivector and scalar
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator+ (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr) -> const Multivector<Scalar_T,LO,HI>;

  /// Geometric sum of scalar and multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator+ (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Geometric sum
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator+ (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Geometric difference of multivector and scalar
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator- (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr) -> const Multivector<Scalar_T,LO,HI>;

  /// Geometric difference of scalar and multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator- (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Geometric difference
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator- (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Product of multivector and scalar
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator* (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr) -> const Multivector<Scalar_T,LO,HI>;

  /// Product of scalar and multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator* (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Geometric product
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator* (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Outer product
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator^ (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Inner product
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator& (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Left contraction
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator% (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Hestenes scalar product
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  star (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs) -> Scalar_T;

  /// Quotient of multivector and scalar
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator/ (const Multivector<Scalar_T,LO,HI>& lhs, const Scalar_T& scr) -> const Multivector<Scalar_T,LO,HI>;

  /// Quotient of scalar and multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator/ (const Scalar_T& scr, const Multivector<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Geometric quotient
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator/ (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Transformation via twisted adjoint action
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  operator| (const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Geometric multiplicative inverse
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  inv(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Integer power of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  pow(const Multivector<Scalar_T,LO,HI>& lhs, int rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Multivector power of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    template<typename, const index_t, const index_t> class RHS,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  pow(const Multivector<Scalar_T,LO,HI>& lhs, const RHS<Scalar_T,LO,HI>& rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Outer product power of multivector
  template< template<typename, const index_t, const index_t> class Multivector,
            typename Scalar_T, const index_t LO, const index_t HI >
  auto
  outer_pow(const Multivector<Scalar_T,LO,HI>& lhs, int rhs) -> const Multivector<Scalar_T,LO,HI>;

  /// Scalar part
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  scalar(const Multivector<Scalar_T,LO,HI>& val) -> Scalar_T;

  /// Real part: synonym for scalar part
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  real(const Multivector<Scalar_T,LO,HI>& val) -> Scalar_T;

  /// Imaginary part: deprecated (always 0)
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  imag(const Multivector<Scalar_T,LO,HI>& val) -> Scalar_T;

  /// Pure part
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  pure(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Even part
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  even(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Odd part
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  odd(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Vector part of multivector, as a vector_t with respect to frame()
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  vector_part(const Multivector<Scalar_T,LO,HI>& val) -> const std::vector<Scalar_T>;

  /// Main involution, each {i} is replaced by -{i} in each term, eg. {1}*{2} -> (-{2})*(-{1})
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  involute(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Reversion, eg. {1}*{2} -> {2}*{1}
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  reverse(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Conjugation, rev o invo == invo o rev
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  conj(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Scalar_T quadratic form == (rev(x)*x)(0)
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  quad(const Multivector<Scalar_T,LO,HI>& val) -> Scalar_T;

  /// Scalar_T norm == sum of norm of coordinates
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  norm(const Multivector<Scalar_T,LO,HI>& val) -> Scalar_T;

  /// Absolute value == sqrt(norm)
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  abs(const Multivector<Scalar_T,LO,HI>& val) -> Scalar_T;

  /// Maximum of absolute values of components of multivector: multivector infinity norm
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  max_abs(const Multivector<Scalar_T,LO,HI>& val) -> Scalar_T;

  /// Square root of -1 which commutes with all members of the frame of the given multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  complexifier(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Square root of -1 which commutes with all members of the frame of the given multivector
  /// The name "elliptic" is now deprecated: use "complexifier" instead.
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  elliptic(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Square root of multivector with specified complexifier
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  sqrt(const Multivector<Scalar_T,LO,HI>& val,
       const Multivector<Scalar_T,LO,HI>& i,
       const bool prechecked = false) -> const Multivector<Scalar_T,LO,HI>;

  /// Square root of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  sqrt(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  // Transcendental functions

  /// Exponential of multivector
  template
  < template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  clifford_exp(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Natural logarithm of multivector with specified complexifier
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  log(const Multivector<Scalar_T,LO,HI>& val,
      const Multivector<Scalar_T,LO,HI>& i,
      const bool prechecked = false) -> const Multivector<Scalar_T,LO,HI>;

  /// Natural logarithm of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  log(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Cosine of multivector with specified complexifier
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  cos(const Multivector<Scalar_T,LO,HI>& val,
      const Multivector<Scalar_T,LO,HI>& i,
      const bool prechecked = false) -> const Multivector<Scalar_T,LO,HI>;

  /// Cosine of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  cos(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Inverse cosine of multivector with specified complexifier
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  acos(const Multivector<Scalar_T,LO,HI>& val,
       const Multivector<Scalar_T,LO,HI>& i,
       const bool prechecked = false) -> const Multivector<Scalar_T,LO,HI>;

  /// Inverse cosine of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  acos(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Hyperbolic cosine of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  cosh(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Inverse hyperbolic cosine of multivector with specified complexifier
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  acosh(const Multivector<Scalar_T,LO,HI>& val,
        const Multivector<Scalar_T,LO,HI>& i,
        const bool prechecked = false) -> const Multivector<Scalar_T,LO,HI>;

  /// Inverse hyperbolic cosine of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  acosh(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Sine of multivector with specified complexifier
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  sin(const Multivector<Scalar_T,LO,HI>& val,
      const Multivector<Scalar_T,LO,HI>& i,
      const bool prechecked = false) -> const Multivector<Scalar_T,LO,HI>;

  /// Sine of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  sin(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Inverse sine of multivector with specified complexifier
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  asin(const Multivector<Scalar_T,LO,HI>& val,
       const Multivector<Scalar_T,LO,HI>& i,
       const bool prechecked = false) -> const Multivector<Scalar_T,LO,HI>;

  /// Inverse sine of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  asin(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Hyperbolic sine of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  sinh(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Inverse hyperbolic sine of multivector with specified complexifier
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  asinh(const Multivector<Scalar_T,LO,HI>& val,
        const Multivector<Scalar_T,LO,HI>& i,
        const bool prechecked = false) -> const Multivector<Scalar_T,LO,HI>;

  /// Inverse hyperbolic sine of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  asinh(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Tangent of multivector with specified complexifier
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  tan(const Multivector<Scalar_T,LO,HI>& val,
      const Multivector<Scalar_T,LO,HI>& i,
      const bool prechecked = false) -> const Multivector<Scalar_T,LO,HI>;

  /// Tangent of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  tan(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Inverse tangent of multivector with specified complexifier
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  atan(const Multivector<Scalar_T,LO,HI>& val,
       const Multivector<Scalar_T,LO,HI>& i,
       const bool prechecked = false) -> const Multivector<Scalar_T,LO,HI>;

  /// Inverse tangent of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  atan(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Hyperbolic tangent of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  tanh(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;

  /// Inverse hyperbolic tangent of multivector with specified complexifier
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  atanh(const Multivector<Scalar_T,LO,HI>& val,
        const Multivector<Scalar_T,LO,HI>& i,
        const bool prechecked = false) -> const Multivector<Scalar_T,LO,HI>;

  /// Inverse hyperbolic tangent of multivector
  template
  <
    template<typename, const index_t, const index_t> class Multivector,
    typename Scalar_T, const index_t LO, const index_t HI
  >
  auto
  atanh(const Multivector<Scalar_T,LO,HI>& val) -> const Multivector<Scalar_T,LO,HI>;
}
#endif  // _GLUCAT_CLIFFORD_ALGEBRA_H
