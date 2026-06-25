#ifndef _GLUCAT_MATRIX_MULTI_H
#define _GLUCAT_MATRIX_MULTI_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_multi.h : Declare a class for the matrix representation of a multivector
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

#include <boost/yap/algorithm.hpp>
#include <boost/yap/expression.hpp>
#include <string>
#include <utility>
#include <vector>

#include "glucat/clifford_algebra.h"
#include "glucat/errors.h"
#include "glucat/framed_multi.h"
#include "glucat/global.h"
#include "glucat/index_set.h"
#include "glucat/matrix.h"
#include "glucat/tuning.h"

namespace glucat
{
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  class matrix_multi;

  template <boost::yap::expr_kind Kind, typename Tuple>
  struct matrix_multi_expr;

  namespace detail
  {
    template <typename U>
    struct identity_wrapper
    {
      using type = U;
    };

    // --- find_scalar_type ---
    template <typename T, typename = void>
    struct find_scalar_type
    {
      using type = T;
    };

    template <typename T>
    struct find_scalar_type<T&> : find_scalar_type<T>
    {
    };

    template <typename T>
    struct find_scalar_type<const T&> : find_scalar_type<T>
    {
    };

    template <typename T>
    struct find_scalar_type<T*> : find_scalar_type<T>
    {
    };

    template <typename T>
    struct find_scalar_type<const T*> : find_scalar_type<T>
    {
    };

    template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
    struct find_scalar_type<matrix_multi<Scalar_T, LO, HI, Tune_P>>
    {
      using type = Scalar_T;
    };

    template <typename T, typename = void>
    struct has_scalar_type : std::false_type
    {
    };

    template <typename T>
    struct has_scalar_type<T, std::void_t<typename find_scalar_type<T>::type>> : std::true_type
    {
    };

    template <typename... Args>
    struct find_scalar_in_list;

    template <>
    struct find_scalar_in_list<>
    {
      using type = double;
    };

    template <typename First, typename... Rest>
    struct find_scalar_in_list<First, Rest...>
    {
      using type =
          typename std::conditional_t<has_scalar_type<First>::value, find_scalar_type<First>, find_scalar_in_list<Rest...>>::type;
    };

    template <typename Tuple>
    struct find_scalar_in_tuple;

    template <typename... Args>
    struct find_scalar_in_tuple<boost::hana::tuple<Args...>> : find_scalar_in_list<Args...>
    {
    };

    template <boost::yap::expr_kind Kind, typename Tuple>
    struct find_scalar_type<matrix_multi_expr<Kind, Tuple>> : find_scalar_in_tuple<std::decay_t<Tuple>>
    {
    };

    // --- find_matrix_multi ---
    template <typename T, typename = void>
    struct find_matrix_multi
    {
    };

    template <typename T>
    struct find_matrix_multi<T&> : find_matrix_multi<T>
    {
    };

    template <typename T>
    struct find_matrix_multi<const T&> : find_matrix_multi<T>
    {
    };

    template <typename T>
    struct find_matrix_multi<T*> : find_matrix_multi<T>
    {
    };

    template <typename T>
    struct find_matrix_multi<const T*> : find_matrix_multi<T>
    {
    };

    template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
    struct find_matrix_multi<matrix_multi<Scalar_T, LO, HI, Tune_P>>
    {
      using type = matrix_multi<Scalar_T, LO, HI, Tune_P>;
    };

    template <typename T, typename = void>
    struct has_matrix_multi : std::false_type
    {
    };

    template <typename T>
    struct has_matrix_multi<T, std::void_t<typename find_matrix_multi<T>::type>> : std::true_type
    {
    };

    template <typename... Args>
    struct find_in_list;

    template <>
    struct find_in_list<>
    {
    };

    template <typename First, typename... Rest>
    struct find_in_list<First, Rest...>
        : std::conditional_t<has_matrix_multi<First>::value, find_matrix_multi<First>, find_in_list<Rest...>>
    {
    };

    template <typename Tuple>
    struct find_in_tuple;

    template <typename... Args>
    struct find_in_tuple<boost::hana::tuple<Args...>> : find_in_list<Args...>
    {
    };

    template <boost::yap::expr_kind Kind, typename Tuple>
    struct find_matrix_multi<matrix_multi_expr<Kind, Tuple>> : find_in_tuple<std::decay_t<Tuple>>
    {
    };

    // --- get_terminal ---
    template <typename T>
    struct get_terminal
    {
      using scalar_type = typename find_scalar_type<T>::type;
      using type = typename std::conditional_t<
          has_matrix_multi<T>::value, find_matrix_multi<T>,
          identity_wrapper<matrix_multi<std::decay_t<scalar_type>, DEFAULT_LO, DEFAULT_HI, tuning<>>>>::type;
    };
  }  // namespace detail

  template <typename Expr>
  auto evaluate_to_matrix_multi(const Expr& expr)
  {
    using mm_t = typename detail::get_terminal<Expr>::type;
    mm_t result;
    result = expr;
    return result;
  }

  // Custom Boost.YAP expression template for matrix_multi AST nodes
  template <boost::yap::expr_kind Kind, typename Tuple>
  struct matrix_multi_expr
  {
    static boost::yap::expr_kind const kind = Kind;
    Tuple elements;

    using mm_t = typename detail::get_terminal<matrix_multi_expr>::type;
    using multivector_t = mm_t;
    using scalar_t = typename mm_t::scalar_t;
    using index_set_t = typename mm_t::index_set_t;

    // Member functions
    auto isnan() const { return evaluate_to_matrix_multi(*this).isnan(); }
    auto isinf() const { return evaluate_to_matrix_multi(*this).isinf(); }
    auto even() const { return evaluate_to_matrix_multi(*this).even(); }
    auto odd() const { return evaluate_to_matrix_multi(*this).odd(); }
    auto frame() const { return evaluate_to_matrix_multi(*this).frame(); }
    auto grade() const { return evaluate_to_matrix_multi(*this).grade(); }
    auto scalar() const { return evaluate_to_matrix_multi(*this).scalar(); }
    auto norm() const { return evaluate_to_matrix_multi(*this).norm(); }
    auto max_abs() const { return evaluate_to_matrix_multi(*this).max_abs(); }
    template <typename RHS>
    auto versor(const RHS& R, const bool prechecked = false) const
    { return evaluate_to_matrix_multi(*this).versor(eval_if_expr(R), prechecked); }
    template <typename RHS>
    auto versor_exp(const RHS& A, const bool prechecked = false) const
    { return evaluate_to_matrix_multi(*this).versor_exp(eval_if_expr(A), prechecked); }
    template <typename... Args>
    decltype(auto) write(Args&&... args) const
    { return evaluate_to_matrix_multi(*this).write(std::forward<Args>(args)...); }
    auto get_matrix() const { return evaluate_to_matrix_multi(*this).get_matrix(); }
  };

  // Define binary operators plus/minus using Boost.YAP macros
  BOOST_YAP_USER_BINARY_OPERATOR(plus, matrix_multi_expr, matrix_multi_expr)
  BOOST_YAP_USER_BINARY_OPERATOR(minus, matrix_multi_expr, matrix_multi_expr)

  template <typename T>
  inline const T& eval_if_expr(const T& val)
  { return val; }

  template <boost::yap::expr_kind Kind, typename Tuple>
  inline auto eval_if_expr(const matrix_multi_expr<Kind, Tuple>& expr)
  { return evaluate_to_matrix_multi(expr); }

#define GLUCAT_EXPR_BINARY_OPERATOR(Op)                                                                             \
  template <boost::yap::expr_kind Kind, typename Tuple, typename RHS>                                               \
  inline auto operator Op(const matrix_multi_expr<Kind, Tuple>& lhs, const RHS& rhs)                                \
  {                                                                                                                 \
    return evaluate_to_matrix_multi(lhs) Op eval_if_expr(rhs);                                                      \
  }                                                                                                                 \
  template <typename LHS, boost::yap::expr_kind Kind, typename Tuple>                                               \
  inline auto operator Op(const LHS& lhs, const matrix_multi_expr<Kind, Tuple>& rhs)                                \
  {                                                                                                                 \
    return eval_if_expr(lhs) Op evaluate_to_matrix_multi(rhs);                                                      \
  }                                                                                                                 \
  template <boost::yap::expr_kind KindL, typename TupleL, boost::yap::expr_kind KindR, typename TupleR>             \
  inline auto operator Op(const matrix_multi_expr<KindL, TupleL>& lhs, const matrix_multi_expr<KindR, TupleR>& rhs) \
  {                                                                                                                 \
    return evaluate_to_matrix_multi(lhs) Op evaluate_to_matrix_multi(rhs);                                          \
  }

  GLUCAT_EXPR_BINARY_OPERATOR(*)
  GLUCAT_EXPR_BINARY_OPERATOR(/)
  GLUCAT_EXPR_BINARY_OPERATOR(^)
  GLUCAT_EXPR_BINARY_OPERATOR(&)
  GLUCAT_EXPR_BINARY_OPERATOR(%)
  GLUCAT_EXPR_BINARY_OPERATOR(|)
  GLUCAT_EXPR_BINARY_OPERATOR(==)
  GLUCAT_EXPR_BINARY_OPERATOR(!=)

#undef GLUCAT_EXPR_BINARY_OPERATOR

#define GLUCAT_EXPR_UNARY_FUNCTION(Fn)                       \
  template <boost::yap::expr_kind Kind, typename Tuple>      \
  inline auto Fn(const matrix_multi_expr<Kind, Tuple>& expr) \
  {                                                          \
    return Fn(evaluate_to_matrix_multi(expr));               \
  }

  GLUCAT_EXPR_UNARY_FUNCTION(exp)
  GLUCAT_EXPR_UNARY_FUNCTION(log)
  GLUCAT_EXPR_UNARY_FUNCTION(sqrt)
  GLUCAT_EXPR_UNARY_FUNCTION(cos)
  GLUCAT_EXPR_UNARY_FUNCTION(sin)
  GLUCAT_EXPR_UNARY_FUNCTION(cosh)
  GLUCAT_EXPR_UNARY_FUNCTION(sinh)
  GLUCAT_EXPR_UNARY_FUNCTION(tan)
  GLUCAT_EXPR_UNARY_FUNCTION(atan)
  GLUCAT_EXPR_UNARY_FUNCTION(acosh)
  GLUCAT_EXPR_UNARY_FUNCTION(asinh)
  GLUCAT_EXPR_UNARY_FUNCTION(atanh)
  GLUCAT_EXPR_UNARY_FUNCTION(complexifier)

#undef GLUCAT_EXPR_UNARY_FUNCTION

  template <boost::yap::expr_kind Kind, typename Tuple, typename Arg2, typename Arg3>
  inline auto log(const matrix_multi_expr<Kind, Tuple>& val, const Arg2& i, const Arg3& prechecked)
  { return log(evaluate_to_matrix_multi(val), eval_if_expr(i), prechecked); }

  template <template <class, int, int, class> class Multivector, typename Scalar_T, int LO, int HI, typename Tune_P,
            boost::yap::expr_kind Kind, typename Tuple>
  inline bool approx_equal(const Multivector<Scalar_T, LO, HI, Tune_P>& lhs, const matrix_multi_expr<Kind, Tuple>& rhs)
  { return approx_equal(lhs, evaluate_to_matrix_multi(rhs)); }

  template <template <class, int, int, class> class RHS, typename Scalar_T, int LO, int HI, typename Tune_P,
            boost::yap::expr_kind Kind, typename Tuple>
  inline bool approx_equal(const matrix_multi_expr<Kind, Tuple>& lhs, const RHS<Scalar_T, LO, HI, Tune_P>& rhs)
  { return approx_equal(evaluate_to_matrix_multi(lhs), rhs); }

  template <boost::yap::expr_kind KindL, typename TupleL, boost::yap::expr_kind KindR, typename TupleR>
  inline bool approx_equal(const matrix_multi_expr<KindL, TupleL>& lhs, const matrix_multi_expr<KindR, TupleR>& rhs)
  { return approx_equal(evaluate_to_matrix_multi(lhs), evaluate_to_matrix_multi(rhs)); }

  template <template <class, int, int, class> class Multivector, typename Scalar_T, int LO, int HI, typename Tune_P,
            boost::yap::expr_kind Kind, typename Tuple>
  inline bool approx_equal(const Multivector<Scalar_T, LO, HI, Tune_P>& lhs, const matrix_multi_expr<Kind, Tuple>& rhs,
                           Scalar_T threshold, Scalar_T tolerance)
  { return approx_equal(lhs, evaluate_to_matrix_multi(rhs), threshold, tolerance); }

  template <template <class, int, int, class> class RHS, typename Scalar_T, int LO, int HI, typename Tune_P,
            boost::yap::expr_kind Kind, typename Tuple>
  inline bool approx_equal(const matrix_multi_expr<Kind, Tuple>& lhs, const RHS<Scalar_T, LO, HI, Tune_P>& rhs, Scalar_T threshold,
                           Scalar_T tolerance)
  { return approx_equal(evaluate_to_matrix_multi(lhs), rhs, threshold, tolerance); }

  template <boost::yap::expr_kind KindL, typename TupleL, boost::yap::expr_kind KindR, typename TupleR, typename Scalar_T>
  inline bool approx_equal(const matrix_multi_expr<KindL, TupleL>& lhs, const matrix_multi_expr<KindR, TupleR>& rhs,
                           Scalar_T threshold, Scalar_T tolerance)
  { return approx_equal(evaluate_to_matrix_multi(lhs), evaluate_to_matrix_multi(rhs), threshold, tolerance); }

  template <boost::yap::expr_kind Kind, typename Tuple>
  inline auto norm(const matrix_multi_expr<Kind, Tuple>& expr)
  { return norm(evaluate_to_matrix_multi(expr)); }

  template <boost::yap::expr_kind Kind, typename Tuple>
  inline auto abs(const matrix_multi_expr<Kind, Tuple>& expr)
  { return abs(evaluate_to_matrix_multi(expr)); }

  // Forward declarations for friends

  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  class framed_multi;  // forward

  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  class matrix_multi;  // forward

  // Geometric product
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  matrix_multi<Scalar_T, LO, HI, Tune_P> operator*(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs,
                                                   const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs);

  // Outer product
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  matrix_multi<Scalar_T, LO, HI, Tune_P> operator^(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs,
                                                   const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs);

  // Inner product
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  matrix_multi<Scalar_T, LO, HI, Tune_P> operator&(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs,
                                                   const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs);

  // Left contraction
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  matrix_multi<Scalar_T, LO, HI, Tune_P> operator%(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs,
                                                   const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs);

  // Scalar product: [HS] (1.44) star(a, b) = scalar(a * b) = <ab>_0
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  Scalar_T star(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs, const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs);

  // Hestenes inner product: [H] (1.10) hstar(a, b) = scalar(reverse(a) * b) = <a†b>_0
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  Scalar_T hstar(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs, const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs);

  // Geometric quotient
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  matrix_multi<Scalar_T, LO, HI, Tune_P> operator/(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs,
                                                   const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs);

  // Transformation via twisted adjoint action
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  matrix_multi<Scalar_T, LO, HI, Tune_P> operator|(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs,
                                                   const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs);

  // Add two concrete matrix_multi instances, returning a Boost.YAP expression template
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  inline auto operator+(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs, const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs)
  { return boost::yap::make_terminal<matrix_multi_expr>(lhs) + boost::yap::make_terminal<matrix_multi_expr>(rhs); }

  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  inline auto operator+(matrix_multi<Scalar_T, LO, HI, Tune_P>&& lhs, const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs)
  { return boost::yap::make_terminal<matrix_multi_expr>(std::move(lhs)) + boost::yap::make_terminal<matrix_multi_expr>(rhs); }

  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  inline auto operator+(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs, matrix_multi<Scalar_T, LO, HI, Tune_P>&& rhs)
  { return boost::yap::make_terminal<matrix_multi_expr>(lhs) + boost::yap::make_terminal<matrix_multi_expr>(std::move(rhs)); }

  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  inline auto operator+(matrix_multi<Scalar_T, LO, HI, Tune_P>&& lhs, matrix_multi<Scalar_T, LO, HI, Tune_P>&& rhs)
  { return boost::yap::make_terminal<matrix_multi_expr>(std::move(lhs))
           + boost::yap::make_terminal<matrix_multi_expr>(std::move(rhs)); }

  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  inline auto operator-(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs, const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs)
  { return boost::yap::make_terminal<matrix_multi_expr>(lhs) - boost::yap::make_terminal<matrix_multi_expr>(rhs); }

  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  inline auto operator-(matrix_multi<Scalar_T, LO, HI, Tune_P>&& lhs, const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs)
  { return boost::yap::make_terminal<matrix_multi_expr>(std::move(lhs)) - boost::yap::make_terminal<matrix_multi_expr>(rhs); }

  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  inline auto operator-(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs, matrix_multi<Scalar_T, LO, HI, Tune_P>&& rhs)
  { return boost::yap::make_terminal<matrix_multi_expr>(lhs) - boost::yap::make_terminal<matrix_multi_expr>(std::move(rhs)); }

  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  inline auto operator-(matrix_multi<Scalar_T, LO, HI, Tune_P>&& lhs, matrix_multi<Scalar_T, LO, HI, Tune_P>&& rhs)
  { return boost::yap::make_terminal<matrix_multi_expr>(std::move(lhs))
           - boost::yap::make_terminal<matrix_multi_expr>(std::move(rhs)); }

  // Read multivector from input

  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  std::istream& operator>>(std::istream& s, matrix_multi<Scalar_T, LO, HI, Tune_P>& val);

  // Write multivector to output
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  std::ostream& operator<<(std::ostream& os, const matrix_multi<Scalar_T, LO, HI, Tune_P>& val);

  // Find a common frame for operands of a binary operator
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  index_set<LO, HI> reframe(const matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs, const matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs,
                            matrix_multi<Scalar_T, LO, HI, Tune_P>& lhs_reframed,
                            matrix_multi<Scalar_T, LO, HI, Tune_P>& rhs_reframed);

  // Square root of multivector with specified complexifier
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  matrix_multi<Scalar_T, LO, HI, Tune_P> sqrt(const matrix_multi<Scalar_T, LO, HI, Tune_P>& val,
                                              const matrix_multi<Scalar_T, LO, HI, Tune_P>& i, bool prechecked);

  // Square root of multivector with specified complexifier
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  matrix_multi<Scalar_T, LO, HI, Tune_P> matrix_sqrt(const matrix_multi<Scalar_T, LO, HI, Tune_P>& val,
                                                     const matrix_multi<Scalar_T, LO, HI, Tune_P>& i, const index_t level);

  // Natural logarithm of multivector with specified complexifier
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  matrix_multi<Scalar_T, LO, HI, Tune_P> log(const matrix_multi<Scalar_T, LO, HI, Tune_P>& val,
                                             const matrix_multi<Scalar_T, LO, HI, Tune_P>& i, bool prechecked);

  // Natural logarithm of multivector with specified complexifier
  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  matrix_multi<Scalar_T, LO, HI, Tune_P> matrix_log(const matrix_multi<Scalar_T, LO, HI, Tune_P>& val,
                                                    const matrix_multi<Scalar_T, LO, HI, Tune_P>& i, const index_t level);

  /// A matrix_multi<Scalar_T,LO,HI,Tune_P> is a matrix approximation to a multivector
  template <typename Scalar_T = double, const index_t LO = DEFAULT_LO, const index_t HI = DEFAULT_HI, typename Tune_P = tuning<>>
  class matrix_multi : public clifford_algebra<Scalar_T, index_set<LO, HI>, matrix_multi<Scalar_T, LO, HI, Tune_P>>
  {
  public:
    using multivector_t = matrix_multi;
    using matrix_multi_t = multivector_t;
    using scalar_t = Scalar_T;
    using tune_p = Tune_P;
    using index_set_t = index_set<LO, HI>;
    using term_t = std::pair<const index_set_t, Scalar_T>;
    using vector_t = std::vector<Scalar_T>;
    using error_t = error<multivector_t>;
    using size_type = std::size_t;
    using framed_multi_t = framed_multi<Scalar_T, LO, HI, Tune_P>;
    template <typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI, typename Other_Tune_P>
    friend class framed_multi;
    template <typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI, typename Other_Tune_P>
    friend class matrix_multi;

  public:
    using basis_matrix_t = matrix::sparse_matrix_t<int>;
    using matrix_t = matrix::matrix_t<Scalar_T>;
    using matrix_index_t = matrix::matrix_index_t;

  public:
    // Class name used in messages
    static std::string_view classname();
    /// Destructor
    ~matrix_multi() noexcept override = default;
    // Default constructor
    matrix_multi();
    // Move constructor
    matrix_multi(matrix_multi&& other) noexcept(std::is_nothrow_move_constructible_v<Scalar_T>);
    /// Default copy constructor
    matrix_multi(const matrix_multi&) = default;
    /// Construct a multivector from a multivector with a different scalar type (explicit)
    template <typename Other_Scalar_T, typename Other_Tune_P>
    explicit matrix_multi(const matrix_multi<Other_Scalar_T, LO, HI, Other_Tune_P>& val);
    /// Construct a multivector, within a given frame, from a given multivector
    template <typename Other_Scalar_T, typename Other_Tune_P>
    matrix_multi(const matrix_multi<Other_Scalar_T, LO, HI, Other_Tune_P>& val, const index_set_t frm,
                 const bool prechecked = false);
    // Construct a multivector, within a given frame, from a given multivector
    matrix_multi(const multivector_t& val, const index_set_t frm, const bool prechecked = false);
    // Construct a multivector from an index set and a scalar coordinate
    matrix_multi(const index_set_t ist, const Scalar_T& crd = Scalar_T(1));
    // Construct a multivector, within a given frame, from an index set and a scalar coordinate
    matrix_multi(const index_set_t ist, const Scalar_T& crd, const index_set_t frm, const bool prechecked = false);
    // Construct a multivector from a scalar (within a frame, if given)
    matrix_multi(const Scalar_T& scr, const index_set_t frm = index_set_t());
    // Construct a multivector from an int (within a frame, if given)
    matrix_multi(const int scr, const index_set_t frm = index_set_t());
    // Construct a multivector, within a given frame, from a given vector
    matrix_multi(const vector_t& vec, const index_set_t frm, const bool prechecked = false);
    // Construct a multivector from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
    explicit matrix_multi(const std::string& str);
    // Construct a multivector, within a given frame, from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
    matrix_multi(const std::string& str, const index_set_t frm, const bool prechecked = false);
    // Construct a multivector from a char*: eg: "3+2{1,2}-6.1e-2{2,3}"
    explicit matrix_multi(const char* str) { *this = matrix_multi(std::string(str)); };
    // Construct a multivector, within a given frame, from a char*: eg: "3+2{1,2}-6.1e-2{2,3}"
    matrix_multi(const char* str, const index_set_t frm, const bool prechecked = false)
    {
      *this = matrix_multi(std::string(str), frm, prechecked);
    };
    /// Construct a multivector from a framed_multi_t (explicit)
    template <typename Other_Scalar_T, typename Other_Tune_P>
    explicit matrix_multi(const framed_multi<Other_Scalar_T, LO, HI, Other_Tune_P>& val);
    /// Construct a multivector, within a given frame, from a framed_multi_t
    template <typename Other_Scalar_T, typename Other_Tune_P>
    matrix_multi(const framed_multi<Other_Scalar_T, LO, HI, Other_Tune_P>& val, const index_set_t frm,
                 const bool prechecked = false);
    // Constructor from a Boost.YAP expression template
    template <typename Expr>
      requires(boost::yap::is_expr<Expr>::value)
    matrix_multi(const Expr& expr);

    // Use generalized FFT to construct a matrix_multi_t
    matrix_multi_t fast_matrix_multi(const index_set_t frm) const;
    // Use inverse generalized FFT to construct a framed_multi_t
    template <typename Other_Scalar_T, typename Other_Tune_P>
    framed_multi<Other_Scalar_T, LO, HI, Other_Tune_P> fast_framed_multi() const;

    /// Project onto grade k
    matrix_multi project(int k) const;
    /// Decompose into a vector of multivectors, one for each grade 0..n
    std::vector<matrix_multi> decompose() const;

  private:
    /// Construct a multivector within a given frame from a given matrix
    template <typename Matrix_T>
      requires(!boost::yap::is_expr<Matrix_T>::value)
    matrix_multi(const Matrix_T& mtx, const index_set_t frm);
    // Construct a multivector within a given frame from a given matrix
    matrix_multi(const matrix_t& mtx, const index_set_t frm);
    // Construct a multivector within a given frame from a Boost.YAP expression
    template <typename Expr>
      requires(boost::yap::is_expr<Expr>::value)
    matrix_multi(const Expr& expr, const index_set_t frm);

  public:
    // Create a basis element matrix within the current frame
    basis_matrix_t basis_element(const index_set<LO, HI>& ist) const;

  public:
    _GLUCAT_CLIFFORD_ALGEBRA_OPERATIONS
    _GLUCAT_CLIFFORD_ALGEBRA_ASSIGNMENT_OPERATIONS

    // Number of terms
    size_type nbr_terms() const;
    // Get the underlying matrix
    const matrix_t& get_matrix() const { return m_matrix; }

  protected:
    // Number of rows
    matrix_index_t nbr_rows() const;
    // Number of columns
    matrix_index_t nbr_cols() const;

  public:
    // Move assignment
    matrix_multi& operator=(matrix_multi&& other) noexcept(std::is_nothrow_move_assignable_v<Scalar_T>);
    /// Default copy assignment
    matrix_multi& operator=(const matrix_multi&) = default;
    // Assignment from a Boost.YAP expression template
    template <typename Expr>
      requires(boost::yap::is_expr<Expr>::value)
    auto operator=(const Expr& expr) -> matrix_multi&;

    // Random multivector within a frame
    static matrix_multi_t random(const index_set_t frm, Scalar_T fill = Scalar_T(1));

    // Friend declarations

    friend matrix_multi_t operator* <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);
    friend matrix_multi_t operator^ <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);
    friend matrix_multi_t operator& <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);
    friend matrix_multi_t operator% <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);
    friend Scalar_T star<>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);
    friend matrix_multi_t operator/ <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);
    friend matrix_multi_t operator| <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);

    friend std::istream& operator>> <>(std::istream& s, multivector_t& val);
    friend std::ostream& operator<< <>(std::ostream& os, const multivector_t& val);
    // Reframe
    template <typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI, typename Other_Tune_P>
    friend index_set<Other_LO, Other_HI> reframe(const matrix_multi<Other_Scalar_T, Other_LO, Other_HI, Other_Tune_P>& lhs,
                                                 const matrix_multi<Other_Scalar_T, Other_LO, Other_HI, Other_Tune_P>& rhs,
                                                 matrix_multi<Other_Scalar_T, Other_LO, Other_HI, Other_Tune_P>& lhs_reframed,
                                                 matrix_multi<Other_Scalar_T, Other_LO, Other_HI, Other_Tune_P>& rhs_reframed);
    template <typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI, typename Other_Tune_P>
    friend matrix_multi<Other_Scalar_T, Other_LO, Other_HI, Other_Tune_P> matrix_sqrt(
        const matrix_multi<Other_Scalar_T, Other_LO, Other_HI, Other_Tune_P>& val,
        const matrix_multi<Other_Scalar_T, Other_LO, Other_HI, Other_Tune_P>& i, const index_t level);
    template <typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI, typename Other_Tune_P>
    friend matrix_multi<Other_Scalar_T, Other_LO, Other_HI, Other_Tune_P> matrix_log(
        const matrix_multi<Other_Scalar_T, Other_LO, Other_HI, Other_Tune_P>& val,
        const matrix_multi<Other_Scalar_T, Other_LO, Other_HI, Other_Tune_P>& i, const index_t level);

    // Add a term, if non-zero
    multivector_t& operator+=(const term_t& rhs);

  private:
    // Data members

    // Index set representing the frame for the subalgebra which contains the multivector
    index_set_t m_frame;
    // Matrix value representing the multivector within the folded frame
    matrix_t m_matrix;
    // Scalar part of geometric product: Trace(A*B) / Dim
    Scalar_T scalar_of_prod(const matrix_multi& other) const;
  };

  // Non-members

  template <typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P>
  matrix_multi<Scalar_T, LO, HI, Tune_P> exp(const matrix_multi<Scalar_T, LO, HI, Tune_P>& val);

}  // namespace glucat

namespace std
{
  /// Numeric limits for matrix_multi inherit limits for the corresponding scalar type
  template <typename Scalar_T, const glucat::index_t LO, const glucat::index_t HI, typename Tune_P>
  struct numeric_limits<glucat::matrix_multi<Scalar_T, LO, HI, Tune_P>> : public numeric_limits<Scalar_T>
  {
  };
}  // namespace std
#endif  // _GLUCAT_MATRIX_MULTI_H
