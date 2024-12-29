#ifndef _GLUCAT_MATRIX_IMP_H
#define _GLUCAT_MATRIX_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_imp.h : Implement common matrix functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2012 by Paul C. Leopardi
                         : uBLAS interface contributed by Joerg Walter
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

#include "glucat/errors.h"
#include "glucat/scalar.h"
#include "glucat/matrix.h"

# if  defined(_GLUCAT_GCC_IGNORE_UNUSED_LOCAL_TYPEDEFS)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-local-typedefs"
# endif
# if  defined(_GLUCAT_HAVE_BOOST_SERIALIZATION_ARRAY_WRAPPER_H)
#  include <boost/serialization/array_wrapper.hpp>
# endif
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

#if defined(_GLUCAT_USE_BLAZE)
#include <blaze/Math.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DynamicVector.h>
#endif

# if defined(_GLUCAT_GCC_IGNORE_UNUSED_LOCAL_TYPEDEFS)
#  pragma GCC diagnostic pop
# endif

#include <set>
#include <vector>

namespace glucat { namespace matrix
{
  // References for algorithms:
  // [v]: C. F. van Loan and N. Pitsianis, "Approximation with Kronecker products",
  // in Linear Algebra for Large Scale and Real-Time Applications, Marc S. Moonen,
  // Gene H. Golub, and Bart L. R. Moor (eds.), 1993, pp. 293--314.

  /// Kronecker tensor product of matrices - as per Matlab kron
  template< typename LHS_T, typename RHS_T >
  auto
  kron(const LHS_T& lhs, const RHS_T& rhs) -> const
  RHS_T
  {
    const auto rhs_s1 = rhs.size1();
    const auto rhs_s2 = rhs.size2();
    auto result = RHS_T(lhs.size1()*rhs_s1, lhs.size2()*rhs_s2);
    result.clear();

    for (auto
        lhs_it1 = lhs.begin1();
        lhs_it1 != lhs.end1();
        ++lhs_it1)
      for (auto
          lhs_it2 = lhs_it1.begin();
          lhs_it2 != lhs_it1.end();
          ++lhs_it2)
      {
        const auto start1 = rhs_s1 * lhs_it2.index1();
        const auto start2 = rhs_s2 * lhs_it2.index2();
        const auto& lhs_val = *lhs_it2;
        for (auto
            rhs_it1 = rhs.begin1();
            rhs_it1 != rhs.end1();
            ++rhs_it1)
          for (auto
              rhs_it2 = rhs_it1.begin();
              rhs_it2 != rhs_it1.end();
              ++rhs_it2)
            result(start1 + rhs_it2.index1(), start2 + rhs_it2.index2()) = lhs_val * *rhs_it2;
      }
    return result;
  }

  /// Sparse Kronecker tensor product of monomial matrices
  template< typename LHS_T, typename RHS_T >
  auto
  mono_kron(const LHS_T& lhs, const RHS_T& rhs) -> const
  RHS_T
  {
    const auto rhs_s1 = rhs.size1();
    const auto rhs_s2 = rhs.size2();
    const auto dim = lhs.size1()*rhs_s1;
    auto result = RHS_T(dim, dim, dim);
    result.clear();

    for (auto
        lhs_it1 = lhs.begin1();
        lhs_it1 != lhs.end1();
        ++lhs_it1)
    {
      const auto lhs_it2 = lhs_it1.begin();
      const auto start1 = rhs_s1 * lhs_it2.index1();
      const auto start2 = rhs_s2 * lhs_it2.index2();
      const auto& lhs_val = *lhs_it2;
      for (auto
          rhs_it1 = rhs.begin1();
          rhs_it1 != rhs.end1();
          ++rhs_it1)
      {
        const auto rhs_it2 = rhs_it1.begin();
        result(start1 + rhs_it2.index1(), start2 + rhs_it2.index2()) = lhs_val * *rhs_it2;
      }
    }
    return result;
  }

  /// Utility routine for nork: calculate result for a range of indices
  template< typename LHS_T, typename RHS_T >
  void
  nork_range(RHS_T& result,
             const typename LHS_T::const_iterator2 lhs_it2,
             const RHS_T& rhs,
             const typename RHS_T::size_type res_s1,
             const typename RHS_T::size_type res_s2)
  {
    // Definition matches [v] Section 4, Theorem 4.1.
    const auto start1 = res_s1 * lhs_it2.index1();
    const auto start2 = res_s2 * lhs_it2.index2();
    using ublas::range;
    const auto& range1 = range(start1, start1 + res_s1);
    const auto& range2 = range(start2, start2 + res_s2);
    using matrix_range_t = ublas::matrix_range<const RHS_T>;
    const auto& rhs_range = matrix_range_t(rhs, range1, range2);
    using Scalar_T = typename RHS_T::value_type;
    const auto lhs_val = numeric_traits<Scalar_T>::to_scalar_t(*lhs_it2);
    for (auto
        rhs_it1 = rhs_range.begin1();
        rhs_it1 != rhs_range.end1();
        ++rhs_it1)
      for (auto
          rhs_it2 = rhs_it1.begin();
          rhs_it2 != rhs_it1.end();
          ++rhs_it2)
          result(rhs_it2.index1(), rhs_it2.index2()) += lhs_val * *rhs_it2;
   }

  /// Left inverse of Kronecker product
  template< typename LHS_T, typename RHS_T >
  auto
  nork(const LHS_T& lhs, const RHS_T& rhs, const bool mono) -> const
  RHS_T
  {
    // nork(A, kron(A, B)) is close to B
    // Definition matches [v] Section 4, Theorem 4.1.
    const auto lhs_s1 = lhs.size1();
    const auto lhs_s2 = lhs.size2();
    const auto rhs_s1 = rhs.size1();
    const auto rhs_s2 = rhs.size2();
    const auto res_s1 = rhs_s1 / lhs_s1;
    const auto res_s2 = rhs_s2 / lhs_s2;
    using Scalar_T = typename RHS_T::value_type;
    const auto norm_frob2_lhs = norm_frob2(lhs);
    if (!mono)
    {
      using error_t = error<RHS_T>;
      if (rhs_s1 == 0)
        throw error_t("matrix", "nork: number of rows must not be 0");
      if (rhs_s2 == 0)
        throw error_t("matrix", "nork: number of cols must not be 0");
      if (res_s1 * lhs_s1 != rhs_s1)
        throw error_t("matrix", "nork: incompatible numbers of rows");
      if (res_s2 * lhs_s2 != rhs_s2)
        throw error_t("matrix", "nork: incompatible numbers of cols");
      if (norm_frob2_lhs == Scalar_T(0))
        throw error_t("matrix", "nork: LHS must not be 0");
    }
    auto result = RHS_T(res_s1, res_s2);
    result.clear();
    for (auto
        lhs_it1 = lhs.begin1();
        lhs_it1 != lhs.end1();
        ++lhs_it1)
      for (auto
          lhs_it2 = lhs_it1.begin();
          lhs_it2 != lhs_it1.end();
          ++lhs_it2)
        if (*lhs_it2  != Scalar_T(0))
          nork_range<LHS_T, RHS_T>(result, lhs_it2, rhs, res_s1, res_s2);
    result /= norm_frob2_lhs;
    return result;
  }

  /// Left inverse of Kronecker product where lhs is a signed permutation matrix
  template< typename LHS_T, typename RHS_T >
  auto
  signed_perm_nork(const LHS_T& lhs, const RHS_T& rhs) -> const
  RHS_T
  {
    // signed_perm_nork(A, kron(A, B)) is close to B
    // Definition matches [v] Section 4, Theorem 4.1.
    const auto lhs_s1 = lhs.size1();
    const auto lhs_s2 = lhs.size2();
    const auto rhs_s1 = rhs.size1();
    const auto rhs_s2 = rhs.size2();
    const auto res_s1 = rhs_s1 / lhs_s1;
    const auto res_s2 = rhs_s2 / lhs_s2;
    using Scalar_T = typename RHS_T::value_type;
    const auto norm_frob2_lhs = Scalar_T( double(lhs_s1) );
    auto result = RHS_T(res_s1, res_s2);
    result.clear();
    for (auto
        lhs_it1 = lhs.begin1();
        lhs_it1 != lhs.end1();
        ++lhs_it1)
    {
      const auto lhs_it2 = lhs_it1.begin();
      nork_range<LHS_T, RHS_T>(result, lhs_it2, rhs, res_s1, res_s2);
    }
    result /= norm_frob2_lhs;
    return result;
  }

  /// Number of non-zeros
  template< typename Matrix_T >
  auto
  nnz(const Matrix_T& m) -> typename Matrix_T::size_type
  {
    using size_t = typename Matrix_T::size_type;
    auto result = size_t(0);
    for (auto
        it1 = m.begin1();
        it1 != m.end1();
        ++it1)
      for (auto& entry : it1)
        if (entry != 0)
          ++result;
    return result;
  }

  /// Infinite
  template< typename Matrix_T >
  auto
  isinf(const Matrix_T& m) -> bool
  {
    using Scalar_T = typename Matrix_T::value_type;
    for (auto
        it1 = m.begin1();
        it1 != m.end1();
        ++it1)
      for (auto& entry : it1)
        if (numeric_traits<Scalar_T>::isInf(entry))
          return true;

    return false;
  }

  /// Not a Number
  template< typename Matrix_T >
  auto
  isnan(const Matrix_T& m) -> bool
  {
    using Scalar_T = typename Matrix_T::value_type;
    for (auto
        it1 = m.begin1();
        it1 != m.end1();
        ++it1)
      for (auto& entry : it1)
        if (numeric_traits<Scalar_T>::isNaN(entry))
          return true;

    return false;
  }

  /// Unit matrix - as per Matlab eye
  template< typename Matrix_T >
  inline
  auto
  unit(const typename Matrix_T::size_type dim) -> const
  Matrix_T
  {
    using Scalar_T = typename Matrix_T::value_type;
    return ublas::identity_matrix<Scalar_T>(dim);
  }

  /// Product of monomial matrices
  template< typename LHS_T, typename RHS_T >
  auto
  mono_prod(const ublas::matrix_expression<LHS_T>& lhs,
            const ublas::matrix_expression<RHS_T>& rhs) -> const typename RHS_T::expression_type
  {
    using rhs_expression_t = const RHS_T;
    using matrix_row_t = typename ublas::matrix_row<rhs_expression_t>;

    const auto dim = lhs().size1();
    // The following assumes that RHS_T is a sparse matrix type.
    auto result = RHS_T(dim, dim, dim);
    for (auto
        lhs_row = lhs().begin1();
        lhs_row != lhs().end1();
        ++lhs_row)
    {
      const auto& lhs_it = lhs_row.begin();
      if (lhs_it != lhs_row.end())
      {
        const auto& rhs_row = matrix_row_t(rhs(), lhs_it.index2());
        const auto& rhs_it = rhs_row.begin();
        if (rhs_it != rhs_row.end())
          result(lhs_it.index1(), rhs_it.index()) = (*lhs_it) * (*rhs_it);
      }
    }
    return result;
  }

  /// Product of sparse matrices
  template< typename LHS_T, typename RHS_T >
  inline
  auto
  sparse_prod(const ublas::matrix_expression<LHS_T>& lhs,
              const ublas::matrix_expression<RHS_T>& rhs) -> const typename RHS_T::expression_type
  {
    using expression_t = typename RHS_T::expression_type;
    return ublas::sparse_prod<expression_t>(lhs(), rhs());
  }

  /// Product of matrices
  template< typename LHS_T, typename RHS_T >
  inline
  auto
  prod(const ublas::matrix_expression<LHS_T>& lhs,
       const ublas::matrix_expression<RHS_T>& rhs) -> const typename RHS_T::expression_type
  {
    const auto dim = lhs().size1();
    RHS_T result(dim, dim);
    ublas::axpy_prod(lhs, rhs, result, true);
    return result;
  }

  /// Inner product: sum(lhs(i,j)*rhs(i,j))/lhs.nrows()
  template< typename Scalar_T, typename LHS_T, typename RHS_T >
  auto
  inner(const LHS_T& lhs, const RHS_T& rhs) -> Scalar_T
  {
    auto result = Scalar_T(0);
    for (auto
        lhs_it1 = lhs.begin1();
        lhs_it1 != lhs.end1();
        ++lhs_it1)
      for (auto
          lhs_it2 = lhs_it1.begin();
          lhs_it2 != lhs_it1.end();
          ++lhs_it2)
      {
        const auto& rhs_val = rhs(lhs_it2.index1(),lhs_it2.index2());
        if (rhs_val != Scalar_T(0))
          result += (*lhs_it2) * rhs_val;
      }
    return result / lhs.size1();
  }

  /// Square of Frobenius norm
  template< typename Matrix_T >
  auto
  norm_frob2(const Matrix_T& val) -> typename Matrix_T::value_type
  {
    using Scalar_T = typename Matrix_T::value_type;

    auto result = Scalar_T(0);
    for (auto
        val_it1 = val.begin1();
        val_it1 != val.end1();
        ++val_it1)
      for (auto& val_entry : val_it1)
      {
        if (numeric_traits<Scalar_T>::isNaN(val_entry))
          return numeric_traits<Scalar_T>::NaN();
        result += val_entry * val_entry;
      }
    return result;
  }

  /// Matrix trace
  template< typename Matrix_T >
  auto
  trace(const Matrix_T& val) -> typename Matrix_T::value_type
  {
    using Scalar_T = typename Matrix_T::value_type;

    auto result = Scalar_T(0);
    auto dim = val.size1();
    for (auto
        ndx = decltype(dim)(0);
        ndx != dim;
        ++ndx)
    {
      const Scalar_T crd = val(ndx, ndx);
      if (numeric_traits<Scalar_T>::isNaN(crd))
        return numeric_traits<Scalar_T>::NaN();
      result += crd;
    }
    return result;
  }

#if defined(_GLUCAT_USE_BLAZE)
  /// Convert matrix to Blaze format
  template< typename Matrix_T >
  static
  auto
  to_blaze(const Matrix_T& val) -> blaze::DynamicMatrix<double,blaze::rowMajor>
  {
    const auto s1 = val.size1();
    const auto s2 = val.size2();

    using blaze_matrix_t = typename blaze::DynamicMatrix<double,blaze::rowMajor>;
    auto result = blaze_matrix_t(s1, s2);

    using Scalar_T = typename Matrix_T::value_type;
    using traits_t = numeric_traits<Scalar_T>;

    for (auto
        val_it1 = val.begin1();
        val_it1 != val.end1();
        ++val_it1)
      for (auto
          val_it2 = val_it1.begin();
          val_it2 != val_it1.end();
          ++val_it2)
        result(val_it2.index1(), val_it2.index2()) = traits_t::to_double(*val_it2);

    return result;
  }

#endif

  /// Eigenvalues of a matrix
  template< typename Matrix_T >
  auto
  eigenvalues(const Matrix_T& val) -> std::vector< std::complex<double> >
  {
    using complex_t = std::complex<double>;
    using complex_vector_t = typename std::vector<complex_t>;

    const auto dim = val.size1();
    auto lambda = complex_vector_t(dim);

#if defined(_GLUCAT_USE_BLAZE)
    using blaze_matrix_t = typename blaze::DynamicMatrix<double, blaze::rowMajor>;
    using complex_t = std::complex<double>;
    using blaze_complex_vector_t = blaze::DynamicVector<complex_t, blaze::columnVector>;

    auto blaze_val = to_blaze(val);
    auto blaze_lambda = blaze_complex_vector_t(dim);
    blaze::geev(blaze_val, blaze_lambda);

    for (auto
        k = decltype(dim)(0);
        k != dim;
        ++k)
      lambda[k] = blaze_lambda[k];
#endif
    return lambda;
  }

  /// Classify the eigenvalues of a matrix
  template< typename Matrix_T >
  auto
  classify_eigenvalues(const Matrix_T& val) -> eig_genus<Matrix_T>
  {
    using Scalar_T = typename Matrix_T::value_type;
    eig_genus<Matrix_T> result;

    using complex_t = std::complex<double>;
    using complex_vector_t = typename std::vector<complex_t>;
    auto lambda = eigenvalues(val);

    std::set<double> arg_set;

    using vector_index_t = typename complex_vector_t::size_type;
    const auto dim = lambda.size();
    static const auto epsilon =
      std::max(std::numeric_limits<double>::epsilon(),
               numeric_traits<Scalar_T>::to_double(std::numeric_limits<Scalar_T>::epsilon()));
    static const auto zero_eig_tol = 4096.0*epsilon;

    bool neg_real_eig_found = false;
    bool imag_eig_found = false;
    bool zero_eig_found = false;

    for (auto
        k = decltype(dim)(0);
        k != dim;
        ++k)
    {
      const auto lambda_k = lambda[k];
      arg_set.insert(std::arg(lambda_k));

      const auto real_lambda_k = std::real(lambda_k);
      const auto imag_lambda_k = std::imag(lambda_k);
      const auto norm_tol = 4096.0*epsilon*std::norm(lambda_k);

      if (!neg_real_eig_found &&
          real_lambda_k < -epsilon &&
         (imag_lambda_k == 0.0 ||
          imag_lambda_k * imag_lambda_k < norm_tol))
        neg_real_eig_found = true;
      if (!imag_eig_found &&
          imag_lambda_k > epsilon &&
         (real_lambda_k == 0.0 ||
          real_lambda_k * real_lambda_k < norm_tol))
        imag_eig_found = true;
      if (!zero_eig_found &&
        std::norm(lambda_k) < zero_eig_tol)
        zero_eig_found = true;
    }

    if (zero_eig_found)
        result.m_is_singular = true;

    static const auto pi = numeric_traits<double>::pi();
    if (neg_real_eig_found)
    {
      if (imag_eig_found)
        result.m_eig_case = both_eigs;
      else
      {
        result.m_eig_case = neg_real_eigs;
        result.m_safe_arg = Scalar_T(-pi / 2.0);
      }
    }

    if (result.m_eig_case == both_eigs)
    {
      auto arg_it = arg_set.begin();
      auto first_arg = *arg_it;
      auto best_arg = first_arg;
      auto best_diff = 0.0;
      auto previous_arg = first_arg;
      for (++arg_it;
          arg_it != arg_set.end();
          ++arg_it)
      {
        const auto arg_diff = *arg_it - previous_arg;
        if (arg_diff > best_diff)
        {
          best_diff = arg_diff;
          best_arg = previous_arg;
        }
        previous_arg = *arg_it;
      }
      const auto arg_diff = first_arg + 2.0*pi - previous_arg;
      if (arg_diff > best_diff)
      {
        best_diff = arg_diff;
        best_arg = previous_arg;
      }
      result.m_safe_arg = Scalar_T(pi - (best_arg + best_diff / 2.0));
    }
    return result;
  }
} }

#endif  // _GLUCAT_MATRIX_IMP_H
