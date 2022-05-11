#ifndef _GLUCAT_MATRIX_H
#define _GLUCAT_MATRIX_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix.h : Declare common matrix functions
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

#include <boost/numeric/ublas/fwd.hpp>

#include <complex>
#include <vector>

namespace glucat
{
  namespace ublas = boost::numeric::ublas;

  namespace matrix
  {
    /// Kronecker tensor product of matrices - as per Matlab kron
    template< typename LHS_T, typename RHS_T >
    auto
    kron(const LHS_T& lhs, const RHS_T& rhs) -> const
    RHS_T;

    /// Sparse Kronecker tensor product of monomial matrices
    template< typename LHS_T, typename RHS_T >
    auto
    mono_kron(const LHS_T& lhs, const RHS_T& rhs) -> const
    RHS_T;

    /// Left inverse of Kronecker product
    template< typename LHS_T, typename RHS_T >
    auto
    nork(const LHS_T& lhs, const RHS_T& rhs, const bool mono = true) -> const
    RHS_T;

    /// Left inverse of Kronecker product where lhs is a signed permutation matrix
    template< typename LHS_T, typename RHS_T >
    auto
    signed_perm_nork(const LHS_T& lhs, const RHS_T& rhs) -> const
    RHS_T;

    /// Number of non-zeros
    template< typename Matrix_T >
    auto
    nnz(const Matrix_T& m) -> typename Matrix_T::size_type;

    /// Infinite
    template< typename Matrix_T >
    auto
    isinf(const Matrix_T& m) -> bool;

    /// Not a Number
    template< typename Matrix_T >
    auto
    isnan(const Matrix_T& m) -> bool;

    /// Unit matrix - as per Matlab eye
    template< typename Matrix_T >
    auto
    unit(const typename Matrix_T::size_type n) -> const
    Matrix_T;

    /// Product of monomial matrices
    template< typename LHS_T, typename RHS_T >
    auto
    mono_prod(const ublas::matrix_expression<LHS_T>& lhs,
              const ublas::matrix_expression<RHS_T>& rhs) -> const
    typename RHS_T::expression_type;

    /// Product of sparse matrices
    template< typename LHS_T, typename RHS_T >
    auto
    sparse_prod(const ublas::matrix_expression<LHS_T>& lhs,
                const ublas::matrix_expression<RHS_T>& rhs) -> const
    typename RHS_T::expression_type;

    /// Product of matrices
    template< typename LHS_T, typename RHS_T >
    auto
    prod(const ublas::matrix_expression<LHS_T>& lhs,
         const ublas::matrix_expression<RHS_T>& rhs) -> const
    typename RHS_T::expression_type;

    /// Inner product: sum(x(i,j)*y(i,j))/x.nrows()
    template< typename Scalar_T, typename LHS_T, typename RHS_T >
    auto
    inner(const LHS_T& lhs, const RHS_T& rhs) -> Scalar_T;

    /// Square of Frobenius norm
    template< typename Matrix_T >
    auto
    norm_frob2(const Matrix_T& val) -> typename Matrix_T::value_type;

    /// Matrix trace
    template< typename Matrix_T >
    auto
    trace(const Matrix_T& val) -> typename Matrix_T::value_type;

    /// Eigenvalues of a matrix
    template< typename Matrix_T >
    auto
    eigenvalues(const Matrix_T& val) -> std::vector< std::complex<double> >;

    /// Classification of eigenvalues of a matrix
    using eig_case_t = enum {
      safe_eigs,
      neg_real_eigs,
      both_eigs};

    ///  Structure containing classification of eigenvalues
    template< typename Matrix_T >
    struct eig_genus
    {
      using Scalar_T = typename Matrix_T::value_type;
      /// Is the matrix singular?
      bool m_is_singular = false;
      /// What kind of eigenvalues does the matrix contain?
      eig_case_t m_eig_case = safe_eigs;
      /// Argument such that exp(pi-m_safe_arg) lies between arguments of eigenvalues
      Scalar_T   m_safe_arg = Scalar_T(0);
    };

    /// Classify the eigenvalues of a matrix
    template< typename Matrix_T >
    auto
    classify_eigenvalues(const Matrix_T& val) -> eig_genus<Matrix_T>;
  }
}

#endif  // _GLUCAT_MATRIX_H
