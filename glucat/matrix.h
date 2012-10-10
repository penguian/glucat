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

namespace glucat
{
  namespace ublas = boost::numeric::ublas;

  namespace matrix
  {
    /// Kronecker tensor product of matrices - as per Matlab kron
    template< typename LHS_T, typename RHS_T >
    const
    RHS_T
    kron(const LHS_T& lhs, const RHS_T& rhs);

    /// Sparse Kronecker tensor product of monomial matrices
    template< typename LHS_T, typename RHS_T >
    const
    RHS_T
    mono_kron(const LHS_T& lhs, const RHS_T& rhs);

    /// Left inverse of Kronecker product
    template< typename LHS_T, typename RHS_T >
    const
    RHS_T
    nork(const LHS_T& lhs, const RHS_T& rhs, const bool mono = true);

    /// Number of non-zeros
    template< typename Matrix_T >
    typename Matrix_T::size_type
    nnz(const Matrix_T& m);

    /// Not a Number
    template< typename Matrix_T >
    bool
    isnan(const Matrix_T& m);

    /// Unit matrix - as per Matlab eye
    template< typename Matrix_T >
    const
    Matrix_T
    unit(const typename Matrix_T::size_type n);

    /// Product of monomial matrices
    template< typename LHS_T, typename RHS_T >
    const
    typename RHS_T::expression_type
    mono_prod(const ublas::matrix_expression<LHS_T>& lhs,
              const ublas::matrix_expression<RHS_T>& rhs);

    /// Product of sparse matrices
    template< typename LHS_T, typename RHS_T >
    const
    typename RHS_T::expression_type
    sparse_prod(const ublas::matrix_expression<LHS_T>& lhs,
                const ublas::matrix_expression<RHS_T>& rhs);

    /// Product of matrices
    template< typename LHS_T, typename RHS_T >
    const
    typename RHS_T::expression_type
    prod(const ublas::matrix_expression<LHS_T>& lhs,
         const ublas::matrix_expression<RHS_T>& rhs);

    /// Inner product: sum(x(i,j)*y(i,j))/x.nrows()
    template< typename Scalar_T, typename LHS_T, typename RHS_T >
    Scalar_T
    inner(const LHS_T& lhs, const RHS_T& rhs);

    /// Square of Frobenius norm
    template< typename Matrix_T >
    typename Matrix_T::value_type
    norm_frob2(const Matrix_T& val);

    /// Matrix trace
    template< typename Matrix_T >
    typename Matrix_T::value_type
    trace(const Matrix_T& val);

    /// Eigenvalues of a matrix
    template< typename Matrix_T >
    ublas::vector< std::complex<double> >
    eigenvalues(const Matrix_T& val);

    /// Classification of eigenvalues of a matrix
    typedef enum {safe_eig_case, negative_eig_case, both_eig_case} eig_case_t;

    ///  Structure containing classification of eigenvalues
    template< typename Matrix_T >
    struct eig_genus
    {
      typedef typename Matrix_T::value_type Scalar_T;
      /// What kind of eigenvalues does the matrix contain?
      eig_case_t m_eig_case;
      /// Argument such that exp(pi-m_safe_arg) lies between arguments of eigenvalues
      Scalar_T   m_safe_arg;
    };

    /// Classify the eigenvalues of a matrix
    template< typename Matrix_T >
    eig_genus<Matrix_T>
    classify_eigenvalues(const Matrix_T& val);
  }
}

#endif  // _GLUCAT_MATRIX_H
