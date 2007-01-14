#ifndef _GLUCAT_MATRIX_H
#define _GLUCAT_MATRIX_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix.h : Declare common matrix functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2007 by Paul C. Leopardi
                         : uBLAS interface contributed by Joerg Walter
 ***************************************************************************
 *   This library is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2.1 of the  *
 *   License, or (at your option) any later version.                       *
 *   See http://www.fsf.org/copyleft/lesser.html for details               *
 ***************************************************************************
 This library is based on a prototype written by Arvind Raja and was
 licensed under the LGPL with permission of the author. See Arvind Raja,
 "Object-oriented implementations of Clifford algebras in C++: a prototype",
 in Ablamowicz, Lounesto and Parra (eds.)
 "Clifford algebras with numeric and symbolic computations", Birkhauser, 1996.
 ***************************************************************************
     See also Arvind Raja's original header comments in glucat.h
 ***************************************************************************/

namespace glucat
{
  namespace ublas = boost::numeric::ublas;

  namespace matrix
  {
    /// Kronecker tensor product of matrices - as per Matlab kron
    template< typename Matrix_T >
    const
    Matrix_T
    kron(const Matrix_T& lhs, const Matrix_T& rhs);

    /// Left inverse of Kronecker product
    template< typename Matrix_T >
    const
    Matrix_T
    nork(const Matrix_T& lhs, const Matrix_T& rhs, const bool mono = true);

    /// Number of non-zeros
    template< typename Matrix_T >
    typename Matrix_T::size_type
    nnz(const Matrix_T& m);

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

    /// Inner product: sum(x(i,j)*y(i,j))/x.nrows()
    template< typename Scalar_T, typename Matrix_T >
    Scalar_T
    inner(const Matrix_T& lhs, const Matrix_T& rhs);
  }
}

#endif  // _GLUCAT_MATRIX_H
