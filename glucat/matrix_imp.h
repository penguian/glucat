#ifndef _GLUCAT_MATRIX_IMP_H
#define _GLUCAT_MATRIX_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_imp.h : Implement common matrix functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001 by Paul C. Leopardi
    email                : leopardi@bigpond.net.au
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

namespace glucat { namespace matrix
{
  /// Kronecker tensor product of matrices - as per Matlab kron
  template< typename Matrix_T >
  const
  Matrix_T
  kron(const Matrix_T& lhs, const Matrix_T& rhs)
  {
    typedef typename Matrix_T::size_type  matrix_index_t;
    typedef typename Matrix_T::value_type scalar_t;
    const matrix_index_t lhs_nrows = lhs.size1();
    const matrix_index_t lhs_ncols = lhs.size2();
    const matrix_index_t rhs_nrows = rhs.size1();
    const matrix_index_t rhs_ncols = rhs.size2();
    Matrix_T result(lhs_nrows*rhs_nrows, lhs_ncols*rhs_ncols);

    for (      typename Matrix_T::const_iterator1 i = lhs.begin1(); i != lhs.end1(); ++i)
      for (    typename Matrix_T::const_iterator2 j = i.begin(); j != i.end(); ++j)
      {
        const matrix_index_t rj1 = j.index1()*rhs_nrows;
        const matrix_index_t cj2 = j.index2()*rhs_ncols;
        const scalar_t lhs_ij = (*j);
        for (  typename Matrix_T::const_iterator1 k = rhs.begin1(); k != rhs.end1(); ++k)
          for (typename Matrix_T::const_iterator2 l = k.begin(); l != k.end(); ++l)
              result(rj1 + l.index1(), cj2 + l.index2()) = lhs_ij * (*l);
      }
    return result;
  }

  /// Unit matrix - as per Matlab eye
  template< typename Matrix_T >
  const
  Matrix_T
  unit(const typename Matrix_T::size_type dim)
  {
    typedef typename Matrix_T::size_type matrix_index_t;
    Matrix_T result(dim, dim);
    for (matrix_index_t k = 0; k != dim; ++k)
      result(k, k) = 1;
    return result;
  }

  /// Equality of matrices
  template< typename LHS_T, typename RHS_T >
  inline
  bool
  operator==(const ublas::matrix_expression<LHS_T>& lhs,
             const ublas::matrix_expression<RHS_T>& rhs)
  { return ublas::equals(lhs, rhs); }

  /// Product of monomial matrices
  template< typename Matrix_T, typename LHS_T, typename RHS_T >
  const Matrix_T
  mono_prod(const ublas::matrix_expression<LHS_T>& lhs,
            const ublas::matrix_expression<RHS_T>& rhs)
  {
    typedef const LHS_T lhs_expression_t;
    typedef const RHS_T rhs_expression_t;
    typedef typename Matrix_T::size_type  matrix_index_t;
    typedef typename Matrix_T::value_type scalar_t;
    typedef typename lhs_expression_t::const_iterator1   lhs_const_iterator1;
    typedef typename lhs_expression_t::const_iterator2   lhs_const_iterator2;
    typedef typename ublas::matrix_row<rhs_expression_t> matrix_row_t;
    typedef typename matrix_row_t::const_iterator        row_const_iterator;

    const matrix_index_t dim = lhs().size1();
    Matrix_T result(dim, dim);
    for (lhs_const_iterator1 lhs_row = lhs().begin1(); lhs_row != lhs().end1(); ++lhs_row)
    {
      const lhs_const_iterator2& lhs_it = lhs_row.begin();
      if (lhs_it != lhs_row.end())
      {
        const matrix_row_t rhs_row(rhs(), lhs_it.index2());
        const row_const_iterator&  rhs_it = rhs_row.begin();
        if (rhs_it != rhs_row.end())
          result(lhs_it.index1(), rhs_it.index()) = (*lhs_it) * (*rhs_it);
      }
    }
    return result;
  }

  /// Product of compressed matrices
  template< typename Matrix_T, typename LHS_T, typename RHS_T >
  inline
  const Matrix_T
  compressed_prod(const ublas::matrix_expression<LHS_T>& lhs,
                  const ublas::matrix_expression<RHS_T>& rhs)
  { return ublas::sparse_prod<Matrix_T>(lhs, rhs, typename Matrix_T::orientation_category()); }

  /// Inner product: sum(lhs(i,j)*rhs(i,j))/lhs.nrows()
  template< typename Matrix_T, typename Scalar_T >
  Scalar_T
  inner(const Matrix_T& lhs, const Matrix_T& rhs)
  {
    Scalar_T result = 0;
    for (typename Matrix_T::const_iterator1 i = lhs.begin1(); i != lhs.end1(); ++i)
      for (typename Matrix_T::const_iterator2 j = i.begin(); j != i.end(); ++j)
      {
        const Scalar_T rhs_j12 = rhs(j.index1(),j.index2());
        if (rhs_j12 != Scalar_T(0))
          result += (*j) * rhs_j12;
      }
    return result / lhs.size1();
  }
} }
#endif  // _GLUCAT_MATRIX_IMP_H
