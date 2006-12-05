#ifndef _GLUCAT_MATRIX_IMP_H
#define _GLUCAT_MATRIX_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_imp.h : Implement common matrix functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001 by Paul C. Leopardi
                         : uBLAS interface contributed by Joerg Walter
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

    for (typename Matrix_T::const_iterator1 
        i = lhs.begin1(); 
        i != lhs.end1(); 
        ++i)
      for (typename Matrix_T::const_iterator2 
          j = i.begin(); 
          j != i.end(); 
          ++j)
      {
        const matrix_index_t rj1 = j.index1()*rhs_nrows;
        const matrix_index_t cj2 = j.index2()*rhs_ncols;
        const scalar_t lhs_ij = (*j);
        for (typename Matrix_T::const_iterator1 
            k = rhs.begin1(); 
            k != rhs.end1(); 
            ++k)
          for (typename Matrix_T::const_iterator2 
              l = k.begin(); 
              l != k.end(); 
              ++l)
            result(rj1 + l.index1(), cj2 + l.index2()) = lhs_ij * (*l);
      }
    return result;
  }

  /// Left inverse of Kronecker product
  template< typename Matrix_T >
  const
  Matrix_T
  nork(const Matrix_T& lhs, const Matrix_T& rhs, const bool mono)
  {
    // nork(A, kron(A, B)) is close to B
    typedef Matrix_T matrix_t;
    typedef typename matrix_t::size_type matrix_index_t;
    const matrix_index_t lhs_s1 = lhs.size1();
    const matrix_index_t lhs_s2 = lhs.size2();
    const matrix_index_t rhs_s1 = rhs.size1();
    const matrix_index_t rhs_s2 = rhs.size2();
    const matrix_index_t res_s1 = rhs_s1 / lhs_s1;
    const matrix_index_t res_s2 = rhs_s2 / lhs_s2;
    typedef error<matrix_t> error_t;
    if (rhs_s1 == 0)
      throw error_t("matrix", "nork: number of rows must not be 0");
    if (rhs_s2 == 0)
      throw error_t("matrix", "nork: number of cols must not be 0");
    if (res_s1 * lhs_s1 != rhs_s1)
      throw error_t("matrix", "nork: incompatible numbers of rows");
    if (res_s2 * lhs_s2 != rhs_s2)
      throw error_t("matrix", "nork: incompatible numbers of cols");
    typedef typename matrix_t::value_type scalar_t;
    const scalar_t nnz_lhs = scalar_t( mono ? lhs_s1 : nnz(lhs) );
    if (nnz_lhs == scalar_t(0))
      throw error_t("matrix", "nork: LHS must not be 0");
    matrix_t result(res_s1, res_s2);
    typedef typename matrix_t::const_iterator1 const_iterator1;
    typedef typename matrix_t::const_iterator2 const_iterator2;
    for (const_iterator1 
        lhs_it1 = lhs.begin1();
        lhs_it1 != lhs.end1(); 
        ++lhs_it1)
      for (const_iterator2 
          lhs_it2 = lhs_it1.begin();
          lhs_it2 != lhs_it1.end(); 
          ++lhs_it2)
        if (*lhs_it2 != scalar_t(0))
        {
          using namespace ublas;
          typedef matrix_range<const matrix_t> matrix_range_t;
          const matrix_index_t start1 = res_s1 * lhs_it2.index1();
          const matrix_index_t start2 = res_s2 * lhs_it2.index2();
          const range& range1 = range(start1, start1+res_s1);
          const range& range2 = range(start2, start2+res_s2);
          const matrix_range_t& rhs_range = matrix_range_t(rhs, range1, range2);
          result += rhs_range / (*lhs_it2 * nnz_lhs);
        }
    return result;
  }

  /// Number of non-zeros
  template< typename Matrix_T >
  inline
  typename Matrix_T::size_type
  nnz(const Matrix_T& m)
  {
    typedef Matrix_T matrix_t;
    typedef typename matrix_t::size_type matrix_index_t;
    typedef typename matrix_t::const_iterator1 const_iterator1;
    typedef typename matrix_t::const_iterator2 const_iterator2;
    matrix_index_t result = 0;
    for (const_iterator1 
        it1 = m.begin1(); 
        it1 != m.end1(); 
        ++it1)
      for (const_iterator2 
          it2 = it1.begin(); 
          it2 != it1.end(); 
          ++it2)
        if (*it2 != 0)
          ++result;
    return result;
  }

  /// Unit matrix - as per Matlab eye
  template< typename Matrix_T >
  inline
  const
  Matrix_T
  unit(const typename Matrix_T::size_type dim)
  {
    typedef typename Matrix_T::value_type scalar_t;
    return ublas::identity_matrix<scalar_t>(dim);
  }

  /// Product of monomial matrices
  template< typename LHS_T, typename RHS_T >
  const typename RHS_T::expression_type
  mono_prod(const ublas::matrix_expression<LHS_T>& lhs,
            const ublas::matrix_expression<RHS_T>& rhs)
  {
    typedef const LHS_T lhs_expression_t;
    typedef const RHS_T rhs_expression_t;
    typedef typename RHS_T::expression_type matrix_t;
    typedef typename matrix_t::size_type  matrix_index_t;
    typedef typename matrix_t::value_type scalar_t;
    typedef typename lhs_expression_t::const_iterator1   lhs_const_iterator1;
    typedef typename lhs_expression_t::const_iterator2   lhs_const_iterator2;
    typedef typename ublas::matrix_row<rhs_expression_t> matrix_row_t;
    typedef typename matrix_row_t::const_iterator        row_const_iterator;

    const matrix_index_t dim = lhs().size1();
    matrix_t result(dim, dim);
    for (lhs_const_iterator1 
        lhs_row = lhs().begin1(); 
        lhs_row != lhs().end1(); 
        ++lhs_row)
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

  /// Product of sparse matrices
  template< typename LHS_T, typename RHS_T >
  inline
  const typename RHS_T::expression_type
  sparse_prod(const ublas::matrix_expression<LHS_T>& lhs,
              const ublas::matrix_expression<RHS_T>& rhs)
  {
    typedef typename RHS_T::expression_type matrix_t;
    return ublas::sparse_prod<matrix_t>(lhs(), rhs());
  }

  /// Inner product: sum(lhs(i,j)*rhs(i,j))/lhs.nrows()
  template< typename Scalar_T, typename Matrix_T >
  Scalar_T
  inner(const Matrix_T& lhs, const Matrix_T& rhs)
  {
    Scalar_T result = 0;
    for (typename Matrix_T::const_iterator1 
        i = lhs.begin1(); 
        i != lhs.end1(); 
        ++i)
      for (typename Matrix_T::const_iterator2 
          j = i.begin(); 
          j != i.end(); 
          ++j)
      {
        const Scalar_T rhs_j12 = rhs(j.index1(),j.index2());
        if (rhs_j12 != Scalar_T(0))
          result += (*j) * rhs_j12;
      }
    return result / lhs.size1();
  }
} }

#endif  // _GLUCAT_MATRIX_IMP_H
