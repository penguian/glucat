#ifndef _GLUCAT_MATRIX_IMP_H
#define _GLUCAT_MATRIX_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_imp.h : Implement common matrix functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2007 by Paul C. Leopardi
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

namespace glucat { namespace matrix
{
  /// Kronecker tensor product of matrices - as per Matlab kron
  template< typename LHS_T, typename RHS_T >
  const
  RHS_T
  kron(const LHS_T& lhs, const RHS_T& rhs)
  {
    typedef RHS_T matrix_t;
    typedef typename matrix_t::size_type matrix_index_t;
    const matrix_index_t rhs_s1 = rhs.size1();
    const matrix_index_t rhs_s2 = rhs.size2();
    matrix_t result(lhs.size1()*rhs_s1, lhs.size2()*rhs_s2);
    result.clear();

    typedef typename matrix_t::value_type Scalar_T;
    typedef typename LHS_T::const_iterator1 lhs_const_iterator1;
    typedef typename LHS_T::const_iterator2 lhs_const_iterator2;
    typedef typename RHS_T::const_iterator1 rhs_const_iterator1;
    typedef typename RHS_T::const_iterator2 rhs_const_iterator2;
    for (lhs_const_iterator1
        lhs_it1 = lhs.begin1();
        lhs_it1 != lhs.end1();
        ++lhs_it1)
      for (lhs_const_iterator2
          lhs_it2 = lhs_it1.begin();
          lhs_it2 != lhs_it1.end();
          ++lhs_it2)
      {
        const matrix_index_t start1 = rhs_s1 * lhs_it2.index1();
        const matrix_index_t start2 = rhs_s2 * lhs_it2.index2();
        const Scalar_T& lhs_val = *lhs_it2;
        for (rhs_const_iterator1
            rhs_it1 = rhs.begin1();
            rhs_it1 != rhs.end1();
            ++rhs_it1)
          for (rhs_const_iterator2
              rhs_it2 = rhs_it1.begin();
              rhs_it2 != rhs_it1.end();
              ++rhs_it2)
            result(start1 + rhs_it2.index1(), start2 + rhs_it2.index2()) = lhs_val * *rhs_it2;
      }
    return result;
  }

  /// Sparse Kronecker tensor product of monomial matrices
  template< typename LHS_T, typename RHS_T >
  const
  RHS_T
  mono_kron(const LHS_T& lhs, const RHS_T& rhs)
  {
    typedef RHS_T matrix_t;
    typedef typename matrix_t::size_type matrix_index_t;
    const matrix_index_t rhs_s1 = rhs.size1();
    const matrix_index_t rhs_s2 = rhs.size2();
    const matrix_index_t dim = lhs.size1()*rhs_s1;
    matrix_t result(dim, dim, dim);
    result.clear();

    typedef typename matrix_t::value_type Scalar_T;
    typedef typename LHS_T::const_iterator1 lhs_const_iterator1;
    typedef typename LHS_T::const_iterator2 lhs_const_iterator2;
    typedef typename RHS_T::const_iterator1 rhs_const_iterator1;
    typedef typename RHS_T::const_iterator2 rhs_const_iterator2;
    for (lhs_const_iterator1
        lhs_it1 = lhs.begin1();
        lhs_it1 != lhs.end1();
        ++lhs_it1)
    {
      const lhs_const_iterator2 lhs_it2 = lhs_it1.begin();
      const matrix_index_t start1 = rhs_s1 * lhs_it2.index1();
      const matrix_index_t start2 = rhs_s2 * lhs_it2.index2();
      const Scalar_T& lhs_val = *lhs_it2;
      for (rhs_const_iterator1
          rhs_it1 = rhs.begin1();
          rhs_it1 != rhs.end1();
          ++rhs_it1)
      {
        const rhs_const_iterator2 rhs_it2 = rhs_it1.begin();
        result(start1 + rhs_it2.index1(), start2 + rhs_it2.index2()) = lhs_val * *rhs_it2;
      }
    }
    return result;
  }

  /// Left inverse of Kronecker product
  template< typename LHS_T, typename RHS_T >
  const
  RHS_T
  nork(const LHS_T& lhs, const RHS_T& rhs, const bool mono)
  {
    // nork(A, kron(A, B)) is close to B
    typedef RHS_T matrix_t;
    typedef typename LHS_T::size_type lhs_index_t;
    typedef typename matrix_t::size_type matrix_index_t;
    const lhs_index_t lhs_s1 = lhs.size1();
    const lhs_index_t lhs_s2 = lhs.size2();
    const matrix_index_t rhs_s1 = rhs.size1();
    const matrix_index_t rhs_s2 = rhs.size2();
    const matrix_index_t res_s1 = rhs_s1 / lhs_s1;
    const matrix_index_t res_s2 = rhs_s2 / lhs_s2;
    typedef typename matrix_t::value_type Scalar_T;
    const Scalar_T& nnz_lhs = Scalar_T( mono ? lhs_s1 : nnz(lhs) );
    if (!mono)
    {
      typedef error<matrix_t> error_t;
      if (rhs_s1 == 0)
        throw error_t("matrix", "nork: number of rows must not be 0");
      if (rhs_s2 == 0)
        throw error_t("matrix", "nork: number of cols must not be 0");
      if (res_s1 * lhs_s1 != rhs_s1)
        throw error_t("matrix", "nork: incompatible numbers of rows");
      if (res_s2 * lhs_s2 != rhs_s2)
        throw error_t("matrix", "nork: incompatible numbers of cols");
      if (nnz_lhs == Scalar_T(0))
        throw error_t("matrix", "nork: LHS must not be 0");
    }
    matrix_t result(res_s1, res_s2);
    result.clear();
    for (typename LHS_T::const_iterator1
        lhs_it1 = lhs.begin1();
        lhs_it1 != lhs.end1();
        ++lhs_it1)
      for (typename LHS_T::const_iterator2
          lhs_it2 = lhs_it1.begin();
          lhs_it2 != lhs_it1.end();
          ++lhs_it2)
        if (*lhs_it2 != Scalar_T(0))
        {
          typedef ublas::matrix_range<const matrix_t> matrix_range_t;
          const Scalar_T& lhs_scale = *lhs_it2 * nnz_lhs;
          const matrix_index_t start1 = res_s1 * lhs_it2.index1();
          const matrix_index_t start2 = res_s2 * lhs_it2.index2();
          using ublas::range;
          const range& range1 = range(start1, start1 + res_s1);
          const range& range2 = range(start2, start2 + res_s2);
          const matrix_range_t& rhs_range = matrix_range_t(rhs, range1, range2);
          for (typename matrix_range_t::const_iterator1
              rhs_it1 = rhs_range.begin1();
              rhs_it1 != rhs_range.end1();
              ++rhs_it1)
            for (typename matrix_range_t::const_iterator2
                rhs_it2 = rhs_it1.begin();
                rhs_it2 != rhs_it1.end();
                ++rhs_it2)
                result(rhs_it2.index1(), rhs_it2.index2()) += *rhs_it2 / lhs_scale;
        }
    return result;
  }

  /// Number of non-zeros
  template< typename Matrix_T >
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

  /// Not a Number
  template< typename Matrix_T >
  bool
  isnan(const Matrix_T& m)
  {
    typedef Matrix_T matrix_t;
    typedef typename matrix_t::value_type Scalar_T;
    typedef typename matrix_t::size_type matrix_index_t;
    typedef typename matrix_t::const_iterator1 const_iterator1;
    typedef typename matrix_t::const_iterator2 const_iterator2;
    for (const_iterator1
        it1 = m.begin1();
        it1 != m.end1();
        ++it1)
      for (const_iterator2
          it2 = it1.begin();
          it2 != it1.end();
          ++it2)
        if (numeric_traits<Scalar_T>::isNaN(*it2))
          return true;

    return false;
  }

  /// Unit matrix - as per Matlab eye
  template< typename Matrix_T >
  inline
  const
  Matrix_T
  unit(const typename Matrix_T::size_type dim)
  {
    typedef typename Matrix_T::value_type Scalar_T;
    return ublas::identity_matrix<Scalar_T>(dim);
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
    typedef typename matrix_t::value_type Scalar_T;
    typedef typename lhs_expression_t::const_iterator1   lhs_const_iterator1;
    typedef typename lhs_expression_t::const_iterator2   lhs_const_iterator2;
    typedef typename ublas::matrix_row<rhs_expression_t> matrix_row_t;
    typedef typename matrix_row_t::const_iterator        row_const_iterator;

    const matrix_index_t dim = lhs().size1();
    // The following assumes that matrix_t is a sparse matrix type.
    matrix_t result = matrix_t(dim, dim, dim);
    for (lhs_const_iterator1
        lhs_row = lhs().begin1();
        lhs_row != lhs().end1();
        ++lhs_row)
    {
      const lhs_const_iterator2& lhs_it = lhs_row.begin();
      if (lhs_it != lhs_row.end())
      {
        const matrix_row_t& rhs_row = matrix_row_t(rhs(), lhs_it.index2());
        const row_const_iterator& rhs_it = rhs_row.begin();
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
    typedef typename RHS_T::expression_type expression_t;
    return ublas::sparse_prod<expression_t>(lhs(), rhs());
  }

  /// Product of matrices
  template< typename LHS_T, typename RHS_T >
  inline
  const typename RHS_T::expression_type
  prod(const ublas::matrix_expression<LHS_T>& lhs,
       const ublas::matrix_expression<RHS_T>& rhs)
  {
#ifdef _GLUCAT_USE_DENSE_MATRICES
    typedef typename RHS_T::size_type matrix_index_t;
    const matrix_index_t dim = lhs().size1();
    RHS_T result(dim, dim);
    ublas::axpy_prod(lhs, rhs, result, true);
    return result;
#else
    typedef typename RHS_T::expression_type expression_t;
    return ublas::sparse_prod<expression_t>(lhs(), rhs());
#endif
  }

  /// Inner product: sum(lhs(i,j)*rhs(i,j))/lhs.nrows()
  template< typename Scalar_T, typename LHS_T, typename RHS_T >
  Scalar_T
  inner(const LHS_T& lhs, const RHS_T& rhs)
  {
    Scalar_T result = Scalar_T(0);
    for (typename LHS_T::const_iterator1
        lhs_it1 = lhs.begin1();
        lhs_it1 != lhs.end1();
        ++lhs_it1)
      for (typename LHS_T::const_iterator2
          lhs_it2 = lhs_it1.begin();
          lhs_it2 != lhs_it1.end();
          ++lhs_it2)
      {
        const Scalar_T& rhs_val = rhs(lhs_it2.index1(),lhs_it2.index2());
        if (rhs_val != Scalar_T(0))
          result += (*lhs_it2) * rhs_val;
      }
    return result / lhs.size1();
  }

  /// Square of Frobenius norm
  template< typename Matrix_T >
  typename Matrix_T::value_type
  norm_frob2(const Matrix_T& val)
  {
    typedef typename Matrix_T::value_type Scalar_T;
    typedef typename Matrix_T::size_type matrix_index_t;

    Scalar_T result = Scalar_T(0);
    for (typename Matrix_T::const_iterator1
        val_it1 = val.begin1();
        val_it1 != val.end1();
        ++val_it1)
      for (typename Matrix_T::const_iterator2
          val_it2 = val_it1.begin();
          val_it2 != val_it1.end();
          ++val_it2)
      {
        if (numeric_traits<Scalar_T>::isNaN(*val_it2))
          return numeric_traits<Scalar_T>::NaN();
        result += (*val_it2) * (*val_it2);
      } 
    return result;
  }

  /// Matrix trace
  template< typename Matrix_T >
  typename Matrix_T::value_type
  trace(const Matrix_T& val)
  {
    typedef typename Matrix_T::value_type Scalar_T;
    typedef typename Matrix_T::size_type matrix_index_t;

    Scalar_T result = Scalar_T(0);
    matrix_index_t dim = val.size1();
    for (matrix_index_t ndx=0;
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
} }

#endif  // _GLUCAT_MATRIX_IMP_H
