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

  /// Product of sparse matrices
  template< typename LHS_T, typename RHS_T >
  inline
  const typename RHS_T::expression_type
  sparse_prod(const ublas::matrix_expression<LHS_T>& lhs,
              const ublas::matrix_expression<RHS_T>& rhs)
  {
    typedef typename RHS_T::expression_type matrix_t;
    typedef typename matrix_t::orientation_category orientation_category;
    return ublas::sparse_prod<matrix_t>(lhs, rhs, orientation_category());
  }

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

  // uBLAS interface by Joerg Walter
  /// Swap rows of a vector
  template<class PM, class MV>
  BOOST_UBLAS_INLINE
  void
  swap_rows (const PM &pm, MV &mv, ublas::vector_tag)
  {
    typedef typename PM::size_type size_type;
    typedef typename MV::value_type value_type;

    size_type size = pm.size ();
    for (size_type i = 0; i < size; ++ i)
    {
      value_type t (mv (i));
      mv (i) = mv (pm (i));
      mv (pm (i)) = t;
    }
  }

  /// Swap rows of a matrix
  template<class PM, class MV>
  BOOST_UBLAS_INLINE
  void
  swap_rows (const PM &pm, MV &mv, ublas::matrix_tag)
  {
    typedef typename PM::size_type size_type;
    typedef typename MV::value_type value_type;

    size_type size = pm.size ();
    for (size_type i = 0; i < size; ++ i)
      for (size_type j = 0; j < mv.size2 (); ++ j)
      {
        value_type t (mv (i, j));
        mv (i, j) = mv (pm (i), j);
        mv (pm (i), j) = t;
      }
  }

  /// Swap rows of a matrix or a vector: dispatcher
  template<class PM, class MV>
  BOOST_UBLAS_INLINE
  void
  swap_rows (const PM &pm, MV &mv)
  { swap_rows (pm, mv, BOOST_UBLAS_TYPENAME MV::type_category ()); }

  /// LU factorize: return LU in place of matrix
  template<class M, class PM>
  typename M::size_type lu_factorize (M &m, PM &pm)
  {
    using namespace ublas;
    typedef M matrix_type;
    typedef BOOST_UBLAS_TYPENAME M::size_type size_type;
    typedef BOOST_UBLAS_TYPENAME M::value_type value_type;

#ifdef BOOST_UBLAS_TYPE_CHECK
    matrix_type cm (m);
#endif
    int singular = 0;
    size_type size1 = m.size1 ();
    size_type size2 = m.size2 ();
    size_type size = std::min (size1, size2);
    for (size_type i = 0; i < size; ++ i) {
      matrix_column<M> mci (column (m, i));
      matrix_row<M> mri (row (m, i));
      size_type i_norm_inf = i + index_norm_inf (project (mci, range (i, size1)));
      BOOST_UBLAS_CHECK (i_norm_inf < m.size1 (), external_logic ());
      if (m (i_norm_inf, i) != value_type ())
      {
        pm (i) = i_norm_inf;
        if (i_norm_inf != i)
            row (m, i_norm_inf).swap (mri);
        project (mci, range (i + 1, size1)) *= value_type (1) / m (i, i);
      } else if (singular == 0)
        singular = i + 1;
      project (m, range (i + 1, size1), range (i + 1, size2)).minus_assign (
        outer_prod (project (mci, range (i + 1, size1)),
                    project (mri, range (i + 1, size2))));
    }
#ifdef BOOST_UBLAS_TYPE_CHECK
    swap_rows (pm, cm);
    BOOST_UBLAS_CHECK (
      singular != 0 ||
      equals (prod (triangular_adaptor<matrix_type, unit_lower> (m),
                    triangular_adaptor<matrix_type, upper> (m)), cm), internal_logic ());
#endif
    return singular;
  }

  /// LU solve: linear solve using LU and permutation
  template<class M, class PM, class MV>
  void
  lu_solve (M &m, const PM &pm, MV &mv)
  {
    using namespace ublas;
    typedef M matrix_type;
    typedef MV matrix_vector_type;

    swap_rows (pm, mv);
#ifdef BOOST_UBLAS_TYPE_CHECK
    matrix_vector_type cmv1 (mv);
#endif
    inplace_solve (triangular_adaptor<matrix_type, unit_lower> (m),
                   mv, unit_lower_tag (),
                   BOOST_UBLAS_TYPENAME MV::type_category ());
#ifdef BOOST_UBLAS_TYPE_CHECK
    BOOST_UBLAS_CHECK (equals
      (prod(triangular_adaptor<matrix_type, unit_lower> (m), mv), cmv1), internal_logic ());
    matrix_vector_type cmv2 (mv);
#endif
    inplace_solve (triangular_adaptor<matrix_type, upper> (m),
                   mv, upper_tag (),
                   BOOST_UBLAS_TYPENAME MV::type_category ());
#ifdef BOOST_UBLAS_TYPE_CHECK
    BOOST_UBLAS_CHECK (equals
      (prod (triangular_adaptor<matrix_type, upper> (m), mv), cmv2), internal_logic ());
#endif
  }

} }
#endif  // _GLUCAT_MATRIX_IMP_H
