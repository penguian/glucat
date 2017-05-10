#ifndef _GLUCAT_MATRIX_MULTI_IMP_H
#define _GLUCAT_MATRIX_MULTI_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_multi_imp.h : Implement the matrix representation of a multivector
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

#include "glucat/matrix_multi.h"

#include "glucat/matrix.h"
#include "glucat/generation.h"

# if  defined(_GLUCAT_GCC_IGNORE_UNUSED_LOCAL_TYPEDEFS)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-local-typedefs"
# endif
# if  defined(_GLUCAT_HAVE_BOOST_SERIALIZATION_ARRAY_WRAPPER_H)
#  include <boost/serialization/array_wrapper.hpp>
# endif
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
# if defined(_GLUCAT_GCC_IGNORE_UNUSED_LOCAL_TYPEDEFS)
#  pragma GCC diagnostic pop
# endif

#include <fstream>
#include <iomanip>

namespace glucat
{
  // References for algorithms:
  // [M]: Scott Meyers, "Effective C++" Second Edition, Addison-Wesley, 1998.
  // [P]: Ian R. Porteous, "Clifford algebras and the classical groups", Cambridge UP, 1995.
  // [L]: Pertti Lounesto, "Clifford algebras and spinors", Cambridge UP, 1997.

  /// Class name used in messages
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const std::string
  matrix_multi<Scalar_T,LO,HI>::
  classname()
  { return "matrix_multi"; }

  /// Determine the log2 dim corresponding to signature p, q
  // Reference: [P] Table 15.27, p 133
  inline
  index_t
  offset_level(const index_t p, const index_t q)
  {
    // Offsets between the log2 of the matrix dimension for the current signature
    // and that of the real superalgebra
    static const int offset_log2_dim[] = {0, 1, 0, 1, 1, 2, 1, 1};
    const index_t bott = pos_mod(p-q, 8);
    return (p+q)/2 + offset_log2_dim[bott];
  }

  /// Determine the matrix dimension of the fold of a subalegbra
  // Reference: [P] Table 15.27, p 133
  template< typename Matrix_Index_T, const index_t LO, const index_t HI >
  inline
  static
  Matrix_Index_T
  folded_dim( const index_set<LO,HI>& sub )
  { return 1 << offset_level(sub.count_pos(), sub.count_neg()); }

  /// Default constructor
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi()
  : m_frame( index_set_t() ),
    m_matrix( matrix_t( 1, 1 ) )
  { this->m_matrix.clear(); }

  /// Construct a multivector from a multivector with a different scalar type
  template< typename Scalar_T, const index_t LO, const index_t HI >
  template< typename Other_Scalar_T >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const matrix_multi<Other_Scalar_T,LO,HI>& val)
  : m_frame( val.m_frame ), m_matrix( val.m_matrix.size1(), val.m_matrix.size2() )
  {
    this->m_matrix.clear();
    typedef typename matrix_multi<Other_Scalar_T,LO,HI>::matrix_t other_matrix_t;
    typedef typename other_matrix_t::const_iterator1 other_const_iterator1;
    typedef typename other_matrix_t::const_iterator2 other_const_iterator2;
    for (other_const_iterator1
        val_it1 = val.m_matrix.begin1();
        val_it1 != val.m_matrix.end1();
        ++val_it1)
      for (other_const_iterator2
          val_it2 = val_it1.begin();
          val_it2 != val_it1.end();
          ++val_it2)
        this->m_matrix(val_it2.index1(), val_it2.index2()) = numeric_traits<Scalar_T>::to_scalar_t(*val_it2);
  }

  /// Construct a multivector, within a given frame, from a given multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  template< typename Other_Scalar_T >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const matrix_multi<Other_Scalar_T,LO,HI>& val, const index_set_t frm, const bool prechecked)
  : m_frame( frm )
  {
    if (frm != val.m_frame)
      *this = multivector_t(framed_multi_t(val), frm);
    else
    {
      const matrix_index_t dim = folded_dim<matrix_index_t>(frm);
      this->m_matrix.resize(dim, dim, false);
      this->m_matrix.clear();
      typedef typename matrix_multi<Other_Scalar_T,LO,HI>::matrix_t other_matrix_t;
      typedef typename other_matrix_t::const_iterator1 other_const_iterator1;
      typedef typename other_matrix_t::const_iterator2 other_const_iterator2;
      for (other_const_iterator1
          val_it1 = val.m_matrix.begin1();
          val_it1 != val.m_matrix.end1();
          ++val_it1)
        for (other_const_iterator2
            val_it2 = val_it1.begin();
            val_it2 != val_it1.end();
            ++val_it2)
          this->m_matrix(val_it2.index1(), val_it2.index2()) = numeric_traits<Scalar_T>::to_scalar_t(*val_it2);
    }
  }

  /// Construct a multivector, within a given frame, from a given multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const multivector_t& val, const index_set_t frm, const bool prechecked)
  : m_frame( frm )
  {
    if (frm != val.m_frame)
      *this = multivector_t(framed_multi_t(val), frm);
    else
      this->m_matrix = val.m_matrix;
  }

  /// Construct a multivector from an index set and a scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const index_set_t ist, const Scalar_T& crd)
  : m_frame( ist )
  {
    const matrix_index_t dim = folded_dim<matrix_index_t>(this->m_frame);
    this->m_matrix.resize(dim, dim, false);
    this->m_matrix.clear();
    *this += term_t(ist, crd);
  }

  /// Construct a multivector, within a given frame, from an index set and a scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const index_set_t ist, const Scalar_T& crd, const index_set_t frm, const bool prechecked)
  : m_frame( frm )
  {
    if (!prechecked && (ist | frm) != frm)
      throw error_t("multivector_t(ist,crd,frm): cannot initialize with value outside of frame");
    const matrix_index_t dim = folded_dim<matrix_index_t>(frm);
    this->m_matrix.resize(dim, dim, false);
    this->m_matrix.clear();
    *this += term_t(ist, crd);
  }

  /// Construct a multivector from a scalar (within a frame, if given)
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const Scalar_T& scr, const index_set_t frm)
  : m_frame( frm )
  {
    const matrix_index_t dim = folded_dim<matrix_index_t>(frm);
    this->m_matrix.resize(dim, dim, false);
    this->m_matrix.clear();
    *this += term_t(index_set_t(), scr);
  }

  /// Construct a multivector from an int (within a frame, if given)
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const int scr, const index_set_t frm)
  { *this = multivector_t(Scalar_T(scr), frm); }

  /// Construct a multivector, within a given frame, from a given vector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const vector_t& vec,
               const index_set_t frm, const bool prechecked)
  : m_frame( frm )
  {
    if (!prechecked && index_t(vec.size()) != frm.count())
      throw error_t("multivector_t(vec,frm): cannot initialize with vector not matching frame");
    const matrix_index_t dim = folded_dim<matrix_index_t>(frm);
    this->m_matrix.resize(dim, dim, false);
    this->m_matrix.clear();
    typename vector_t::const_iterator vec_it = vec.begin();
    const index_t begin_index = frm.min();
    const index_t end_index = frm.max()+1;

    for (index_t
        idx = begin_index;
        idx != end_index;
        ++idx)
      if (frm[idx])
      {
        *this += term_t(index_set_t(idx), *vec_it);
        ++vec_it;
      }
  }

  /// Construct a multivector from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const std::string& str)
  { *this = framed_multi_t(str); }

  /// Construct a multivector, within a given frame, from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const std::string& str, const index_set_t frm, const bool prechecked)
  { *this = multivector_t(framed_multi_t(str), frm, prechecked); }

  /// Construct a multivector from a framed_multi_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  template< typename Other_Scalar_T >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const framed_multi<Other_Scalar_T,LO,HI>& val)
  : m_frame( val.frame() )
  {
    if (val.size() >= Tune_P::fast_size_threshold)
      try
      {
        *this = val.template fast_matrix_multi<Scalar_T>(this->m_frame);
        return;
      }
      catch (const glucat_error& e)
      { }
    const matrix_index_t dim = folded_dim<matrix_index_t>(this->m_frame);
    this->m_matrix.resize(dim, dim, false);
    this->m_matrix.clear();

    typedef framed_multi<Other_Scalar_T,LO,HI> framed_multi_t;
    for (typename framed_multi_t::const_iterator
        val_it = val.begin();
        val_it != val.end();
        ++val_it)
      *this += *val_it;
  }

  /// Construct a multivector, within a given frame, from a framed_multi_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  template< typename Other_Scalar_T >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const framed_multi<Other_Scalar_T,LO,HI>& val, const index_set_t frm, const bool prechecked)
  {
    const index_set_t our_frame = val.frame() | frm;
    if (val.size() >= Tune_P::fast_size_threshold)
      try
      {
        *this = val.template fast_matrix_multi<Scalar_T>(our_frame);
        return;
      }
      catch (const glucat_error& e)
      { }
    this->m_frame = our_frame;
    const matrix_index_t dim = folded_dim<matrix_index_t>(our_frame);
    this->m_matrix.resize(dim, dim, false);
    this->m_matrix.clear();

    typedef framed_multi<Other_Scalar_T,LO,HI> framed_multi_t;
    for (typename framed_multi_t::const_iterator
        val_it = val.begin();
        val_it != val.end();
        ++val_it)
      *this += *val_it;
  }

  /// Construct a multivector within a given frame from a given matrix
  template< typename Scalar_T, const index_t LO, const index_t HI >
  template< typename Matrix_T >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const Matrix_T& mtx, const index_set_t frm)
  : m_frame( frm ), m_matrix( mtx.size1(), mtx.size2() )
  {
    this->m_matrix.clear();

    typedef typename Matrix_T::const_iterator1 const_iterator1;
    typedef typename Matrix_T::const_iterator2 const_iterator2;
    for (const_iterator1
        mtx_it1 = mtx.begin1();
        mtx_it1 != mtx.end1();
        ++mtx_it1)
      for (const_iterator2
          mtx_it2 = mtx_it1.begin();
          mtx_it2 != mtx_it1.end();
          ++mtx_it2)
        this->m_matrix(mtx_it2.index1(), mtx_it2.index2()) = numeric_traits<Scalar_T>::to_scalar_t(*mtx_it2);
  }

  /// Construct a multivector within a given frame from a given matrix
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const matrix_t& mtx, const index_set_t frm)
  : m_frame( frm ), m_matrix( mtx )
  { }

  /// Assignment operator
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator= (const matrix_multi_t& rhs)
  {
    // Check for assignment to self
    if (this == &rhs)
      return *this;
    this->m_frame = rhs.m_frame;
    this->m_matrix = rhs.m_matrix;
    return *this;
  }

  /// Find a common frame for operands of a binary operator
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const index_set<LO,HI>
  reframe (const matrix_multi<Scalar_T,LO,HI>& lhs,    const matrix_multi<Scalar_T,LO,HI>& rhs,
                 matrix_multi<Scalar_T,LO,HI>& lhs_reframed, matrix_multi<Scalar_T,LO,HI>& rhs_reframed)
  {
    typedef index_set<LO,HI> index_set_t;
    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    typedef typename multivector_t::framed_multi_t framed_multi_t;
    // Determine the initial common frame
    index_set_t our_frame = lhs.m_frame | rhs.m_frame;
    framed_multi_t framed_lhs;
    framed_multi_t framed_rhs;
    if ((lhs.m_frame != our_frame) || (rhs.m_frame != our_frame))
    {
      // The common frame may expand as a result of the transform to framed_multi_t
      framed_lhs = framed_multi_t(lhs);
      framed_rhs = framed_multi_t(rhs);
      our_frame |= framed_lhs.frame() | framed_rhs.frame();
    }
    // Do the reframing only where necessary
    if (lhs.m_frame != our_frame)
      lhs_reframed = multivector_t(framed_lhs, our_frame, true);
    if (rhs.m_frame != our_frame)
      rhs_reframed = multivector_t(framed_rhs, our_frame, true);
    return our_frame;
  }

  /// Test for equality of multivectors
  template< typename Scalar_T, const index_t LO, const index_t HI >
  bool
  matrix_multi<Scalar_T,LO,HI>::
  operator==  (const multivector_t& rhs) const
  {
    // Ensure that there is no aliasing
    if (this == &rhs)
      return true;

    // Operate only within a common frame
    multivector_t lhs_reframed;
    multivector_t rhs_reframed;
    const index_set_t our_frame = reframe(*this, rhs, lhs_reframed, rhs_reframed);
    const multivector_t& lhs_ref = (this->m_frame == our_frame)
      ? *this
      : lhs_reframed;
    const multivector_t& rhs_ref = (rhs.m_frame == our_frame)
      ? rhs
      : rhs_reframed;

#if defined(_GLUCAT_USE_DENSE_MATRICES)
    return ublas::norm_inf(lhs_ref.m_matrix - rhs_ref.m_matrix) == 0;
#else
    typedef typename matrix_t::const_iterator1 const_iterator1;
    typedef typename matrix_t::const_iterator2 const_iterator2;
    // If either matrix contains zero entries,
    // compare using subtraction and ublas::norm_inf
    for (const_iterator1
        it1 =  lhs_ref.m_matrix.begin1();
        it1 != lhs_ref.m_matrix.end1();
        ++it1)
      for (const_iterator2
          it2 =  it1.begin();
          it2 != it1.end();
          ++it2)
        if (*it2 == 0)
          return ublas::norm_inf(lhs_ref.m_matrix - rhs_ref.m_matrix) == 0;
    for (const_iterator1
        it1 =  rhs_ref.m_matrix.begin1();
        it1 != rhs_ref.m_matrix.end1();
        ++it1)
      for (const_iterator2
          it2 =  it1.begin();
          it2 != it1.end();
          ++it2)
        if (*it2 == 0)
          return ublas::norm_inf(lhs_ref.m_matrix - rhs_ref.m_matrix) == 0;
    // Neither matrix contains zero entries.
    // Compare by iterating over both matrices in lock step.
    const_iterator1 this_it1 = lhs_ref.m_matrix.begin1();
    const_iterator1 it1 =      rhs_ref.m_matrix.begin1();
    for (;
        (this_it1 != lhs_ref.m_matrix.end1()) && (it1 != rhs_ref.m_matrix.end1());
        ++this_it1, ++it1)
    {
      if ( this_it1.index1() != it1.index1() )
        return false;
      const_iterator2 this_it2 = this_it1.begin();
      const_iterator2 it2 = it1.begin();
      for (;
          (this_it2 != this_it1.end()) && (it2 != it1.end());
          ++this_it2, ++it2)
        if ( (this_it2.index2() != it2.index2()) || (*this_it2 != *it2) )
          return false;
      if ( (this_it2 != this_it1.end()) || (it2 != it1.end()) )
        return false;
    }
    return (this_it1 == lhs_ref.m_matrix.end1()) && (it1 == rhs_ref.m_matrix.end1());
#endif
  }

  // Test for equality of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  bool
  matrix_multi<Scalar_T,LO,HI>::
  operator==  (const Scalar_T& scr) const
  {
    if (scr != Scalar_T(0))
      return *this == multivector_t(framed_multi_t(scr), this->m_frame, true);
    else if (ublas::norm_inf(this->m_matrix) != 0)
      return false;
    else
    {
      const matrix_index_t dim = this->m_matrix.size1();
      return !(dim == 1 && this->isnan());
    }
  }

  /// Geometric sum of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator+= (const Scalar_T& scr)
  { return *this += term_t(index_set_t(), scr); }

  /// Geometric sum
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator+= (const multivector_t& rhs)
  {
    // Ensure that there is no aliasing
    if (this == &rhs)
      return *this *= Scalar_T(2);

    // Operate only within a common frame
    multivector_t rhs_reframed;
    const index_set_t our_frame = reframe(*this, rhs, *this, rhs_reframed);
    const multivector_t& rhs_ref = (rhs.m_frame == our_frame)
      ? rhs
      : rhs_reframed;

    noalias(this->m_matrix) += rhs_ref.m_matrix;
    return *this;
  }

  /// Geometric difference
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator-= (const multivector_t& rhs)
  {
    // Ensure that there is no aliasing
    if (this == &rhs)
      return *this = Scalar_T(0);

    // Operate only within a common frame
    multivector_t rhs_reframed;
    const index_set_t our_frame = reframe(*this, rhs, *this, rhs_reframed);
    const multivector_t& rhs_ref = (rhs.m_frame == our_frame)
      ? rhs
      : rhs_reframed;

    noalias(this->m_matrix) -= rhs_ref.m_matrix;
    return *this;
  }

  /// Unary -
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  operator- () const
  { return multivector_t(-(this->m_matrix), this->m_frame); }

  /// Product of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator*= (const Scalar_T& scr)
  { // multiply coordinates of all terms by scalar

    typedef numeric_traits<Scalar_T> traits_t;
    if (traits_t::isNaN_or_isInf(scr) || this->isnan())
      return *this = traits_t::NaN();
    if (scr == Scalar_T(0))
      *this = Scalar_T(0);
    else
      this->m_matrix *= scr;
    return *this;
  }

  /// Geometric product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  operator* (const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs)
  {
    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    typedef typename multivector_t::index_set_t index_set_t;

#if defined(_GLUCAT_CHECK_ISNAN)
    if (lhs.isnan() || rhs.isnan())
      return numeric_traits<Scalar_T>::NaN();
#endif

    // Operate only within a common frame
    multivector_t lhs_reframed;
    multivector_t rhs_reframed;
    const index_set_t our_frame = reframe(lhs, rhs, lhs_reframed, rhs_reframed);
    const multivector_t& lhs_ref = (lhs.m_frame == our_frame)
      ? lhs
      : lhs_reframed;
    const multivector_t& rhs_ref = (rhs.m_frame == our_frame)
      ? rhs
      : rhs_reframed;

    typedef typename multivector_t::matrix_t matrix_t;
#if defined(_GLUCAT_USE_DENSE_MATRICES)
    typedef typename matrix_t::size_type matrix_index_t;

    const matrix_index_t dim = lhs_ref.m_matrix.size1();
    multivector_t result = multivector_t(matrix_t(dim, dim), our_frame);
    result.m_matrix.clear();
    ublas::axpy_prod(lhs_ref.m_matrix, rhs_ref.m_matrix, result.m_matrix, true);
    return result;
#else
    typedef typename matrix_t::expression_type expression_t;

    return
      multivector_t(ublas::sparse_prod<expression_t>(lhs_ref.m_matrix, rhs_ref.m_matrix), our_frame);
#endif
  }

  /// Geometric product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator*= (const multivector_t& rhs)
  { return *this = *this * rhs; }

  /// Outer product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  operator^ (const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs)
  {
    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    typedef typename multivector_t::framed_multi_t framed_multi_t;
    return framed_multi_t(lhs) ^ framed_multi_t(rhs);
  }

  /// Outer product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator^= (const multivector_t& rhs)
  { return *this = *this ^ rhs; }

  /// Inner product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  operator& (const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs)
  {
    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    typedef typename multivector_t::framed_multi_t framed_multi_t;
    return framed_multi_t(lhs) & framed_multi_t(rhs);
  }

  /// Inner product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator&= (const multivector_t& rhs)
  { return *this = *this & rhs; }

  /// Left contraction
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  operator% (const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs)
  {
    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    typedef typename multivector_t::framed_multi_t framed_multi_t;
    return framed_multi_t(lhs) % framed_multi_t(rhs);
  }

  /// Left contraction
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator%= (const multivector_t& rhs)
  { return *this = *this % rhs; }

  /// Hestenes scalar product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  star(const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs)
  { return (lhs * rhs).scalar(); }

  /// Quotient of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator/= (const Scalar_T& scr)
  { return *this *= Scalar_T(1)/scr; }

  /// Geometric quotient
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  operator/ (const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs)
  {
    typedef numeric_traits<Scalar_T> traits_t;

#if defined(_GLUCAT_CHECK_ISNAN)
    if (lhs.isnan() || rhs.isnan())
      return traits_t::NaN();
#endif

    if (rhs == Scalar_T(0))
      return traits_t::NaN();

    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    typedef typename multivector_t::index_set_t index_set_t;

    // Operate only within a common frame
    multivector_t lhs_reframed;
    multivector_t rhs_reframed;
    const index_set_t our_frame = reframe(lhs, rhs, lhs_reframed, rhs_reframed);
    const multivector_t& lhs_ref = (lhs.m_frame == our_frame)
      ? lhs
      : lhs_reframed;
    const multivector_t& rhs_ref = (rhs.m_frame == our_frame)
      ? rhs
      : rhs_reframed;

    // Solve result == lhs_ref/rhs_ref <=> result*rhs_ref == lhs_ref
    // We now solve X == B/A
    // (where X == result, B == lhs_ref.m_matrix and A == rhs_ref.m_matrix)
    // X == B/A <=> X*A == B <=> AT*XT == BT
    // So, we solve AT*XT == BT

    typedef typename multivector_t::matrix_t matrix_t;
    typedef typename matrix_t::size_type matrix_index_t;

    const matrix_t& AT = ublas::trans(rhs_ref.m_matrix);
    matrix_t LU = AT;

    typedef ublas::permutation_matrix<matrix_index_t> permutation_t;

    permutation_t pvector(AT.size1());
    if (! ublas::lu_factorize(LU, pvector))
    {
      const matrix_t& BT = ublas::trans(lhs_ref.m_matrix);
      matrix_t XT = BT;
      ublas::lu_substitute(LU, pvector, XT);
#if defined(_GLUCAT_CHECK_ISNAN)
      if (matrix::isnan(XT))
        return traits_t::NaN();
#endif

      // Iterative refinement.
      // Reference: Nicholas J. Higham, "Accuracy and Stability of Numerical Algorithms",
      // SIAM, 1996, ISBN 0-89871-355-2, Chapter 11
      if (Tune_P::div_max_steps > 0)
      {
        // matrix_t R = ublas::prod(AT, XT) - BT;
        matrix_t R = -BT;
        ublas::axpy_prod(AT, XT, R, false);
#if defined(_GLUCAT_CHECK_ISNAN)
        if (matrix::isnan(R))
          return traits_t::NaN();
#endif

        Scalar_T nr = ublas::norm_inf(R);
        if ( nr != Scalar_T(0) && !traits_t::isNaN_or_isInf(nr) )
        {
          matrix_t XTnew = XT;
          Scalar_T nrold = nr + Scalar_T(1);
          for (int
              step = 0;
              step != Tune_P::div_max_steps &&
              nr < nrold &&
              nr != Scalar_T(0) &&
              nr == nr;
              ++step)
          {
            nrold = nr;
            if (step != 0)
              XT = XTnew;
            matrix_t& D = R;
            ublas::lu_substitute(LU, pvector, D);
            XTnew -= D;
            // noalias(R) = ublas::prod(AT, XTnew) - BT;
            R = -BT;
            ublas::axpy_prod(AT, XTnew, R, false);
            nr = ublas::norm_inf(R);
          }
        }
      }
      return multivector_t(ublas::trans(XT), our_frame);
    }
    else
      // AT is singular. Return NaN
      return traits_t::NaN();
  }

  /// Geometric quotient
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator/= (const multivector_t& rhs)
  { return *this = *this / rhs; }

  /// Transformation via twisted adjoint action
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  operator| (const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs)
  { return rhs * lhs / rhs.involute(); }

  /// Transformation via twisted adjoint action
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator|= (const multivector_t& rhs)
  { return *this = rhs * *this / rhs.involute(); }

  /// Clifford multiplicative inverse
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  inv() const
  { return multivector_t(Scalar_T(1), this->m_frame) / *this; }

  /// Integer power of multivector: *this to the m
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  pow(int m) const
  { return glucat::pow(*this, m); }

  /// Outer product power of multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  outer_pow(int m) const
  {
    if (m < 0)
      throw error_t("outer_pow(m): negative exponent");
    framed_multi_t result = Scalar_T(1);
    framed_multi_t a = *this;
    for (;
        m != 0;
        m >>= 1)
    {
      if (m & 1)
        result ^= a;
      a ^= a;
    }
    return result;
  }

  /// Grade of multivector: maximum of the grades of each term
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  index_t
  matrix_multi<Scalar_T,LO,HI>::
  grade() const
  { return framed_multi_t(*this).grade(); }

  /// Frame of multivector: union of index sets of terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const index_set<LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  frame() const
  { return this->m_frame; }

  /// Subscripting: map from index set to scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  matrix_multi<Scalar_T,LO,HI>::
  operator[] (const index_set_t ist) const
  {
    // Use matrix inner product only if ist is in frame
    if ( (ist | this->m_frame) == this->m_frame)
      return matrix::inner<Scalar_T>(this->basis_element(ist), this->m_matrix);
    else
      return Scalar_T(0);
  }

  /// Grading: part where each term is a grade-vector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  operator() (index_t grade) const
  {
    if ((grade < 0) || (grade > HI-LO))
      return 0;
    else
      return (framed_multi_t(*this))(grade);
  }

  /// Scalar part
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  matrix_multi<Scalar_T,LO,HI>::
  scalar() const
  {
    const matrix_index_t dim = this->m_matrix.size1();
    return matrix::trace(this->m_matrix) / Scalar_T( double(dim) );
  }

  /// Pure part
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  pure() const
  { return *this - this->scalar(); }

  /// Even part, sum of the even grade terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  even() const
  { return framed_multi_t(*this).even(); }

  /// Odd part, sum of the odd grade terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  odd() const
  { return framed_multi_t(*this).odd(); }

  /// Vector part of multivector, as a vector_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const typename matrix_multi<Scalar_T,LO,HI>::vector_t
  matrix_multi<Scalar_T,LO,HI>::
  vector_part() const
  { return this->vector_part(this->frame(), true); }

  /// Vector part of multivector, as a vector_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const typename matrix_multi<Scalar_T,LO,HI>::vector_t
  matrix_multi<Scalar_T,LO,HI>::
  vector_part(const index_set_t frm, const bool prechecked) const
  {
    if (!prechecked && (this->frame() | frm) != frm)
      throw error_t("vector_part(frm): value is outside of requested frame");
    vector_t result;
    // If we need to enlarge the frame we may as well use a framed_multi_t
    if (this->frame() != frm)
      return framed_multi_t(*this).vector_part(frm, true);

    const index_t begin_index = frm.min();
    const index_t end_index = frm.max()+1;
    for (index_t
        idx = begin_index;
        idx != end_index;
        ++idx)
      if (frm[idx])
        // Frame may contain indices which do not correspond to a grade 1 term but
        // frame cannot omit any index corresponding to a grade 1 term
        result.push_back(
          matrix::inner<Scalar_T>(this->basis_element(index_set_t(idx)),
          this->m_matrix));
    return result;
  }

  /// Main involution, each {i} is replaced by -{i} in each term
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  involute() const
  { return framed_multi_t(*this).involute(); }

  /// Reversion, order of {i} is reversed in each term
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  reverse() const
  { return framed_multi_t(*this).reverse(); }

  /// Conjugation, conj == reverse o involute == involute o reverse
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  conj() const
  { return framed_multi_t(*this).conj(); }

  /// Quadratic form := scalar part of rev(x)*x
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  matrix_multi<Scalar_T,LO,HI>::
  quad() const
  { // scalar(conj(x)*x) = 2*quad(even(x)) - quad(x)
    // Arvind Raja ref: "old clical: quadfunction(p:pter):pterm in file compmod.pas"
    return framed_multi_t(*this).quad();
  }

  /// Scalar_T norm squared= sum of norm squared of coordinates
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  matrix_multi<Scalar_T,LO,HI>::
  norm() const
  {
    const matrix_index_t dim = this->m_matrix.size1();
    return matrix::norm_frob2(this->m_matrix) / Scalar_T( double(dim) );
  }

  /// Maximum of absolute values of components of multivector: multivector infinity norm
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  matrix_multi<Scalar_T,LO,HI>::
  max_abs() const
  { return framed_multi_t(*this).max_abs(); }

  /// Random multivector within a frame
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  random(const index_set<LO,HI> frm, Scalar_T fill)
  {
    return framed_multi<Scalar_T,LO,HI>::random(frm, fill);
  }

  /// Write multivector to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  void
  matrix_multi<Scalar_T,LO,HI>::
  write(const std::string& msg) const
  { framed_multi_t(*this).write(msg); }

  /// Write out multivector to file
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  void
  matrix_multi<Scalar_T,LO,HI>::
  write(std::ofstream& ofile, const std::string& msg) const
  {
    if (!ofile)
      throw error_t("write(ofile,msg): cannot write to output file");
    framed_multi_t(*this).write(ofile, msg);
  }

  /// Write multivector to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  std::ostream&
  operator<< (std::ostream& os, const matrix_multi<Scalar_T,LO,HI>& val)
  {
    os << typename matrix_multi<Scalar_T,LO,HI>::framed_multi_t(val);
    return os;
  }

  /// Read multivector from input
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  std::istream&
  operator>> (std::istream& s, matrix_multi<Scalar_T,LO,HI>& val)
  { // Input looks like 1.0-2.0{1,2}+3.2{3,4}
    framed_multi<Scalar_T,LO,HI> local;
    s >> local;
    // If s.bad() then we have a corrupt input
    // otherwise we are fine and can copy the resulting matrix_multi
    if (!s.bad())
      val = local;
    return s;
  }

  /// Check if a multivector contains any IEEE NaN values
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  bool
  matrix_multi<Scalar_T,LO,HI>::
  isnan() const
  {
    if (std::numeric_limits<Scalar_T>::has_quiet_NaN)
      return matrix::isnan(this->m_matrix);
    else
      return false;
  }

  /// Remove all terms with relative size smaller than limit
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  truncated(const Scalar_T& limit) const
  { return framed_multi_t(*this).truncated(limit); }

  /// Add a term, if non-zero
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator+= (const term_t& term)
  {
    if (term.second != Scalar_T(0))
      this->m_matrix.plus_assign(matrix_t(this->basis_element(term.first)) * term.second);
    return *this;
  }

  /// Inverse generalized Fast Fourier Transform
  template< typename Multivector_T, typename Matrix_T, typename Basis_Matrix_T >
  static
  Multivector_T
  fast(const Matrix_T& X, index_t level)
  {
    typedef Multivector_T framed_multi_t;

    typedef typename framed_multi_t::index_set_t index_set_t;
    typedef typename framed_multi_t::scalar_t Scalar_T;
    typedef Matrix_T matrix_t;
    typedef Basis_Matrix_T basis_matrix_t;
    typedef typename basis_matrix_t::value_type  basis_scalar_t;
    typedef numeric_traits<Scalar_T> traits_t;

    if (level == 0)
      return framed_multi_t(traits_t::to_scalar_t(X(0,0)));

    if (ublas::norm_inf(X) == 0)
      return Scalar_T(0);

    const basis_matrix_t&  I = matrix::unit<basis_matrix_t>(2);
    basis_matrix_t J(2,2,2);
    J.clear();
    J(0,1)  = basis_scalar_t(-1);
    J(1,0)  = basis_scalar_t( 1);
    basis_matrix_t K = J;
    K(0,1)  = basis_scalar_t( 1);
    basis_matrix_t JK = I;
    JK(0,0) = basis_scalar_t(-1);

    using matrix::signed_perm_nork;
    const index_set_t ist_mn   = index_set_t(-level);
    const index_set_t ist_pn   = index_set_t(level);
    const index_set_t ist_mnpn = ist_mn | ist_pn;
    if (level == 1)
    {
      typedef typename framed_multi_t::term_t term_t;
      const Scalar_T i_x  = traits_t::to_scalar_t(signed_perm_nork(I, X)(0, 0));
      const Scalar_T j_x  = traits_t::to_scalar_t(signed_perm_nork(J, X)(0, 0));
      const Scalar_T k_x  = traits_t::to_scalar_t(signed_perm_nork(K, X)(0, 0));
      const Scalar_T jk_x = traits_t::to_scalar_t(signed_perm_nork(JK,X)(0, 0));
      framed_multi_t
             result  = i_x;
             result += term_t(ist_mn,   j_x);  // j_x *  mn;
             result += term_t(ist_pn,   k_x);  // k_x *  pn;
      return result += term_t(ist_mnpn, jk_x); // jk_x * mnpn;
    }
    else
    {
      const framed_multi_t& mn   = framed_multi_t(ist_mn);
      const framed_multi_t& pn   = framed_multi_t(ist_pn);
      const framed_multi_t& mnpn = framed_multi_t(ist_mnpn);
      const framed_multi_t& i_x  = fast<framed_multi_t, matrix_t, basis_matrix_t>
                                       (signed_perm_nork(I, X), level-1);
      const framed_multi_t& j_x  = fast<framed_multi_t, matrix_t, basis_matrix_t>
                                       (signed_perm_nork(J, X), level-1);
      const framed_multi_t& k_x  = fast<framed_multi_t, matrix_t, basis_matrix_t>
                                       (signed_perm_nork(K, X), level-1);
      const framed_multi_t& jk_x = fast<framed_multi_t, matrix_t, basis_matrix_t>
                                       (signed_perm_nork(JK,X), level-1);
      framed_multi_t
             result  =  i_x.even() - jk_x.odd();
             result += (j_x.even() - k_x.odd()) * mn;
             result += (k_x.even() - j_x.odd()) * pn;
      return result += (jk_x.even() - i_x.odd()) * mnpn;
    }
  }

  /// Use generalized FFT to construct a matrix_multi_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  fast_matrix_multi(const index_set_t frm) const
  {
    if (this->m_frame == frm)
      return *this;
    else
      return (this->template fast_framed_multi<Scalar_T>()).template fast_matrix_multi<Scalar_T>(frm);
  }

  /// Use inverse generalized FFT to construct a framed_multi_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  template <typename Other_Scalar_T>
  const framed_multi<Other_Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  fast_framed_multi() const
  {
    // Determine the amount of off-centering needed
    index_t p = this->m_frame.count_pos();
    index_t q = this->m_frame.count_neg();

    const index_t bott = pos_mod(p-q, 8);
    p += std::max(gen::offset_to_super[bott],index_t(0));
    q -= std::min(gen::offset_to_super[bott],index_t(0));

    const index_t orig_p = p;
    const index_t orig_q = q;
    while (p-q > 4)
      { p -= 4; q += 4; }
    while (p-q < -3)
      { p += 4; q -= 4; }
    if (p-q > 1)
    {
      index_t old_p = p;
      p = q+1;
      q = old_p-1;
    }
    const index_t level = (p+q)/2;

    // Do the inverse fast transform
    typedef framed_multi<Other_Scalar_T,LO,HI> framed_multi_t;
    framed_multi_t val = fast<framed_multi_t, matrix_t, basis_matrix_t>(this->m_matrix, level);

    // Off-centre val
    switch (pos_mod(orig_p-orig_q, 8))
    {
    case 2:
    case 3:
    case 4:
      val.centre_qp1_pm1(p, q);
      break;
    default:
      break;
    }
    if (orig_p-orig_q > 4)
      while (p != orig_p)
        val.centre_pp4_qm4(p, q);
    if (orig_p-orig_q < -3)
      while (p != orig_p)
        val.centre_pm4_qp4(p, q);

    // Return unfolded val
    return val.unfold(this->m_frame);
  }

  /// Table of basis elements used as a cache by basis_element()
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Matrix_T >
  class basis_table :
  public std::map< std::pair< const index_set<LO,HI>, const index_set<LO,HI> >,
                   Matrix_T* >
  {
  public:
    /// Single instance of basis table
    static basis_table& basis() { static basis_table b; return b;}
  private:
    // Enforce singleton
    // Reference: A. Alexandrescu, "Modern C++ Design", Chapter 6
    basis_table() {}
    ~basis_table() {}
    basis_table(const basis_table&);
    basis_table& operator= (const basis_table&);

    /// Friend declaration to avoid compiler warning:
    /// "... only defines a private destructor and has no friends"
    /// Ref: Carlos O'Ryan, ACE http://doc.ece.uci.edu
    friend class friend_for_private_destructor;
  };

  /// Create a basis element matrix within the current frame
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const typename matrix_multi<Scalar_T,LO,HI>::basis_matrix_t
  matrix_multi<Scalar_T,LO,HI>::
  basis_element(const index_set_t& ist) const
  {
    typedef std::pair<const index_set_t, const index_set_t>     index_set_pair_t;
    const index_set_pair_t& unfolded_pair = index_set_pair_t(ist, this->m_frame);

    typedef basis_table<Scalar_T,LO,HI,basis_matrix_t>          basis_table_t;
    typedef typename basis_table_t::const_iterator              basis_const_iterator_t;
    basis_table_t& basis_cache = basis_table_t::basis();

    const index_t frame_count = this->m_frame.count();
    const bool use_cache = frame_count <= index_t(Tune_P::basis_max_count);

    if (use_cache)
    {
      const basis_const_iterator_t basis_it = basis_cache.find(unfolded_pair);
      if (basis_it != basis_cache.end())
        return *(basis_it->second);
    }
    const index_set_t folded_set = ist.fold(this->m_frame);
    const index_set_t folded_frame = this->m_frame.fold();
    const index_set_pair_t& folded_pair = index_set_pair_t(folded_set, folded_frame);
    typedef std::pair<const index_set_pair_t, basis_matrix_t*>  basis_pair_t;
    if (use_cache)
    {
      const basis_const_iterator_t basis_it = basis_cache.find(folded_pair);
      if (basis_it != basis_cache.end())
      {
        basis_matrix_t* result_ptr = basis_it->second;
        basis_cache.insert(basis_pair_t(unfolded_pair, result_ptr));
        return *result_ptr;
      }
    }
    const index_t folded_max = folded_frame.max();
    const index_t folded_min = folded_frame.min();
    const index_t p = std::max(folded_max,           index_t(0));
    const index_t q = std::max(index_t(-folded_min), index_t(0));
    const basis_matrix_t* e = (gen::generator_table<basis_matrix_t>::generator())(p, q);
    const matrix_index_t dim = 1 << offset_level(p, q);
    basis_matrix_t result = matrix::unit<basis_matrix_t>(dim);
    for (index_t
        k = folded_min;
        k <= folded_max;
        ++k)
      if (folded_set[k])
        result = matrix::mono_prod(result, e[k]);
    if (use_cache)
    {
      basis_matrix_t* result_ptr = new basis_matrix_t(result);
      basis_cache.insert(basis_pair_t(folded_pair, result_ptr));
      basis_cache.insert(basis_pair_t(unfolded_pair, result_ptr));
    }
    return result;
  }

  /// Pade' approximation
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  static
  const matrix_multi<Scalar_T,LO,HI>
  pade_approx(const int array_size, const Scalar_T a[], const Scalar_T b[], const matrix_multi<Scalar_T,LO,HI>& X)
  {
    // Pade' approximation
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [GL], Section 11.3, p572-576.

    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    typedef numeric_traits<Scalar_T> traits_t;

    if (X.isnan())
      return traits_t::NaN();

    // Array size is assumed to be even
    const int nbr_even_powers = array_size/2 - 1;

    // Create an array of even powers
    std::vector<multivector_t> XX(nbr_even_powers);
    XX[0] = X * X;
    XX[1] = XX[0] * XX[0];
    for (int
      k = 2;
      k != nbr_even_powers;
      ++k)
      XX[k] = XX[k-2] * XX[1];

    // Calculate numerator N and denominator D
    multivector_t N = a[1];
    for (int
        k = 0;
        k != nbr_even_powers;
        ++k)
      N += XX[k] * a[2*k + 3];
    N *= X;
    N += a[0];
    for (int
        k = 0;
        k != nbr_even_powers;
        ++k)
      N += XX[k] * a[2*k + 2];
    multivector_t D = b[1];
    for (int
        k = 0;
        k != nbr_even_powers;
        ++k)
      D += XX[k] * b[2*k + 3];
    D *= X;
    D += b[0];
    for (int
        k = 0;
        k != nbr_even_powers;
        ++k)
      D += XX[k] * b[2*k + 2];
    return N / D;
  }

  /// Single step of product form of Denman-Beavers square root iteration
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  static
  void
  db_step(matrix_multi<Scalar_T,LO,HI>& M, matrix_multi<Scalar_T,LO,HI>& Y)
  {
    // Reference: [CHKL]
    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    const multivector_t& invM = inv(M);
    M = ((M + invM)/Scalar_T(2) + Scalar_T(1)) / Scalar_T(2);
    Y *= (invM + Scalar_T(1)) / Scalar_T(2);
  }

  /// Product form of Denman-Beavers square root iteration
  template< typename Scalar_T, const index_t LO, const index_t HI >
  static
  const matrix_multi<Scalar_T,LO,HI>
  db_sqrt(const matrix_multi<Scalar_T,LO,HI>& val)
  {
    // Reference: [CHKL]
    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;

    if (val == Scalar_T(0))
      return val;

    typedef std::numeric_limits<Scalar_T> limits_t;
    static const Scalar_T tol = std::pow(limits_t::epsilon(), 2);
    static const Scalar_T tol2 = tol * tol;
    static const int sqrt_max_steps = Tune_P::sqrt_max_steps;
    multivector_t M = val;
    multivector_t Y = val;
    Scalar_T norm_M_1 = norm(M - Scalar_T(1));
    typedef numeric_traits<Scalar_T> traits_t;

    for (int step = 0;
        step != sqrt_max_steps && norm_M_1 > tol2;
        ++step)
    {
      if (Y.isnan())
        return traits_t::NaN();
      db_step(M, Y);
      norm_M_1 = norm(M - Scalar_T(1));
    }
    if (norm_M_1 > tol2)
      return traits_t::NaN();
    else
      return Y;
  }
}

namespace {
  /// Coefficients of numerator polynomials of Pade approximations produced by Pade1(sqrt(1+x),x,n,n)
  // Reference: [Z], Pade1
  template< typename Scalar_T >
  struct pade_sqrt_a
  {
    static const int array_size = 14;
    static const Scalar_T array[array_size];
  };

  /// Coefficients of denominator polynomials of Pade approximations produced by Pade1(sqrt(1+x),x,n,n)
  // Reference: [Z], Pade1
  template< typename Scalar_T >
  struct pade_sqrt_b
  {
    static const int array_size = 14;
    static const Scalar_T array[array_size];
  };

  template< typename Scalar_T >
  const Scalar_T pade_sqrt_a<Scalar_T>::array[pade_sqrt_a<Scalar_T>::array_size] =
  {
        1.0,               27.0/4.0,         81.0/4.0,       2277.0/64.0,
    10395.0/256.0,      32319.0/1024.0,    8721.0/512.0,    26163.0/4096.0,
    53703.0/32768.0,    36465.0/131072.0,  3861.0/131072.0,  7371.0/4194304.0,
      819.0/16777216.0,    27.0/67108864.0
  };
  template< typename Scalar_T >
  const Scalar_T pade_sqrt_b<Scalar_T>::array[pade_sqrt_b<Scalar_T>::array_size] =
  {
        1.0,               25.0/4.0,         69.0/4.0,       1771.0/64.0,
     7315.0/256.0,      20349.0/1024.0,    4845.0/512.0,    12597.0/4096.0,
    21879.0/32768.0,    12155.0/131072.0,  1001.0/131072.0,  1365.0/4194304.0,
       91.0/16777216.0,     1.0/67108864.0
  };

  template< >
  struct pade_sqrt_a<float>
  {
    static const int array_size = 10;
    static const float array[array_size];
  };
  template< >
  struct pade_sqrt_b<float>
  {
    static const int array_size = 10;
    static const float array[array_size];
  };
  const float pade_sqrt_a<float>::array[pade_sqrt_a<float>::array_size] =
  {
       1.0,            19.0/4.0,        19.0/2.0,      665.0/64.0,
    1729.0/256.0,    2717.0/1024.0,    627.0/1024.0,   627.0/8192.0,
     285.0/65536.0,    19.0/262144.0
  };
  const float pade_sqrt_b<float>::array[pade_sqrt_a<float>::array_size] =
  {
       1.0,            17.0/4.0,        15.0/2.0,      455.0/64.0,
    1001.0/256.0,    1287.0/1024.0,    231.0/1024.0,   165.0/8192.0,
      45.0/65536,       1.0/262144.0
  };

  template< >
  struct pade_sqrt_a<long double>
  {
    static const int array_size = 18;
    static const  long double array[array_size];
  };
  template< >
  struct pade_sqrt_b<long double>
  {
    static const int array_size = 18;
    static const long double array[array_size];
  };
  const long double pade_sqrt_a<long double>::array[pade_sqrt_a<long double>::array_size] =
  {
        1.0L,                   35.0L/4.0L,             35.0L,               5425.0L/64.0L,
    35525.0L/256.0L,        166257.0L/1024.0L,      143325.0L/1024.0L,     740025.0L/8192.0L,
  2877875.0L/65536.0L,     4206125.0L/262144.0L,    572033.0L/131072.0L,  1820105.0L/2097152.0L,
  1028755.0L/8388608.0L,    395675.0L/33554432.0L,   24225.0L/33554432.0L,   6783.0L/268435456.0L,
     1785.0L/4294967296.0L,     35.0L/17179869184.0L
  };
  const long double pade_sqrt_b<long double>::array[pade_sqrt_a<long double>::array_size] =
  {
        1.0L,                   33.0L/4.0L,             31.0L,               4495.0L/64.0L,
    27405.0L/256.0L,        118755.0L/1024.0L,       94185.0L/1024.0L,     444015.0L/8192.0L,
  1562275.0L/65536.0L,     2042975.0L/262144.0L,    245157.0L/131072.0L,   676039.0L/2097152.0L,
   323323.0L/8388608.0L,    101745.0L/33554432.0L,    4845.0L/33554432.0L,    969.0L/268435456.0L,
      153.0L/4294967296.0L,      1.0L/17179869184.0L
  };

#if defined(_GLUCAT_USE_QD)
  template< >
  struct pade_sqrt_a<dd_real>
  {
    static const int array_size = 22;
    static const  dd_real array[array_size];
  };
  template< >
  struct pade_sqrt_b<dd_real>
  {
    static const int array_size = 22;
    static const dd_real array[array_size];
  };
  const dd_real pade_sqrt_a<dd_real>::array[pade_sqrt_a<dd_real>::array_size] =
  {
          dd_real("1"),                               dd_real("43")/dd_real("4"),
        dd_real("215")/dd_real("4"),               dd_real("10621")/dd_real("64"),
      dd_real("90687")/dd_real("256"),            dd_real("567987")/dd_real("1024"),
     dd_real("168861")/dd_real("256"),           dd_real("1246355")/dd_real("2048"),
    dd_real("7228859")/dd_real("16384"),        dd_real("16583853")/dd_real("65536"),
    dd_real("7538115")/dd_real("65536"),       dd_real("173376645")/dd_real("4194304"),
  dd_real("195747825")/dd_real("16777216"),    dd_real("171655785")/dd_real("67108864"),
   dd_real("14375115")/dd_real("33554432"),     dd_real("14375115")/dd_real("268435456"),
   dd_real("20764055")/dd_real("4294967296"),    dd_real("5167525")/dd_real("17179869184"),
     dd_real("206701")/dd_real("17179869184"),     dd_real("76153")/dd_real("274877906944"),
       dd_real("3311")/dd_real("1099511627776") ,     dd_real("43")/dd_real("4398046511104")
  };
  const dd_real pade_sqrt_b<dd_real>::array[pade_sqrt_a<dd_real>::array_size] =
{
          dd_real("1"),                               dd_real("41")/dd_real("4"),
        dd_real("195")/dd_real("4"),                dd_real("9139")/dd_real("64"),
      dd_real("73815")/dd_real("256"),            dd_real("435897")/dd_real("1024"),
     dd_real("121737")/dd_real("256"),            dd_real("840565")/dd_real("2048"),
    dd_real("4539051")/dd_real("16384"),         dd_real("9641775")/dd_real("65536"),
    dd_real("4032015")/dd_real("65536"),        dd_real("84672315")/dd_real("4194304"),
   dd_real("86493225")/dd_real("16777216"),     dd_real("67863915")/dd_real("67108864"),
    dd_real("5014575")/dd_real("33554432"),      dd_real("4345965")/dd_real("268435456"),
    dd_real("5311735")/dd_real("4294967296"),    dd_real("1081575")/dd_real("17179869184"),
      dd_real("33649")/dd_real("17179869184"),      dd_real("8855")/dd_real("274877906944"),
        dd_real("231")/dd_real("1099511627776"),       dd_real("1")/dd_real("4398046511104")
  };

  template< >
  struct pade_sqrt_a<qd_real>
  {
    static const int array_size = 34;
    static const  qd_real array[array_size];
  };
  template< >
  struct pade_sqrt_b<qd_real>
  {
    static const int array_size = 34;
    static const qd_real array[array_size];
  };
  const qd_real pade_sqrt_a<qd_real>::array[pade_sqrt_a<qd_real>::array_size] =
  {
              qd_real("1"),                                            qd_real("67")/qd_real("4"),
            qd_real("134"),                                         qd_real("43617")/qd_real("64"),
         qd_real("633485")/qd_real("256"),                        qd_real("6992857")/qd_real("1024"),
       qd_real("15246721")/qd_real("1024"),                     qd_real("215632197")/qd_real("8192"),
     qd_real("2518145487")/qd_real("65536"),                  qd_real("12301285425")/qd_real("262144"),
     qd_real("6344873535")/qd_real("131072"),                 qd_real("89075432355")/qd_real("2097152"),
   qd_real("267226297065")/qd_real("8388608"),               qd_real("687479618945")/qd_real("33554432"),
   qd_real("379874182975")/qd_real("33554432"),             qd_real("1443521895305")/qd_real("268435456"),
  qd_real("9425348845815")/qd_real("4294967296"),          qd_real("13195488384141")/qd_real("17179869184"),
   qd_real("987417498133")/qd_real("4294967296"),           qd_real("8055248011085")/qd_real("137438953472"),
  qd_real("6958363175533")/qd_real("549755813888"),         qd_real("5056698705201")/qd_real("2199023255552"),
   qd_real("766166470485")/qd_real("2199023255552"),         qd_real("766166470485")/qd_real("17592186044416"),
   qd_real("623623871325")/qd_real("140737488355328"),       qd_real("203123203803")/qd_real("562949953421312"),
     qd_real("6478601247")/qd_real("281474976710656"),         qd_real("5038912081")/qd_real("4503599627370496"),
      qd_real("719844583")/qd_real("18014398509481984"),         qd_real("71853815")/qd_real("72057594037927936"),
        qd_real("1165197")/qd_real("72057594037927936"),            qd_real("87703")/qd_real("576460752303423488"),
          qd_real("12529")/qd_real("18446744073709551616"),            qd_real("67")/qd_real("73786976294838206464")
  };
  const qd_real pade_sqrt_b<qd_real>::array[pade_sqrt_a<qd_real>::array_size] =
  {
            qd_real("1"),                                              qd_real("65")/qd_real("4"),
          qd_real("126"),                                           qd_real("39711")/qd_real("64"),
       qd_real("557845")/qd_real("256"),                          qd_real("5949147")/qd_real("1024"),
     qd_real("12515965")/qd_real("1024"),                       qd_real("170574723")/qd_real("8192"),
   qd_real("1916797311")/qd_real("65536"),                     qd_real("8996462475")/qd_real("262144"),
   qd_real("4450881435")/qd_real("131072"),                   qd_real("59826782925")/qd_real("2097152"),
 qd_real("171503444385")/qd_real("8388608"),                 qd_real("420696483235")/qd_real("33554432"),
 qd_real("221120793075")/qd_real("33554432"),                qd_real("797168807855")/qd_real("268435456"),
qd_real("4923689695575")/qd_real("4294967296"),             qd_real("6499270398159")/qd_real("17179869184"),
 qd_real("456864812569")/qd_real("4294967296"),             qd_real("3486599885395")/qd_real("137438953472"),
qd_real("2804116503573")/qd_real("549755813888"),           qd_real("1886827875075")/qd_real("2199023255552"),
 qd_real("263012370465")/qd_real("2199023255552"),           qd_real("240141729555")/qd_real("17592186044416"),
 qd_real("176848560525")/qd_real("140737488355328"),          qd_real("51538723353")/qd_real("562949953421312"),
   qd_real("1450433115")/qd_real("281474976710656"),            qd_real("977699359")/qd_real("4503599627370496"),
    qd_real("118183439")/qd_real("18014398509481984"),            qd_real("9652005")/qd_real("72057594037927936"),
       qd_real("121737")/qd_real("72057594037927936"),               qd_real("6545")/qd_real("576460752303423488"),
          qd_real("561")/qd_real("18446744073709551616"),               qd_real("1")/qd_real("73786976294838206464")
  };
#endif
}

namespace glucat
{
  /// Square root of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  sqrt(const matrix_multi<Scalar_T,LO,HI>& val, const matrix_multi<Scalar_T,LO,HI>& i, bool prechecked)
  {
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [Z], Pade1

    typedef numeric_traits<Scalar_T> traits_t;

    if (val.isnan())
      return traits_t::NaN();

    check_complex(val, i, prechecked);

    switch (Tune_P::function_precision)
    {
    case precision_demoted:
      {
        typedef typename traits_t::demoted::type demoted_scalar_t;
        typedef matrix_multi<demoted_scalar_t,LO,HI> demoted_multivector_t;

        const demoted_multivector_t& demoted_val = demoted_multivector_t(val);
        const demoted_multivector_t& demoted_i = demoted_multivector_t(i);

        return matrix_sqrt(demoted_val, demoted_i);
      }
      break;
    case precision_promoted:
      {
        typedef typename traits_t::promoted::type promoted_scalar_t;
        typedef matrix_multi<promoted_scalar_t,LO,HI> promoted_multivector_t;

        const promoted_multivector_t& promoted_val = promoted_multivector_t(val);
        const promoted_multivector_t& promoted_i = promoted_multivector_t(i);

        return matrix_sqrt(promoted_val, promoted_i);
      }
      break;
    default:
      return matrix_sqrt(val, i);
    }
  }

  /// Square root of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  matrix_sqrt(const matrix_multi<Scalar_T,LO,HI>& val, const matrix_multi<Scalar_T,LO,HI>& i)
  {
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [Z], Pade1

    typedef numeric_traits<Scalar_T> traits_t;

    if (val.isnan())
      return traits_t::NaN();

    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;

    const Scalar_T realval = val.scalar();
    if (val == realval)
    {
      if (realval < Scalar_T(0))
        return i * traits_t::sqrt(-realval);
      else
        return traits_t::sqrt(realval);
    }

    static const Scalar_T sqrt_2 = traits_t::sqrt(Scalar_T(2));

#if !defined(_GLUCAT_USE_EIGENVALUES)
    const multivector_t val2 = val*val;
    const Scalar_T real_val2 = val2.scalar();
    if (val2 == real_val2 && real_val2 > Scalar_T(0))
      return matrix_sqrt(-i * val, i) * (i + Scalar_T(1)) / sqrt_2;
#endif

    // Scale val towards abs(A) == 1 or towards A == 1 as appropriate
    const Scalar_T scale =
      (realval != Scalar_T(0) && norm(val/realval - Scalar_T(1)) < Scalar_T(1))
      ? realval
      : (realval < Scalar_T(0))
        ? -abs(val)
        :  abs(val);
    const Scalar_T sqrt_scale = traits_t::sqrt(traits_t::abs(scale));
    if (traits_t::isNaN_or_isInf(sqrt_scale))
      return traits_t::NaN();

    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    multivector_t rescale = sqrt_scale;
    if (scale < Scalar_T(0))
      rescale = i * sqrt_scale;

    const multivector_t& unitval = val / scale;
    const Scalar_T max_norm = Scalar_T(1.0/4.0);

#if defined(_GLUCAT_USE_EIGENVALUES)
    multivector_t scaled_result;
    typedef typename multivector_t::matrix_t matrix_t;

    // What kind of eigenvalues does the matrix contain?
    matrix::eig_genus<matrix_t> genus = matrix::classify_eigenvalues(unitval.m_matrix);
    switch (genus.m_eig_case)
    {
    case matrix::negative_eig_case:
      scaled_result = matrix_sqrt(-i * unitval, i) * (i + Scalar_T(1)) / sqrt_2;
      break;
    case matrix::both_eig_case:
      {
        const Scalar_T safe_arg = genus.m_safe_arg;
        scaled_result = matrix_sqrt(exp(i*safe_arg) * unitval, i) * exp(-i*safe_arg/Scalar_T(2));
      }
      break;
    default:
      scaled_result =
        (norm(unitval - Scalar_T(1)) < max_norm)
          // Pade' approximation of square root
        ? pade_approx(pade_sqrt_a<Scalar_T>::array_size,
                      pade_sqrt_a<Scalar_T>::array,
                      pade_sqrt_b<Scalar_T>::array,
                      unitval - Scalar_T(1))
          // Product form of Denman-Beavers square root iteration
        : db_sqrt(unitval);
      break;
    }
    if (scaled_result.isnan())
      return traits_t::NaN();
    else
      return scaled_result * rescale;
#else
    const multivector_t& scaled_result =
      (norm(unitval - Scalar_T(1)) < max_norm)
        // Pade' approximation of square root
      ? pade_approx(pade_sqrt_a<Scalar_T>::array_size,
                    pade_sqrt_a<Scalar_T>::array,
                    pade_sqrt_b<Scalar_T>::array,
                    unitval - Scalar_T(1))
        // Product form of Denman-Beavers square root iteration
      : db_sqrt(unitval);
    if (scaled_result.isnan())
    {
      if (inv(unitval).isnan())
        return traits_t::NaN();

      const multivector_t& mi_unitval = -i * unitval;

      const multivector_t& scaled_mi_result =
        (norm(mi_unitval - Scalar_T(1)) < max_norm)
          // Pade' approximation of square root
        ? pade_approx(pade_sqrt_a<Scalar_T>::array_size,
                      pade_sqrt_a<Scalar_T>::array,
                      pade_sqrt_b<Scalar_T>::array,
                      mi_unitval - Scalar_T(1))
          // Product form of Denman-Beavers square root iteration
        : db_sqrt(mi_unitval);
      if (scaled_mi_result.isnan())
        return traits_t::NaN();
      else
        return scaled_mi_result * rescale * (i + Scalar_T(1)) / sqrt_2;
    }
    else
      return scaled_result * rescale;
#endif
  }
}

namespace {
  /// Coefficients of numerator polynomials of Pade approximations produced by Pade1(log(1+x),x,n,n)
  // Reference: [Z], Pade1
  template< typename Scalar_T >
  struct pade_log_a
  {
    static const int array_size = 14;
    static const Scalar_T array[array_size];
  };

  /// Coefficients of denominator polynomials of Pade approximations produced by Pade1(log(1+x),x,n,n)
  // Reference: [Z], Pade1
  template< typename Scalar_T >
  struct pade_log_b
  {
    static const int array_size = 14;
    static const Scalar_T array[array_size];
  };
  template< typename Scalar_T >
  const Scalar_T pade_log_a<Scalar_T>::array[pade_log_a<Scalar_T>::array_size] =
  {
         0.0,                     1.0,                    6.0,            4741.0/300.0,
      1441.0/60.0,           107091.0/4600.0,          8638.0/575.0,    263111.0/40250.0,
    153081.0/80500.0,        395243.0/1101240.0,      28549.0/688275.0, 605453.0/228813200.0,
    785633.0/10296594000.0, 1145993.0/1873980108000.0
  };
  template< typename Scalar_T >
  const Scalar_T pade_log_b<Scalar_T>::array[pade_log_b<Scalar_T>::array_size] =
  {
         1.0,                    13.0/2.0,             468.0/25.0,        1573.0/50.0,
      1573.0/46.0,            11583.0/460.0,         10296.0/805.0,       2574.0/575.0,
     11583.0/10925.0,           143.0/874.0,           572.0/37145.0,      117.0/148580.0,
        13.0/742900.0,            1.0/10400600.0
  };

  template< >
  struct pade_log_a<float>
  {
    static const int array_size = 10;
    static const float array[array_size];
  };
  template< >
  struct pade_log_b<float>
  {
    static const int array_size = 10;
    static const float array[array_size];
  };
  const float pade_log_a<float>::array[pade_log_a<float>::array_size] =
  {
      0.0,            1.0,             4.0,       1337.0/204.0,
    385.0/68.0,    1879.0/680.0,     193.0/255.0,  197.0/1820.0,
    419.0/61880.0, 7129.0/61261200.0
  };
  const float pade_log_b<float>::array[pade_log_a<float>::array_size] =
  {
      1.0,            9.0/2.0,       144.0/17.0,   147.0/17.0,
    441.0/85.0,      63.0/34.0,       84.0/221.0,    9.0/221.0,
      9.0/4862.0,     1.0/48620.0
  };

  template< >
  struct pade_log_a<long double>
  {
    static const int array_size = 18;
    static const  long double array[array_size];
  };
  template< >
  struct pade_log_b<long double>
  {
    static const int array_size = 18;
    static const long double array[array_size];
  };
  const long double pade_log_a<long double>::array[pade_log_a<long double>::array_size] =
  {
         0.0L,                       1.0L,                           8.0L,                    3835.0L/132.0L,
      8365.0L/132.0L,         11363807.0L/122760.0L,            162981.0L/1705.0L,         9036157.0L/125860.0L,
  18009875.0L/453096.0L,      44211925.0L/2718576.0L,          4149566.0L/849555.0L,      16973929.0L/16020180.0L,
    172459.0L/1068012.0L,    116317061.0L/7025382936.0L,      19679783.0L/18441630207.0L, 23763863.0L/614721006900.0L,
     50747.0L/79318839600.0L, 42142223.0L/14295951736466400.0L
  };
  const long double pade_log_b<long double>::array[pade_log_a<long double>::array_size] =
  {
         1.0L,                      17.0L/2.0L,                   1088.0L/33.0L,               850.0L/11.0L,
     41650.0L/341.0L,           140777.0L/1023.0L,             1126216.0L/9889.0L,           63206.0L/899.0L,
    790075.0L/24273.0L,          60775.0L/5394.0L,               38896.0L/13485.0L,          21658.0L/40455.0L,
     21658.0L/310155.0L,          4165.0L/682341.0L,               680.0L/2047023.0L,           34.0L/3411705.0L,
        17.0L/129644790.0L,          1.0L/2333606220
  };
#if defined(_GLUCAT_USE_QD)
  template< >
  struct pade_log_a<dd_real>
  {
    static const int array_size = 22;
    static const  dd_real array[array_size];
  };
  template< >
  struct pade_log_b<dd_real>
  {
    static const int array_size = 22;
    static const dd_real array[array_size];
  };
  const dd_real pade_log_a<dd_real>::array[pade_log_a<dd_real>::array_size] =
  {
          dd_real("0"),                                  dd_real("1"),
         dd_real("10"),                              dd_real("22781")/dd_real("492"),
      dd_real("21603")/dd_real("164"),             dd_real("5492649")/dd_real("21320"),
     dd_real("978724")/dd_real("2665"),            dd_real("4191605")/dd_real("10619"),
   dd_real("12874933")/dd_real("39442"),          dd_real("11473457")/dd_real("54612"),
    dd_real("2406734")/dd_real("22755"),         dd_real("166770367")/dd_real("4004880"),
   dd_real("30653165")/dd_real("2402928"),       dd_real("647746389")/dd_real("215195552"),
   dd_real("25346331")/dd_real("47074027"),      dd_real("278270613")/dd_real("3900419380"),
  dd_real("105689791")/dd_real("15601677520"),   dd_real("606046475")/dd_real("1379188292768"),
     dd_real("969715")/dd_real("53502994116"),    dd_real("11098301")/dd_real("26204577562592"),
     dd_real("118999")/dd_real("26204577562592"), dd_real("18858053")/dd_real("1392249205900512960")
  };
  const dd_real pade_log_b<dd_real>::array[pade_log_a<dd_real>::array_size] =
  {
          dd_real("1"),                                 dd_real("21")/dd_real("2"),
       dd_real("2100")/dd_real("41"),                dd_real("12635")/dd_real("82"),
     dd_real("341145")/dd_real("1066"),            dd_real("1037799")/dd_real("2132"),
   dd_real("11069856")/dd_real("19721"),           dd_real("9883800")/dd_real("19721"),
    dd_real("6918660")/dd_real("19721"),            dd_real("293930")/dd_real("1517"),
    dd_real("1410864")/dd_real("16687"),             dd_real("88179")/dd_real("3034"),
     dd_real("734825")/dd_real("94054"),            dd_real("305235")/dd_real("188108"),
     dd_real("348840")/dd_real("1363783"),           dd_real("40698")/dd_real("1363783"),
       dd_real("6783")/dd_real("2727566"),            dd_real("9975")/dd_real("70916716"),
        dd_real("266")/dd_real("53187537"),              dd_real("7")/dd_real("70916716"),
          dd_real("7")/dd_real("8155422340"),            dd_real("1")/dd_real("538257874440")
  };

  template< >
  struct pade_log_a<qd_real>
  {
    static const int array_size = 34;
    static const  qd_real array[array_size];
  };
  template< >
  struct pade_log_b<qd_real>
  {
    static const int array_size = 34;
    static const qd_real array[array_size];
  };
  const qd_real pade_log_a<qd_real>::array[pade_log_a<qd_real>::array_size] =
  {
                qd_real("0"),                                                          qd_real("1"),
                qd_real("16"),                                                     qd_real("95201")/qd_real("780"),
             qd_real("30721")/qd_real("52"),                                     qd_real("7416257")/qd_real("3640"),
           qd_real("1039099")/qd_real("195"),                                 qd_real("6097772319")/qd_real("555100"),
        qd_real("1564058073")/qd_real("85400"),                              qd_real("30404640205")/qd_real("1209264"),
         qd_real("725351278")/qd_real("25193"),                            qd_real("4092322670789")/qd_real("147429436"),
     qd_real("4559713849589")/qd_real("201040140"),                        qd_real("5049361751189")/qd_real("320023080"),
       qd_real("74979677195")/qd_real("8000577"),                         qd_real("16569850691873")/qd_real("3481514244"),
     qd_real("1065906022369")/qd_real("515779888"),                      qd_real("335956770855841")/qd_real("438412904800"),
  qd_real("1462444287585964")/qd_real("6041877844275"),                  qd_real("397242326339851")/qd_real("6122436215532"),
    qd_real("64211291334131")/qd_real("4373168725380"),                  qd_real("142322343550859")/qd_real("51080680851480"),
   qd_real("154355972958659")/qd_real("351179680853925"),                qd_real("167483568676259")/qd_real("2937139148960100"),
     qd_real("4230788929433")/qd_real("704913395750424"),                qd_real("197968763176019")/qd_real("392923948371995600"),
    qd_real("10537522306718")/qd_real("319250708052246425"),             qd_real("236648286272519")/qd_real("144249197475035425500"),
   qd_real("260715545088119")/qd_real("4375558990076074573500"),         qd_real("289596255666839")/qd_real("192874640282553367199880"),
     qd_real("8802625510547")/qd_real("361639950529787563499775"),       qd_real("373831661521439")/qd_real("1659204093030665341336967700"),
   qd_real("446033437968239")/qd_real("464577146048586295574350956000"),  qd_real("53676090078349")/qd_real("47386868896955802148583797512000")
      };
  const qd_real pade_log_b<qd_real>::array[pade_log_a<qd_real>::array_size] =
  {
                 qd_real("1"),                                                        qd_real("33")/qd_real("2"),
              qd_real("8448")/qd_real("65"),                                       qd_real("42284")/qd_real("65"),
            qd_real("211420")/qd_real("91"),                                      qd_real("573562")/qd_real("91"),
          qd_real("32119472")/qd_real("2379"),                                  qd_real("92917044")/qd_real("3965"),
         qd_real("603960786")/qd_real("17995"),                                qd_real("144626625")/qd_real("3599"),
        qd_real("2776831200")/qd_real("68381"),                              qd_real("16692542100")/qd_real("478667"),
       qd_real("12241197540")/qd_real("478667"),                              qd_real("1098569010")/qd_real("68381"),
       qd_real("31387686000")/qd_real("3624193"),                             qd_real("9939433900")/qd_real("2479711"),
       qd_real("67091178825")/qd_real("42155087"),                            qd_real("2683647153")/qd_real("4959422"),
       qd_real("19083713088")/qd_real("121505839"),                           qd_real("4708152900")/qd_real("121505839"),
         qd_real("941630580")/qd_real("116546417"),                             qd_real("88704330")/qd_real("62755763"),
          qd_real("12902448")/qd_real("62755763"),                               qd_real("1542684")/qd_real("62755763"),
           qd_real("6427850")/qd_real("2698497809"),                             qd_real("3471039")/qd_real("18889484663"),
           qd_real("8544096")/qd_real("774468871183"),                             qd_real("39556")/qd_real("79027435835"),
            qd_real("118668")/qd_real("7191496660985"),                            qd_real("10230")/qd_real("27327687311743"),
              qd_real("5456")/qd_real("1011124430534491"),                            qd_real("44")/qd_real("1011124430534491"),
                qd_real("11")/qd_real("70778710137414370"),                            qd_real("1")/qd_real("7219428434016265740")
  };
#endif
}

namespace glucat{
  /// Pade' approximation of log
  template< typename Scalar_T, const index_t LO, const index_t HI >
  static
  const matrix_multi<Scalar_T,LO,HI>
  pade_log(const matrix_multi<Scalar_T,LO,HI>& val)
  {
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [CHKL]
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [Z], Pade1

    typedef numeric_traits<Scalar_T> traits_t;
    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();
    else
      return pade_approx(pade_log_a<Scalar_T>::array_size,
                         pade_log_a<Scalar_T>::array,
                         pade_log_b<Scalar_T>::array,
                         val - Scalar_T(1));
  }

  /// Incomplete square root cascade and Pade' approximation of log
  template< typename Scalar_T, const index_t LO, const index_t HI >
  static
  const matrix_multi<Scalar_T,LO,HI>
  cascade_log(const matrix_multi<Scalar_T,LO,HI>& val)
  {
    // Reference: [CHKL]
    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    typedef numeric_traits<Scalar_T> traits_t;
    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();

    typedef std::numeric_limits<Scalar_T> limits_t;
    static const Scalar_T epsilon = limits_t::epsilon();
    static const Scalar_T max_inner_norm = traits_t::pow(epsilon, 2);
    static const Scalar_T max_outer_norm = Scalar_T(6.0/limits_t::digits);
    multivector_t Y = val;
    multivector_t E = Scalar_T(0);
    Scalar_T norm_Y_1;
    int outer_step;
    Scalar_T pow_2_outer_step = Scalar_T(1);
    Scalar_T pow_4_outer_step = Scalar_T(1);
    for (outer_step = 0, norm_Y_1 = norm(Y - Scalar_T(1));
        outer_step != Tune_P::log_max_outer_steps && norm_Y_1 * pow_2_outer_step > max_outer_norm;
        ++outer_step,    norm_Y_1 = norm(Y - Scalar_T(1)))
    {
      if (Y == Scalar_T(0) || Y.isnan())
        return traits_t::NaN();

      // Incomplete product form of Denman-Beavers square root iteration
      multivector_t M = Y;
      for (int
          inner_step = 0;
          inner_step != Tune_P::log_max_inner_steps &&
            norm(M - Scalar_T(1)) * pow_4_outer_step > max_inner_norm;
          ++inner_step)
        db_step(M, Y);

      E += (M - Scalar_T(1)) * pow_2_outer_step;
      pow_2_outer_step *= Scalar_T(2);
      pow_4_outer_step *= Scalar_T(4);
    }
    if (outer_step == Tune_P::log_max_outer_steps && norm_Y_1 * pow_2_outer_step > max_outer_norm)
      return traits_t::NaN();
    else
      return pade_log(Y) * pow_2_outer_step - E;
  }

  /// Natural logarithm of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  log(const matrix_multi<Scalar_T,LO,HI>& val, const matrix_multi<Scalar_T,LO,HI>& i, bool prechecked)
  {
    typedef numeric_traits<Scalar_T> traits_t;

    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();

    check_complex(val, i, prechecked);

    switch (Tune_P::function_precision)
    {
    case precision_demoted:
      {
        typedef typename traits_t::demoted::type demoted_scalar_t;
        typedef matrix_multi<demoted_scalar_t,LO,HI> demoted_multivector_t;

        const demoted_multivector_t& demoted_val = demoted_multivector_t(val);
        const demoted_multivector_t& demoted_i = demoted_multivector_t(i);

        return matrix_log(demoted_val, demoted_i);
      }
      break;
    case precision_promoted:
      {
        typedef typename traits_t::promoted::type promoted_scalar_t;
        typedef matrix_multi<promoted_scalar_t,LO,HI> promoted_multivector_t;

        const promoted_multivector_t& promoted_val = promoted_multivector_t(val);
        const promoted_multivector_t& promoted_i = promoted_multivector_t(i);

        return matrix_log(promoted_val, promoted_i);
      }
      break;
    default:
      return matrix_log(val, i);
    }
  }

  /// Natural logarithm of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  matrix_log(const matrix_multi<Scalar_T,LO,HI>& val, const matrix_multi<Scalar_T,LO,HI>& i)
  {
    // Scaled incomplete square root cascade and scaled Pade' approximation of log
    // Reference: [CHKL]

    typedef numeric_traits<Scalar_T> traits_t;
    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();

    static const Scalar_T pi = traits_t::pi();
    const Scalar_T realval = val.scalar();
    if (val == realval)
    {
      if (realval < Scalar_T(0))
        return i * pi + traits_t::log(-realval);
      else
        return traits_t::log(realval);
    }
    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
#if !defined(_GLUCAT_USE_EIGENVALUES)
    const multivector_t val2 = val*val;
    const Scalar_T real_val2 = val2.scalar();
    if (val2 == real_val2 && real_val2 > 0)
      return matrix_log(-i * val, i) + i * pi/Scalar_T(2);
#endif
    // Scale val towards abs(A) == 1 or towards A == 1 as appropriate
    const Scalar_T max_norm = Scalar_T(1.0/9.0);
    const Scalar_T scale =
      (realval != Scalar_T(0) && norm(val/realval - Scalar_T(1)) < max_norm)
      ? realval
      : (realval < Scalar_T(0))
        ? -abs(val)
        :  abs(val);
    if (scale == Scalar_T(0))
      return traits_t::NaN();

    const Scalar_T log_scale = traits_t::log(traits_t::abs(scale));
    multivector_t rescale = log_scale;
    if (scale < Scalar_T(0))
      rescale = i * pi + log_scale;
    const multivector_t unitval = val/scale;
    if (inv(unitval).isnan())
      return traits_t::NaN();
#if defined(_GLUCAT_USE_EIGENVALUES)
    multivector_t scaled_result;
    typedef typename multivector_t::matrix_t matrix_t;

    // What kind of eigenvalues does the matrix contain?
    matrix::eig_genus<matrix_t> genus = matrix::classify_eigenvalues(unitval.m_matrix);
    switch (genus.m_eig_case)
    {
    case matrix::negative_eig_case:
      scaled_result = matrix_log(-i * unitval, i) + i * pi/Scalar_T(2);
      break;
    case matrix::both_eig_case:
      {
        const Scalar_T safe_arg = genus.m_safe_arg;
        scaled_result = matrix_log(exp(i*safe_arg) * unitval, i) - i * safe_arg;
      }
      break;
    default:
      scaled_result = cascade_log(unitval);
      break;
    }
#else
    multivector_t scaled_result = cascade_log(unitval);
#endif
    if (scaled_result.isnan())
      return traits_t::NaN();
    else
      return scaled_result + rescale;
  }

  /// Exponential of multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  exp(const matrix_multi<Scalar_T,LO,HI>& val)
  {
    typedef numeric_traits<Scalar_T> traits_t;
    if (val.isnan())
      return traits_t::NaN();

    const Scalar_T s = scalar(val);
    if (val == s)
      return traits_t::exp(s);

    switch (Tune_P::function_precision)
    {
    case precision_demoted:
      {
        typedef typename traits_t::demoted::type demoted_scalar_t;
        typedef matrix_multi<demoted_scalar_t,LO,HI> demoted_multivector_t;

        const demoted_multivector_t& demoted_val = demoted_multivector_t(val);
        return clifford_exp(demoted_val);
      }
      break;
    case precision_promoted:
      {
        typedef typename traits_t::promoted::type promoted_scalar_t;
        typedef matrix_multi<promoted_scalar_t,LO,HI> promoted_multivector_t;

        const promoted_multivector_t& promoted_val = promoted_multivector_t(val);
        return clifford_exp(promoted_val);
      }
      break;
    default:
      return clifford_exp(val);
    }
  }
}
#endif  // _GLUCAT_MATRIX_MULTI_IMP_H
