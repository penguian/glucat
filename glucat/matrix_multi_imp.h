#ifndef _GLUCAT_MATRIX_MULTI_IMP_H
#define _GLUCAT_MATRIX_MULTI_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_multi_imp.h : Implement the matrix representation of a multivector
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2010 by Paul C. Leopardi
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

  /// Construct a multivector, within a given frame, from a given multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const multivector_t& val, const index_set_t frm, const bool prechecked)
  : m_frame( frm )
  {
    if (!prechecked && (val.m_frame | frm) != frm)
      throw error_t("multivector_t(val,frm): cannot initialize with value outside of frame");
    if (frm == val.m_frame)
      this->m_matrix = val.m_matrix;
    else
      *this = multivector_t(framed_multi_t(val), frm, true);
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
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const framed_multi_t& val)
  : m_frame( val.frame() )
  {
    if (val.size() >= Tune_P::fast_size_threshold)
      try
      {
        *this = val.fast_matrix_multi(this->m_frame);
        return;
      }
      catch (const glucat_error& e)
      { }
    const matrix_index_t dim = folded_dim<matrix_index_t>(this->m_frame);
    this->m_matrix.resize(dim, dim, false);
    this->m_matrix.clear();

    for (typename framed_multi_t::const_iterator
        val_it = val.begin();
        val_it != val.end();
        ++val_it)
      *this += *val_it;
  }

  /// Construct a multivector, within a given frame, from a framed_multi_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const framed_multi_t& val, const index_set_t frm, const bool prechecked)
  : m_frame( frm )
  {
    if (!prechecked && (val.frame() | frm) != frm)
      throw error_t("multivector_t(val,frm): cannot initialize with value outside of frame");

    if (val.size() >= Tune_P::fast_size_threshold)
      try
      {
        *this = val.fast_matrix_multi(frm);
        return;
      }
      catch (const glucat_error& e)
      { }
    const matrix_index_t dim = folded_dim<matrix_index_t>(frm);
    this->m_matrix.resize(dim, dim, false);
    this->m_matrix.clear();

    for (typename framed_multi_t::const_iterator
        val_it = val.begin();
        val_it != val.end();
        ++val_it)
      *this += *val_it;
  }

  /// Construct a multivector within a given frame from a given matrix
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const matrix_t& mtx, const index_set_t frm)
  : m_frame( frm )
  {
    const matrix_index_t dim = folded_dim<matrix_index_t>(frm);
    this->m_matrix.resize(dim, dim, false);
    noalias(this->m_matrix) = mtx;
  }

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
    const index_set_t our_frame = this->m_frame | rhs.m_frame;
    const multivector_t& lhs_ref =
      (this->m_frame == our_frame)
      ? *this
      : multivector_t(framed_multi_t(*this), our_frame, true);
    const multivector_t& rhs_ref =
      (rhs.m_frame == our_frame)
      ? rhs
      : multivector_t(framed_multi_t(rhs), our_frame, true);

#ifdef _GLUCAT_USE_DENSE_MATRICES
    return ublas::norm_inf(lhs_ref.m_matrix - rhs_ref.m_matrix) == 0;
#else
    typedef typename matrix_t::const_iterator1 const_iterator1;
    typedef typename matrix_t::const_iterator2 const_iterator2;
    // If either matrix contains zero entries,
    // compare using subtraction and ublas::norm_inf
    for (const_iterator1
        it1 = lhs_ref.m_matrix.begin1();
        it1 != lhs_ref.m_matrix.end1();
        ++it1)
      for (const_iterator2
          it2 = it1.begin();
          it2 != it1.end();
          ++it2)
        if (*it2 == 0)
          return ublas::norm_inf(lhs_ref.m_matrix - rhs_ref.m_matrix) == 0;
    for (const_iterator1
        it1 = rhs.m_matrix.begin1();
        it1 != rhs.m_matrix.end1();
        ++it1)
      for (const_iterator2
          it2 = it1.begin();
          it2 != it1.end();
          ++it2)
        if (*it2 == 0)
          return ublas::norm_inf(lhs_ref.m_matrix - rhs_ref.m_matrix) == 0;
    // Neither matrix contains zero entries.
    // Compare by iterating over both matrices in lock step.
    const_iterator1 this_it1 = lhs_ref.m_matrix.begin1();
    const_iterator1 it1 = rhs_ref.m_matrix.begin1();
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
    const index_set_t our_frame = this->m_frame | rhs.m_frame;
    if (this->m_frame != our_frame)
      // Represent *this in our_frame via conversion through framed_multi_t
      *this = multivector_t(framed_multi_t(*this), our_frame, true);
    const multivector_t& rhs_ref =
      (rhs.m_frame == our_frame)
      ? rhs
      : multivector_t(framed_multi_t(rhs), our_frame, true);

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
    const index_set_t our_frame = this->m_frame | rhs.m_frame;
    if (this->m_frame != our_frame)
      // Represent *this in our_frame via conversion through framed_multi_t
      *this = multivector_t(framed_multi_t(*this), our_frame, true);
    const multivector_t& rhs_ref =
      (rhs.m_frame == our_frame)
      ? rhs
      : multivector_t(framed_multi_t(rhs), our_frame, true);

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
    if (numeric_traits<Scalar_T>::isNaN_or_isInf(scr) || this->isnan())
      return *this = numeric_traits<Scalar_T>::NaN();
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
    typedef typename multivector_t::framed_multi_t framed_multi_t;

#ifndef _GLUCAT_USE_DENSE_MATRICES
    if (lhs.isnan() || rhs.isnan())
      return numeric_traits<Scalar_T>::NaN();
#endif

    // Operate only within a common frame
    const index_set_t our_frame = lhs.m_frame | rhs.m_frame;
    const multivector_t& lhs_ref = (lhs.m_frame == our_frame)
                                 ? lhs
                                 : multivector_t(framed_multi_t(lhs), our_frame, true);
    const multivector_t& rhs_ref = (rhs.m_frame == our_frame)
                                 ? rhs
                                 : multivector_t(framed_multi_t(rhs), our_frame, true);

    typedef typename multivector_t::matrix_t matrix_t;
#ifdef _GLUCAT_USE_DENSE_MATRICES
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
#ifndef _GLUCAT_USE_DENSE_MATRICES
    if (lhs.isnan() || rhs.isnan())
      return numeric_traits<Scalar_T>::NaN();
#endif

    if (rhs == Scalar_T(0))
      return numeric_traits<Scalar_T>::NaN();

    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    typedef typename multivector_t::index_set_t index_set_t;
    typedef typename multivector_t::framed_multi_t framed_multi_t;

    // Operate only within a common frame
    const index_set_t our_frame = lhs.m_frame | rhs.m_frame;
    const multivector_t& lhs_ref = (lhs.m_frame == our_frame)
                                 ? lhs
                                 : multivector_t(framed_multi_t(lhs), our_frame, true);
    const multivector_t& rhs_ref = (rhs.m_frame == our_frame)
                                 ? rhs
                                 : multivector_t(framed_multi_t(rhs), our_frame, true);

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
#ifndef _GLUCAT_USE_DENSE_MATRICES
      if (matrix::isnan(XT))
        return numeric_traits<Scalar_T>::NaN();
#endif

      // Iterative refinement.
      // Reference: Nicholas J. Higham, "Accuracy and Stability of Numerical Algorithms",
      // SIAM, 1996, ISBN 0-89871-355-2, Chapter 11
      if (Tune_P::div_max_steps > 0)
      {
        // matrix_t R = ublas::prod(AT, XT) - BT;
        matrix_t R = -BT;
        ublas::axpy_prod(AT, XT, R, false);
#ifndef _GLUCAT_USE_DENSE_MATRICES
        if (matrix::isnan(R))
          return numeric_traits<Scalar_T>::NaN();
#endif

        Scalar_T nr = ublas::norm_inf(R);
        if ( nr != Scalar_T(0) && !numeric_traits<Scalar_T>::isNaN_or_isInf(nr) )
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
      return numeric_traits<Scalar_T>::NaN();
  }

  /// Geometric quotient
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator/= (const multivector_t& rhs)
  { return *this = *this / rhs; }

  /// Clifford multiplicative inverse
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  inv() const
  { return multivector_t(Scalar_T(1), this->m_frame) / *this; }

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

  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  pow(int m) const
  { return glucat::pow(*this, m); }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  outer_pow(int m) const
  { // outer product power
    if (m < 0)
      throw error_t("outer_pow(m): negative exponent");
    multivector_t result = Scalar_T(1);
    multivector_t a = *this;
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
  {
    vector_t result;
    const index_t begin_index = this->m_frame.min();
    const index_t end_index = this->m_frame.max()+1;
    for (index_t
        idx = begin_index;
        idx != end_index;
        ++idx)
      if (this->m_frame[idx])
        // Frame may contain indices which do not correspond to a grade 1 term but
        // frame cannot omit any index corresponding to a grade 1 term
        result.push_back(
          matrix::inner<Scalar_T>(this->basis_element(index_set_t(idx)),
          this->m_matrix));
    return result;
  }

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

  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  Scalar_T
  matrix_multi<Scalar_T,LO,HI>::
  max_abs() const
  { return framed_multi_t(*this).max_abs(); }

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

  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  truncated(const Scalar_T& limit) const
  { return framed_multi_t(*this).truncated(limit); }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const index_set<LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  frame() const
  { return this->m_frame; }

  /// Add a term, if non-zero
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator+= (const term_t& term)
  {
    if (term.second != Scalar_T(0))
      this->m_matrix.plus_assign(this->basis_element(term.first) * term.second);
    return *this;
  }

  /// Inverse generalized Fast Fourier Transform
  template< typename Multivector_T, typename Matrix_T, typename Basis_Matrix_T >
  static
  Multivector_T
  fast(const Matrix_T& X, index_t level)
  {
    typedef Multivector_T framed_multi_t;

    if (level == 0)
      return framed_multi_t(X(0,0));

    typedef typename framed_multi_t::matrix_multi_t matrix_multi_t;
    typedef typename framed_multi_t::index_set_t index_set_t;
    typedef typename framed_multi_t::scalar_t Scalar_T;
    typedef Matrix_T matrix_t;
    typedef Basis_Matrix_T basis_matrix_t;

    if (ublas::norm_inf(X) == 0)
      return Scalar_T(0);

    const basis_matrix_t&  I = matrix::unit<basis_matrix_t>(2);
    basis_matrix_t J(2,2,2);
    J.clear();
    J(0,1)  = Scalar_T(-1);
    J(1,0)  = Scalar_T( 1);
    basis_matrix_t K = J;
    K(0,1)  = Scalar_T( 1);
    basis_matrix_t JK = I;
    JK(0,0) = Scalar_T(-1);

    using matrix::nork;
    const bool mono = true;
    const index_set_t ist_mn   = index_set_t(-level);
    const index_set_t ist_pn   = index_set_t(level);
    const index_set_t ist_mnpn = ist_mn | ist_pn;
    if (level == 1)
    {
      typedef typename framed_multi_t::term_t term_t;
      const Scalar_T i_x  = nork(I, X, mono)(0, 0);
      const Scalar_T j_x  = nork(J, X, mono)(0, 0);
      const Scalar_T k_x  = nork(K, X, mono)(0, 0);
      const Scalar_T jk_x = nork(JK,X, mono)(0, 0);
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
      const framed_multi_t& i_x  = fast<framed_multi_t, matrix_t, basis_matrix_t>(nork(I, X, mono), level-1);
      const framed_multi_t& j_x  = fast<framed_multi_t, matrix_t, basis_matrix_t>(nork(J, X, mono), level-1);
      const framed_multi_t& k_x  = fast<framed_multi_t, matrix_t, basis_matrix_t>(nork(K, X, mono), level-1);
      const framed_multi_t& jk_x = fast<framed_multi_t, matrix_t, basis_matrix_t>(nork(JK,X, mono), level-1);
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
      return (this->fast_framed_multi()).fast_matrix_multi(frm);
  }

  /// Use inverse generalized FFT to construct a framed_multi_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const typename matrix_multi<Scalar_T,LO,HI>::framed_multi_t
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
                   Matrix_T >
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
    const index_set_t folded_set = ist.fold(this->m_frame);
    const index_set_t folded_frame = this->m_frame.fold();
    const index_t folded_max = folded_frame.max();
    const index_t folded_min = folded_frame.min();
    const matrix_index_t dim = folded_dim<matrix_index_t>(this->m_frame);

    const index_t p = std::max(folded_max,           index_t(0));
    const index_t q = std::max(index_t(-folded_min), index_t(0));

    typedef basis_table<Scalar_T,LO,HI,basis_matrix_t>      basis_table_t;
    typedef std::pair<const index_set_t, const index_set_t> index_set_pair_t;
    typedef std::pair<const index_set_pair_t, basis_matrix_t>   basis_pair_t;
    typedef typename basis_table_t::const_iterator          basis_const_iterator_t;
    basis_table_t& basis_cache = basis_table_t::basis();
    const index_set_pair_t& folded_pair = index_set_pair_t(folded_set, folded_frame);
    if (p+q <= Tune_P::basis_max_count)
    {
      const basis_const_iterator_t basis_it = basis_cache.find(folded_pair);
      if (basis_it != basis_cache.end())
        return basis_it->second;
    }
    const basis_matrix_t* e = (gen::generator_table<basis_matrix_t>::generator())(p,q);
    basis_matrix_t result = matrix::unit<basis_matrix_t>(dim);
    for (index_t
        k = folded_min;
        k <= folded_max;
        ++k)
      if (folded_set[k])
        result = matrix::mono_prod(result, e[k]);
    if (p+q <= Tune_P::basis_max_count)
      basis_cache.insert(basis_pair_t(folded_pair, result));
    return result;
  }

  /// Pade' approximation
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  static
  const matrix_multi<Scalar_T,LO,HI>
  pade_approx(const Scalar_T a[], const Scalar_T b[], const matrix_multi<Scalar_T,LO,HI>& X)
  {
    // Pade' approximation
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [GL], Section 11.3, p572-576.

    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    typedef numeric_traits<Scalar_T> traits_t;

    if (X.isnan())
      return traits_t::NaN();

    const multivector_t& X2 = X*X;
    const multivector_t& X4 = X2*X2;
    const multivector_t& X6 = X4*X2;
    const multivector_t& N = a[0]+X2*a[2]+X4*a[4]+X6*(a[6]+X2*a[8]+X4*a[10]+X6*a[12]) +
                          X*(a[1]+X2*a[3]+X4*a[5]+X6*(a[7]+X2*a[9]+X4*a[11]+X6*a[13]));
    const multivector_t& D = b[0]+X2*b[2]+X4*b[4]+X6*(b[6]+X2*b[8]+X4*b[10]+X6*b[12]) +
                          X*(b[1]+X2*b[3]+X4*b[5]+X6*(b[7]+X2*b[9]+X4*b[11]+X6*b[13]));
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
    const multivector_t& iM = inv(M);
    M = ((M + iM)/Scalar_T(2) + Scalar_T(1)) / Scalar_T(2);
    Y *= (iM + Scalar_T(1)) / Scalar_T(2);
  }

  /// Product form of Denman-Beavers square root iteration
  template< typename Scalar_T, const index_t LO, const index_t HI >
  static
  const matrix_multi<Scalar_T,LO,HI>
  db_sqrt(const matrix_multi<Scalar_T,LO,HI>& val)
  {
    // Reference: [CHKL]
    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    typedef numeric_traits<Scalar_T> traits_t;

    if (val == Scalar_T(0))
      return val;

    multivector_t M = val;
    multivector_t Y = val;
    static const Scalar_T tol = std::numeric_limits<Scalar_T>::epsilon() / Scalar_T(1024);
    static const Scalar_T tol2 = tol * tol;
    static const int sqrt_max_steps = Tune_P::sqrt_max_steps;
    Scalar_T norm_M_1 = norm(M - Scalar_T(1));

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

  /// Square root of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  static
  const matrix_multi<Scalar_T,LO,HI>
  sqrt(const matrix_multi<Scalar_T,LO,HI>& val, const matrix_multi<Scalar_T,LO,HI>& i, bool prechecked)
  {
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [Z], Pade1
    //std::cout << "Sqrt called" << std::endl;

    // Pade approximation produced by Pade1(sqrt(1+x),x,13,13):
    // numer := 1+27/4*x+81/4*x^2+2277/64*x^3
    //          +10395/256*x^4+32319/1024*x^5+8721/512*x^6+26163/4096*x^7
    //          +53703/32768*x^8+36465/131072*x^9+3861/131072*x^10
    //          +7371/4194304*x^11+819/16777216*x^12+27/67108864*x^13;
    // denom := 1+25/4*x+69/4*x^2+1771/64*x^3
    //          +7315/256*x^4+20349/1024*x^5+4845/512*x^6+12597/4096*x^7
    //          +21879/32768*x^8+12155/131072*x^9+1001/131072*x^10
    //          +1365/4194304*x^11+91/16777216*x^12+1/67108864*x^13;

    static const Scalar_T a[] =
    {
          1.0,             27.0/4.0,         81.0/4.0,      2277.0/64.0,
      10395.0/256.0,    32319.0/1024.0,    8721.0/512.0,   26163.0/4096.0,
      53703.0/32768.0,  36465.0/131072.0,  3861.0/131072.0,
       7371.0/4194304.0,  819.0/16777216.0,  27.0/67108864.0
    };
    static const Scalar_T b[] =
    {
          1.0,             25.0/4.0,         69.0/4.0,      1771.0/64.0,
       7315.0/256.0,    20349.0/1024.0,    4845.0/512.0,   12597.0/4096.0,
      21879.0/32768.0,  12155.0/131072.0,  1001.0/131072.0,
       1365.0/4194304.0,   91.0/16777216.0,   1.0/67108864.0
    };

    typedef numeric_traits<Scalar_T> traits_t;
    if (val.isnan())
      return traits_t::NaN();

    check_complex(val, i, prechecked);

    const Scalar_T realval = real(val);
    if (val == realval)
    {
      if (realval < Scalar_T(0))
        return i * traits_t::sqrt(-realval);
      else
        return traits_t::sqrt(realval);
    }

    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    static const Scalar_T sqrt_2 = traits_t::sqrt(Scalar_T(2));

#if !defined(_GLUCAT_USE_EIGENVALUES)
    const multivector_t val2 = val*val;
    const Scalar_T real_val2 = real(val2);
    if (val2 == real_val2 && real_val2 > Scalar_T(0))
      return sqrt(-i * val, i, prechecked) * (i + Scalar_T(1)) / sqrt_2;
#endif

    // Scale val towards abs(A) == 1 or towards A == 1 as appropriate
    const Scalar_T max_norm = Scalar_T(1.0/9.0);
    const Scalar_T scale =
      (realval != Scalar_T(0) && norm(val/realval - Scalar_T(1)) < Scalar_T(1))
        ? realval
        : (realval < Scalar_T(0))
          ? -abs(val)
          : abs(val);
    const Scalar_T sqrt_scale = traits_t::sqrt(traits_t::abs(scale));
    if (traits_t::isNaN_or_isInf(sqrt_scale))
      return traits_t::NaN();

    multivector_t rescale = sqrt_scale;
    if (scale < Scalar_T(0))
      rescale = i * sqrt_scale;

    const multivector_t& unitval = val / scale;

#if defined(_GLUCAT_USE_EIGENVALUES)
    multivector_t scaled_result;
    typedef typename multivector_t::matrix_t matrix_t;

    // What kind of eigenvalues does the matrix contain?
    matrix::eig_genus<matrix_t> genus = matrix::classify_eigenvalues(unitval.m_matrix);
    switch (genus.m_eig_case)
    {
    case matrix::negative_eig_case:
      scaled_result = sqrt(-i * unitval, i, prechecked) * (i + Scalar_T(1)) / sqrt_2;
      break;
    case matrix::both_eig_case:
      {
        const Scalar_T safe_arg = genus.m_safe_arg;
        scaled_result = sqrt(exp(i*safe_arg) * unitval, i, prechecked) * exp(-i*safe_arg/Scalar_T(2));
      }
      break;
    default:
      scaled_result =
          (norm(unitval - Scalar_T(1)) < max_norm)
            // Pade' approximation of square root
            ? pade_approx(a, b, unitval - Scalar_T(1))
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
            ? pade_approx(a, b, unitval - Scalar_T(1))
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
          ? pade_approx(a, b, mi_unitval - Scalar_T(1))
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

  /// Exponential of multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  exp(const matrix_multi<Scalar_T,LO,HI>& val)
  {
    // Scaling and squaring Pade' approximation of matrix exponential
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [H]

    typedef numeric_traits<Scalar_T> traits_t;

    if (val.isnan())
      return traits_t::NaN();

    const Scalar_T scalar_val = scalar(val);
    const Scalar_T scalar_exp = traits_t::exp(scalar_val);
    if (traits_t::isNaN_or_isInf(scalar_exp))
      return traits_t::NaN();
    if (val == scalar_val)
      return scalar_exp;

    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
    multivector_t A = val - scalar_val;
    const Scalar_T pure_scale = A.norm();

    if (traits_t::isNaN_or_isInf(pure_scale))
      return traits_t::NaN();
    if (pure_scale == Scalar_T(0))
      return scalar_exp;

    const int ilog2_scale =
      std::max(0, traits_t::to_int(ceil((log2(pure_scale) + Scalar_T(A.frame().count()))/Scalar_T(2)) - 3));
    const Scalar_T i_scale = traits_t::pow(Scalar_T(2), ilog2_scale);
    if (traits_t::isNaN_or_isInf(i_scale))
      return traits_t::NaN();

    A /= i_scale;
    multivector_t pure_exp;
    {
      const int q = 13;
      static Scalar_T c[q+1];
      if (c[0] != Scalar_T(1))
      {
        c[0] = Scalar_T(1);
        for (int
            k = 0;
            k != q;
            ++k)
          c[k+1] = c[k]*(q-k) / ((2*q-k)*(k+1));
      }
      const multivector_t& A2 = A*A;
      const multivector_t& A4 = A2*A2;
      const multivector_t& A6 = A4*A2;
      const multivector_t& U =     c[0]+A2*c[2]+A4*c[4]+A6*(c[6]+A2*c[8]+A4*c[10]+A6*c[12]);
      const multivector_t& AV = A*(c[1]+A2*c[3]+A4*c[5]+A6*(c[7]+A2*c[9]+A4*c[11]+A6*c[13]));
      pure_exp = (U+AV) / (U-AV);
    }
    for (int
        k = 0;
        k != ilog2_scale;
        ++k)
      pure_exp *= pure_exp;
    return pure_exp * scalar_exp;
  }

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

    // Pade approximation produced by Pade1(log(1-x),x,13,13):
    // numer := -x+6*x^2-4741/300*x^3
    //          +1441/60*x^4-107091/4600*x^5+8638/575*x^6-263111/40250*x^7
    //          +153081/80500*x^8-395243/1101240*x^9+28549/688275*x^10
    //          -605453/228813200*x^11+785633/10296594000*x^12
    //          -1145993/1873980108000*x^13;
    // denom := 1-13/2*x+468/25*x^2-1573/50*x^3
    //          +1573/46*x^4-11583/460*x^5+10296/805*x^6-2574/575*x^7
    //          +11583/10925*x^8-143/874*x^9+572/37145*x^10
    //          -117/148580*x^11+13/742900*x^12-1/10400600*x^13;

    static const Scalar_T a[] =
    {
            0.0,                  -1.0,                      6.0,            -4741.0/300.0,
         1441.0/60.0,        -107091.0/4600.0,            8638.0/575.0,    -263111.0/40250.0,
       153081.0/80500.0,     -395243.0/1101240.0,        28549.0/688275.0,
      -605453.0/228813200.0,  785633.0/10296594000.0, -1145993.0/1873980108000.0
    };
    static const Scalar_T b[] =
    {
         1.0,                    -13.0/2.0,                468.0/25.0,       -1573.0/50.0,
      1573.0/46.0,            -11583.0/460.0,            10296.0/805.0,      -2574.0/575.0,
     11583.0/10925.0,           -143.0/874.0,              572.0/37145.0,
      -117.0/148580.0,            13.0/742900.0,            -1.0/10400600.0
    };
    typedef numeric_traits<Scalar_T> traits_t;

    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();
    else
      return pade_approx(a, b, -val + Scalar_T(1));
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

    multivector_t Y = val;
    multivector_t E = Scalar_T(0);
    static const Scalar_T epsilon = std::numeric_limits<Scalar_T>::epsilon();
    static const Scalar_T max_inner_norm = traits_t::pow(epsilon, 4);
    static const Scalar_T max_outer_norm = Scalar_T(1.0/9.0);
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
    // Scaled incomplete square root cascade and scaled Pade' approximation of log
    // Reference: [CHKL]
    //std::cout << "Log called" << std::endl;
    typedef numeric_traits<Scalar_T> traits_t;

    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();

    static const Scalar_T pi = traits_t::pi();
    const Scalar_T realval = real(val);
    if (val == realval)
    {
      if (realval < Scalar_T(0))
      {
        check_complex(val, i, prechecked);
        return i * pi + traits_t::log(-realval);
      }
      else
        return traits_t::log(realval);
    }
    typedef matrix_multi<Scalar_T,LO,HI> multivector_t;
#if !defined(_GLUCAT_USE_EIGENVALUES)
    const multivector_t val2 = val*val;
    const Scalar_T real_val2 = real(val2);
    if (val2 == real_val2 && real_val2 > 0)
    {
      check_complex(val, i, prechecked);
      return log(-i * val, i, prechecked) + i * pi/Scalar_T(2);
    }
#endif
    // Scale val towards abs(A) == 1 or towards A == 1 as appropriate
    const Scalar_T max_norm = Scalar_T(1.0/9.0);
    const Scalar_T scale =
      (realval != Scalar_T(0) && norm(val/realval - Scalar_T(1)) < max_norm)
        ? realval
        : (realval < Scalar_T(0))
          ? -abs(val)
          : abs(val);
    if (scale == Scalar_T(0))
      return traits_t::NaN();

    const Scalar_T log_scale = traits_t::log(traits_t::abs(scale));
    multivector_t rescale = log_scale;
    if (scale < Scalar_T(0))
    {
      check_complex(val, i, prechecked);
      rescale = i * pi + log_scale;
    }
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
      scaled_result = log(-i * unitval, i, prechecked) + i * pi/Scalar_T(2);
      break;
    case matrix::both_eig_case:
      {
        const Scalar_T safe_arg = genus.m_safe_arg;
        scaled_result = log(exp(i*safe_arg) * unitval, i, prechecked) - i * safe_arg;
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
}
#endif  // _GLUCAT_MATRIX_MULTI_IMP_H
