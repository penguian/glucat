#ifndef _GLUCAT_MATRIX_MULTI_IMP_H
#define _GLUCAT_MATRIX_MULTI_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_multi_imp.h : Implement the matrix representation of a multivector
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

namespace glucat
{
  /// Class name used in messages
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const char*
  matrix_multi<Scalar_T,LO,HI>::
  classname()
  { return "matrix_multi"; }

  /// Default constructor
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi()
  : m_frame( index_set_t() ),
    m_matrix( matrix_t( 1, 1 ) )
  { }

  /// Copy constructor
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const multivector_t& val)
  : m_frame( val.m_frame ),
    m_matrix( matrix_t( val.m_matrix.nrows(), val.m_matrix.ncols()) )
  {
    // This copy constructor is needed because the MTL copy constructor does a shallow copy
    // on a class containing pointers to dynamically allocated memory.
    // Here, instead, we do a deep copy to avoid all possibility of aliasing.
    // Reference: Scott Meyers, "Effective C++", 2nd Edition, Addison-Wesley, 1998, Item 11.
    mtl::copy(val.m_matrix, m_matrix);
  }

	/// Construct a multivector, within a given frame, from a given multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const multivector_t& val, const index_set_t& frm, const bool prechecked)
  : m_frame( frm )
  {
    if (!prechecked && (val.frame() | frm) != frm)
      throw error_t("multivector_t(val,frm): cannot initialize with value outside of frame");
		*this = multivector_t(framed_multi_t(val), frm, true);
	}

  /// Construct a multivector from an index set and a scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const index_set_t& ist, const Scalar_T& crd)
  : m_frame( ist )
  {
    const matrix_index_t dim = folded_dim<Scalar_T,LO,HI>(m_frame);
    m_matrix = matrix_t( dim, dim );
    *this += pair_t(ist, crd);
  }

	/// Construct a multivector, within a given frame, from an index set and a scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const index_set_t& ist, const Scalar_T& crd, const index_set_t& frm, const bool prechecked)
  : m_frame( frm )
  {
    if (!prechecked && (ist | frm) != frm)
      throw error_t("multivector_t(ist,crd,frm): cannot initialize with value outside of frame");
    const matrix_index_t dim = folded_dim<Scalar_T,LO,HI>(m_frame);
    m_matrix = matrix_t( dim, dim );
    *this += pair_t(ist, crd);
	}

  /// Construct a multivector from a scalar (within a frame, if given)
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const Scalar_T& scr, const index_set_t& frm)
  : m_frame( frm )
  {
    const matrix_index_t dim = folded_dim<Scalar_T,LO,HI>(m_frame);
    m_matrix = matrix_t( dim, dim );
    *this += pair_t(index_set_t(), scr);
  }

  /// Construct a multivector from an int (within a frame, if given)
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const int scr, const index_set_t& frm)
  { *this = multivector_t(Scalar_T(scr), frm); }

	/// Construct a multivector, within a given frame, from a given vector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const vector_t& vec,
 							 const index_set_t& frm, const bool prechecked = false)
  : m_frame( frm )
  {
    if (!prechecked && index_t(vec.size()) != frm.count())
      throw error_t("multivector_t(vec,frm): cannot initialize with vector not matching frame");
    const matrix_index_t dim = folded_dim<Scalar_T,LO,HI>(m_frame);
    m_matrix = matrix_t( dim, dim );
    vector_t::const_iterator scvec = vec.begin();
    const index_t begin_index = m_frame.min();
    const index_t end_index = m_frame.max()+1;
    for (index_t idx = begin_index; idx != end_index; ++idx)
      if (m_frame[idx])
      {
        *this += pair_t(index_set_t(idx), *scvec);
        ++scvec;
      }
  }

  /// Construct a multivector from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const string& str)
  { *this = framed_multi_t(str); }

  /// Construct a multivector, within a given frame, from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const string& str, const index_set_t& frm, const bool prechecked)
  { *this = multivector_t(framed_multi_t(str), frm, prechecked); }

	/// Construct a multivector from a framed_multi_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const framed_multi_t& val)
  : m_frame( val.frame() )
  {
    const matrix_index_t dim = folded_dim<Scalar_T,LO,HI>(m_frame);
    m_matrix = matrix_t( dim, dim );

    for(framed_multi_t::const_iterator
      scan = val.begin(); scan != val.end(); ++scan)
      *this += *scan;
  }

	/// Construct a multivector, within a given frame, from a framed_multi_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const framed_multi_t& val, const index_set_t& frm, const bool prechecked)
  : m_frame( frm )
  {
    if (!prechecked && (val.frame() | frm) != frm)
      throw error_t("multivector_t(val,frm): cannot initialize with value outside of frame");
    const matrix_index_t dim = folded_dim<Scalar_T,LO,HI>(m_frame);
    m_matrix = matrix_t( dim, dim );

    for(framed_multi_t::const_iterator
      scan = val.begin(); scan != val.end(); ++scan)
      *this += *scan;
  }

  /// Construct a multivector within a given frame from a given matrix
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>::
  matrix_multi(const matrix_t& mtx, const index_set_t& frm)
  : m_frame( frm )
  {
    const matrix_index_t dim = folded_dim<Scalar_T,LO,HI>(m_frame);
    m_matrix = matrix_t( dim, dim );
    mtl::copy(mtx, m_matrix);
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
    m_frame = rhs.m_frame;
    m_matrix = matrix_t( rhs.m_matrix.nrows(), rhs.m_matrix.ncols() );
    // This assignment operator is needed because the MTL assignment operator does a shallow copy
    // on a class containing pointers to dynamically allocated memory.
    // Here, instead, we do a deep copy to avoid all possibility of aliasing.
    // Reference: Scott Meyers, "Effective C++", 2nd Edition, Addison-Wesley, 1998, Item 11.
    mtl::copy(rhs.m_matrix, m_matrix);
    return *this;
  }

  /// Test for equality of multivectors
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const bool
  matrix_multi<Scalar_T,LO,HI>::
  operator==  (const multivector_t& rhs) const
  {
    // Compare only within a common frame
    if (m_frame != rhs.m_frame)
      return framed_multi_t(*this) == framed_multi_t(rhs);

    typedef typename matrix_t::const_iterator matrix_iterator;
    matrix_iterator ithis =          m_matrix.begin();
    const matrix_iterator this_end = m_matrix.end();
    matrix_iterator ithat =          rhs.m_matrix.begin();
    const matrix_iterator that_end = rhs.m_matrix.end();
    for (; ithis != this_end && ithat != that_end; ++ithis, ++ithat)
    {
      typedef typename matrix_t::Row::const_iterator row_iterator;
      row_iterator jthis =           (*ithis).begin();
      const row_iterator jthis_end = (*ithis).end();
      row_iterator jthat =           (*ithat).begin();
      const row_iterator jthat_end = (*ithat).end();
      for (; jthis != jthis_end && jthat != jthat_end; ++jthis, ++jthat)
        if (*jthis != *jthat)
          return false;
      if (jthis != jthis_end || jthat != jthat_end)
        return false;
    }
    return ithis == this_end && ithat == that_end;
  }

  // Test for equality of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const bool
  matrix_multi<Scalar_T,LO,HI>::
  operator==  (const Scalar_T& scr) const
  { return (*this) == multivector_t(framed_multi_t(scr), m_frame, true); }

  /// Geometric sum of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator+= (const Scalar_T& scr)
  {
    return *this += pair_t(index_set_t(), scr);
  }

  /// Geometric sum
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator+= (const multivector_t& rhs)
  {
    // Operate only within a common frame
    const index_set_t our_frame = m_frame | rhs.m_frame;
    if (m_frame != our_frame)
      // Represent *this in our_frame via conversion through framed_multi_t
      *this = multivector_t(framed_multi_t(*this), our_frame, true);
    if (rhs.m_frame == our_frame)
      mtl::add(rhs.m_matrix, m_matrix);
    else
    { // Represent rhs in our_frame via conversion through framed_multi_t
      const multivector_t our_rhs = multivector_t(framed_multi_t(rhs), our_frame, true);
      mtl::add(our_rhs.m_matrix, m_matrix);
    }
    return *this;
  }

  /// Geometric difference
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator-= (const multivector_t& rhs)
  {
    // Operate only within a common frame
    const index_set_t our_frame = m_frame | rhs.m_frame;
    if (m_frame != our_frame)
      // Represent *this in our_frame via conversion through framed_multi_t
      *this = multivector_t(framed_multi_t(*this), our_frame, true);
    if (rhs.m_frame == our_frame)
      mtl::add(scaled(rhs.m_matrix, Scalar_T(-1)), m_matrix);
    else
    { // Represent rhs in our_frame via conversion through framed_multi_t
      const multivector_t our_rhs = multivector_t(framed_multi_t(rhs), our_frame, true);
      mtl::add(scaled(our_rhs.m_matrix, Scalar_T(-1)), m_matrix);
    }
    return *this;
  }

  /// Unary -
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const
  matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  operator- () const
  {
    multivector_t result(m_matrix, m_frame);
    mtl::scale(result.m_matrix, Scalar_T(-1));
    return result;
  }

  /// Product of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator*= (const Scalar_T& scr)
  { // multiply coordinates of all terms by scalar
    if (scr == Scalar_T(0))
      *this = 0;
    else
      mtl::scale(this->m_matrix, scr);
    return *this;
  }

  /// Geometric product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator*= (const multivector_t& rhs)
  { // MTL quirk: Need to have enough workspace for sparse multiply
    // Work around: Use dense result when necessary, which is almost always
    const bool dense_mult = !(is_perm_shaped(m_matrix) && is_perm_shaped(rhs.m_matrix));
    // Operate only within a common frame
    const index_set_t our_frame = m_frame | rhs.m_frame;
    if (m_frame != our_frame)
      // Represent *this in our_frame via conversion through framed_multi_t
      *this = multivector_t(framed_multi_t(*this), our_frame, true);
    if (dense_mult)
    {
      dense_matrix_t result_matrix( m_matrix.nrows(), m_matrix.ncols());
      if (rhs.m_frame == our_frame)
      {
        mtl::mult(m_matrix, rhs.m_matrix, result_matrix);
        mtl::copy(result_matrix, m_matrix);
      }
      else
      { // Represent rhs in our_frame via conversion through framed_multi_t
        const multivector_t our_rhs = multivector_t(framed_multi_t(rhs), our_frame, true);
        mtl::mult(m_matrix, our_rhs.m_matrix, result_matrix);
        mtl::copy(result_matrix, m_matrix);
      }
    }
    else
    {
      matrix_t result_matrix( m_matrix.nrows(), m_matrix.ncols());
      if (rhs.m_frame == our_frame)
      {
        mtl::mult(m_matrix, rhs.m_matrix, result_matrix);
        mtl::copy(result_matrix, m_matrix);
      }
      else
      { // Represent rhs in our_frame via conversion through framed_multi_t
        const multivector_t our_rhs = multivector_t(framed_multi_t(rhs), our_frame, true);
        mtl::mult(m_matrix, our_rhs.m_matrix, result_matrix);
        mtl::copy(result_matrix, m_matrix);
      }
    }
    return *this;
  }

  /// Contraction
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator%= (const multivector_t& rhs)
  {
    framed_multi_t lhs(*this);
    return *this = lhs %= framed_multi_t(rhs);
  }

  /// Inner product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator&= (const multivector_t& rhs)
  { // Arvind Raja's original reference:
    // "old clical, innerproduct(p,q:pterm):pterm in file compmod.pas"
    framed_multi_t lhs(*this);
    return *this = lhs &= framed_multi_t(rhs);
  }

  /// Outer product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator^= (const multivector_t& rhs)
  { // Arvind Raja's original reference:
    // "old clical, outerproduct(p,q:pterm):pterm in file compmod.pas"
    framed_multi_t lhs(*this);
    return *this = lhs ^= framed_multi_t(rhs);
  }

  /// Quotient of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator/= (const Scalar_T& scr)
  {
    return *this *= Scalar_T(1)/scr;
  }

  // Geometric quotient
  template< typename Scalar_T, const index_t LO, const index_t HI >
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator/= (const multivector_t& rhs)
  {
    if( m_frame != rhs.m_frame )
      return *this = framed_multi_t(*this) /
                     framed_multi_t(rhs);
    else
    { // Solve result == *this/rhs <=> result*rhs == *this
      // Define references so we now solve X == B/A
      // (where X == result, B == m_matrix and A == rhs.m_matrix)
      // X == B/A <=> X*A == B <=> AT*XT == BT

      // We solve AT*XT == BT one column at a time as
      // AT*XT(:,i) == BT(:,i), remembering that
      // columns of BT are rows of B and
      // columns of XT are rows of X, so we solve
      // AT*X(i,:) == B(i,:), ie. solve
      // AT*x == B(i,:); X(i,:) = x;

      typedef mtl::dense1D<Scalar_T> vector_t;
      matrix_t AT(rhs.m_matrix.ncols(), rhs.m_matrix.nrows());
      mtl::transpose(rhs.m_matrix, AT);
      dense_matrix_t LU(AT.nrows(), AT.ncols());
      mtl::copy(AT, LU);
      mtl::dense1D<int> pvector(AT.nrows());
      if (mtl::lu_factor(LU, pvector) != 0)
        // AT is singular. Return IEEE NaN
        return *this = Scalar_T(std::log(0.0)); // This actually returns -Inf
      for (typename matrix_t::iterator i = m_matrix.begin(); i != m_matrix.end(); ++i)
      { // Solve AT*x == b, where b == B(i,:)
        vector_t x(AT.ncols(), Scalar_T(0));
        mtl::lu_solve(LU, pvector, *i, x);
        // Iterative refinement.
        // Reference: Nicholas J. Higham, "Accuracy and Stability of Numerical Algorithms",
        // SIAM, 1996, ISBN 0-89871-355-2, Chapter 11
        if (Tune_P::iterative_refinement_max_steps > 0)
        {
          vector_t mb(AT.nrows());
          mtl::copy(scaled(*i, Scalar_T(-1)), mb);
          vector_t r(AT.nrows(), Scalar_T(0));
          // r = AT*x - b;
          mtl::mult(AT, x, mb, r);
          // nr = |r|;
          Scalar_T nr = mtl::infinity_norm(r);
          // (nr == nr) below is to exclude NaN
          if ( nr != Scalar_T(0) && nr == nr )
          {
            vector_t xnew(AT.ncols());
            // xnew = x;
            mtl::copy(x, xnew);
            Scalar_T nrold = nr + Scalar_T(1);
            for (int step = 0;
                 step != Tune_P::iterative_refinement_max_steps &&
                 nr < nrold &&
                 nr != Scalar_T(0) &&
                 nr == nr;
                  ++step )
            {
              nrold = nr;
              if (step != 0)
                // x = xnew;
                mtl::copy(xnew, x);
              vector_t d(AT.ncols());
              // d = r / A;
              mtl::lu_solve(LU, pvector, r, d);
              // xnew -= d;
              mtl::add(scaled(d, Scalar_T(-1)), xnew);
              // r = A*xnew - b;
              mtl::mult(AT, xnew, mb, r);
              // nr = |r|;
              nr = mtl::infinity_norm(r);
            }
          }
        }
        // Set X(i,:) = x
        const matrix_index_t i_index = i.index();
        mtl::set(*i, 0);
        for (vector_t::const_iterator j = x.begin(); j != x.end(); ++j)
          if (*j !=  Scalar_T(0))
            m_matrix(i_index, j.index()) = *j;
      }
      return *this;
    }
  }

  /// Clifford multiplicative inverse
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const
  matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  inv() const
  {
    multivector_t result(1, m_frame);
    return result /= *this;
  }

  /// Subscripting: map from index set to scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const Scalar_T
  matrix_multi<Scalar_T,LO,HI>::
  operator[] (const index_set_t& ist) const
  {
    matrix_t e;
    basis_element<Scalar_T,LO,HI>(ist, m_frame, e);
    return inner<matrix_t,Scalar_T>(e, m_matrix);
  }

  /// Main involution, each {i} is replaced by -{i} in each term
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  involute() const
  { return framed_multi_t(*this).involute(); }

  /// Reversion, order of {i} is reversed in each term
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  reverse() const
  { return framed_multi_t(*this).reverse(); }

  /// Conjugation, conj == reverse o involute == involute o reverse
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  conj() const
  { return framed_multi_t(*this).conj(); }

  /// Quadratic form := scalar part of rev(x)*x
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const Scalar_T
  matrix_multi<Scalar_T,LO,HI>::
  quad() const
  { // scalar(conj(x)*x) = 2*quad(even(x)) - quad(x)
    // Arvind Raja ref: "old clical: quadfunction(p:pter):pterm in file compmod.pas"
    return framed_multi_t(*this).quad();
  }

  /// Scalar_T norm squared= sum of norm squared of coordinates
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const Scalar_T
  matrix_multi<Scalar_T,LO,HI>::
  norm() const
  { return inner<matrix_t,Scalar_T>(m_matrix, m_matrix); }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  pow(int m) const
  {
    multivector_t a;
    if (m < 0)
    {
      m = -m;
      a = (*this).inv();
    }
    else
      a = *this;
    multivector_t result = 1;
    for (; m != 0; m >>= 1)
    {
      if (m & 1)
        result *= a;
      a *= a;
    }
    return result;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  outer_pow(int m) const
  { // outer product power
    if (m < 0)
      throw error_t("outer_pow(m): negative exponent");
    multivector_t result = 1;
    multivector_t a = *this;
    for (; m != 0; m >>= 1)
    {
      if (m & 1)
        result ^= a;
      a ^= a;
    }
    return result;
  }

  /// Grading: part where each term is a grade-vector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  operator() (index_t grade) const
  {
    if ((grade < 0) || (grade > HI-LO))
  		return 0;
    else
      return (framed_multi_t(*this))(grade);
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  even() const
  { // even part of x, sum of the pure(count) with even count
    return framed_multi_t(*this).even();
  }

  /// Vector part of multivector, as a vector_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  matrix_multi<Scalar_T,LO,HI>::vector_t
  matrix_multi<Scalar_T,LO,HI>::
  vector_part() const
  {
    vector_t result;
    const index_t begin_index = m_frame.min();
    const index_t end_index = m_frame.max()+1;
    for (index_t idx = begin_index; idx != end_index; ++idx)
      if (m_frame[idx])
      {
        // Frame may contain indices which do not correspond to a grade 1 term but
        // frame cannot omit any index corresponding to a grade 1 term
        matrix_t e;
        basis_element<Scalar_T,LO,HI>(index_set_t(idx), m_frame, e);
        result.push_back(inner<matrix_t,Scalar_T>(e, m_matrix));
      }
    return result;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  ostream&
  operator<< (ostream& os, const matrix_multi<Scalar_T,LO,HI>& val)
  {
    os << matrix_multi<Scalar_T,LO,HI>::framed_multi_t(val);
    return os;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  istream&
  operator>> (istream& s, matrix_multi<Scalar_T,LO,HI>& val)
  { // Input looks like 1.0-2.0{1,2}+3.2{3,4}
    framed_multi<Scalar_T,LO,HI> local;
    s >> local;
    // If s.bad() then we have a corrupt input
    // otherwise we are fine and can copy the resulting matrix_multi
    if (!s.bad())
      val = local;
    return s;
  }

  /// Write out multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  void
  matrix_multi<Scalar_T,LO,HI>::
  write(char* msg) const
  { framed_multi_t(*this).write(msg); }

  /// Write out multivector to file
  template< typename Scalar_T, const index_t LO, const index_t HI >
  void
  matrix_multi<Scalar_T,LO,HI>::
  write(ofstream& ofile, char* msg) const
  {
    if (!ofile)
      throw error_t("write(ofile,msg): cannot write to output file");
    framed_multi_t(*this).write(ofile, msg);
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  Scalar_T
  matrix_multi<Scalar_T,LO,HI>::
  max_abs() const
  { return framed_multi_t(*this).max_abs(); }

  /// Check if a multivector contains any IEEE NaN values
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const bool
  matrix_multi<Scalar_T,LO,HI>::
  isnan() const
  { // The distinguishing feature is that NaN != NaN
    for (typename matrix_t::const_iterator i = m_matrix.begin(); i != m_matrix.end(); ++i)
      for (typename matrix_t::Row::const_iterator j = (*i).begin(); j != (*i).end(); ++j)
        if (*j != *j)
          return true;
    return false;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  matrix_multi<Scalar_T,LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  truncated(const Scalar_T& limit) const
  { return framed_multi_t(*this).truncated(limit); }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  const index_set<LO,HI>
  matrix_multi<Scalar_T,LO,HI>::
  frame() const
  { return m_frame; }

  /// Add a term, if non-zero
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  matrix_multi<Scalar_T,LO,HI>&
  matrix_multi<Scalar_T,LO,HI>::
  operator+= (const pair_t& term)
  {
    if (term.second != Scalar_T(0))
    {
      matrix_t e;
      basis_element<Scalar_T,LO,HI>(term.first, m_frame, e);
      mtl::add( scaled(e, term.second), m_matrix );
    }
    return *this;
  }

  /// Create a basis element matrix within a frame
  template< typename Scalar_T, const index_t LO, const index_t HI >
  void
  basis_element( const index_set<LO,HI>& ist, const index_set<LO,HI>& m_frame,
                 matrix_multi<Scalar_T,LO,HI>::matrix_t& result )
  {
    typedef matrix_multi<Scalar_T,LO,HI>  multivector_t;
    typedef multivector_t::matrix_t       matrix_t;
    typedef multivector_t::matrix_index_t matrix_index_t;
    typedef index_set<LO,HI>              index_set_t;

    static vector< vector< matrix_t > > generators;

    const index_set_t folded_set = ist.fold(m_frame);
    const index_set_t folded_frame = m_frame.fold();
    const index_t n = max(index_t(-folded_frame.min()), folded_frame.max());
    const matrix_index_t dim = 1 << n;
    result = matrix_t(dim, dim);
    unit(dim, result);
    const vector<matrix_t>& e = gengen( n, generators );
    for (index_t k = max(LO,index_t(-n)); k <= min(n,HI); ++k )
      if ( folded_set[k] )
      {
        matrix_t y( dim, dim );
        mtl::mult(result, e[k+n], y);
        mtl::copy(y, result);
      }
  }

  /// Determine the matrix dimension of the fold of a subalegbra
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  matrix_multi<Scalar_T,LO,HI>::matrix_index_t
  folded_dim( const index_set<LO,HI>& sub )
  {
    const index_set<LO,HI> folded = sub.fold();
    return 1 << max(index_t(-folded.min()), folded.max());
  }
}
#endif  // _GLUCAT_MATRIX_MULTI_IMP_H
