#ifndef _GLUCAT_FRAMED_MULTI_H
#define _GLUCAT_FRAMED_MULTI_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    framed_multi.h : Declare a class for the framed representation of a multivector
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
  template< typename Scalar_T, const index_t LO, const index_t HI >
  class matrix_multi; // forward

  /// A framed_multi<Scalar_T,LO,HI> is a framed approximation to a multivector
  template< typename Scalar_T, const index_t LO = DEFAULT_LO, const index_t HI = DEFAULT_HI >
  class framed_multi :
  public clifford_algebra< Scalar_T, index_set<LO,HI>, framed_multi<Scalar_T,LO,HI> >,
  private std::map< const index_set<LO,HI>, Scalar_T >
  {
  public:
    typedef framed_multi                  multivector_t;
    typedef Scalar_T                      scalar_t;
    typedef index_set<LO,HI>              index_set_t;
    typedef std::pair< const index_set_t, Scalar_T > pair_t;
    typedef std::vector<Scalar_T>         vector_t;
    typedef error<multivector_t>          error_t;
    typedef matrix_multi<Scalar_T,LO,HI>  matrix_multi_t;
    typedef multivector_t                 framed_multi_t;
    friend class matrix_multi_t;
  private:
    typedef std::map< const index_set_t, Scalar_T > map_t;
    typedef typename map_t::iterator       iterator;
    typedef typename map_t::const_iterator const_iterator;
    typedef typename matrix_multi_t::matrix_t matrix_t;
  public:
    /// Class name used in messages
    static const char* classname();
		/// Destructor
		~framed_multi() {};
    /// Default constructor
    framed_multi();
  	/// Construct a multivector, within a given frame, from a given multivector
    framed_multi(const multivector_t& val,
								 const index_set_t& frm, const bool prechecked = false);
    /// Construct a multivector from an index set and a scalar coordinate
    framed_multi(const index_set_t& ist, const Scalar_T& crd);
  	/// Construct a multivector, within a given frame, from an index set and a scalar coordinate
    framed_multi(const index_set_t& ist, const Scalar_T& crd,
								 const index_set_t& frm, const bool prechecked = false);
    /// Construct a multivector from a scalar (within a frame, if given)
    framed_multi(const Scalar_T& scr, const index_set_t& frm = index_set_t());
    /// Construct a multivector from an int (within a frame, if given)
    framed_multi(const int scr, const index_set_t& frm = index_set_t());
  	/// Construct a multivector, within a given frame, from a given vector
    framed_multi(const vector_t& vec,
								 const index_set_t& frm, const bool prechecked = false);
    /// Construct a multivector from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
    framed_multi(const std::string& str);
    /// Construct a multivector, within a given frame, from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
    framed_multi(const std::string& str,
								 const index_set_t& frm, const bool prechecked = false);
    /// Construct a multivector from a matrix_multi_t
    framed_multi        (const matrix_multi_t& val);
    _GLUCAT_CLIFFORD_ALGEBRA_OPERATIONS
    friend std::istream&
      operator>> <>(std::istream& s, multivector_t& val);
    friend std::ostream&
      operator<< <>(std::ostream& os, const multivector_t& val);
    friend std::ostream&
      operator<< <>(std::ostream& os, const pair_t& term);
  private:
    /// Add a term, if non-zero
    multivector_t&      operator+= (const pair_t& term);
  };
  // non-members

  /// Read multivector from input
  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::istream&
  operator>> (std::istream& s, framed_multi<Scalar_T,LO,HI>& val);

  /// Write multivector to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::ostream&
  operator<< (std::ostream& os, const framed_multi<Scalar_T,LO,HI>& val);

  /// Write term to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::ostream&
  operator<< (std::ostream& os, const std::pair< const index_set<LO,HI>, Scalar_T >& term);

  /// Product of terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const std::pair<const index_set<LO,HI>, Scalar_T>
  operator*
   (const std::pair<const index_set<LO,HI>, Scalar_T>& lhs,
    const std::pair<const index_set<LO,HI>, Scalar_T>& rhs);
}
#endif  // _GLUCAT_FRAMED_MULTI_H
