#ifndef _GLUCAT_INDEX_SET_H
#define _GLUCAT_INDEX_SET_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    index_set.h : Declare a class for a set of non-zero integer indices
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
  /// Index set class based on std::bitset<> in Gnu standard C++ library
  template<const index_t LO, const index_t HI>
  class index_set :
  private std::bitset<HI-LO+1>
  {
  private:
    typedef std::bitset<HI-LO+1>  bitset_t;
    typedef error<index_set>      error_t;
  public:
    typedef index_set             index_set_t;
    static const index_t          v_lo = LO;
    static const index_t          v_hi = HI;
    
    static const std::string  classname();
    /// Default constructor creates an empty set
    index_set     () { }
    /// Constructor from index
    index_set     (const index_t& idx);
    /// Constructor from set value of an index set folded within the given frame
    index_set     (const set_value_t& folded_val, const index_set& frm, const bool prechecked = false);
    /// Constructor from string
    index_set     (const std::string& str);

    /// Equality
    bool          operator==  (const index_set& rhs) const;
    /// Inequality
    bool          operator!=  (const index_set& rhs) const;
    /// Symmetric set difference: exclusive or
    index_set&     operator^=  (const index_set& rhs);
    /// Set union: or
    index_set&     operator|=  (const index_set& rhs);
    /// Subscripting: Test idx for membership: test value of bit idx
    bool          operator[]  (index_t idx) const;
    /// Test idx for membership: test value of bit idx
    bool          test(index_t idx) const;
    /// Include all indices except 0: set all bits except 0
    index_set&    set();
    /// Include idx: Set bit at idx if idx != 0
    index_set&    set(index_t idx);
    /// Set membership of idx to val if idx != 0: Set bit at idx to val if idx != 0
    index_set&    set(index_t idx, int val);
    /// Make set empty: Set all bits to 0
    index_set&    reset();
    /// Exclude idx:  Set bit at idx to 0
    index_set&    reset(index_t idx);
    /// Set complement, except 0: flip all bits, except 0
    index_set&    flip();
    /// Complement membership of idx if idx != 0: flip bit at idx if idx != 0
    index_set&    flip(index_t idx);
    /// Cardinality: Number of indices included in set
    index_t       count() const;
    /// Number of positive indices included in set
    index_t       count_pos() const;
    /// Number of negative indices included in set
    index_t       count_neg() const;
    /// Minimum member
    index_t       min() const;
    /// Maximum member
    index_t       max() const;

  // Functions which support Clifford algebra operations
    /// Fold this index set within itself as a frame
    const index_set     fold          () const;
    /// Fold this index set within the given frame
    const index_set     fold          (const index_set& frm, const bool prechecked = false) const;
    /// Unfold this index set within the given frame
    const index_set     unfold        (const index_set& frm, const bool prechecked = false) const;
    /// The set value of the fold of this index set within the given frame
    set_value_t   value_of_fold (const index_set& frm) const;
    int           sign_of_mult  (const index_set& ist) const;

  // Member reference:
    class reference;
    friend class  reference;

    /// Index set member reference
    class reference {
      friend class index_set;

      /// Private default constructor is left undefined
      reference();
    public:
      reference   ( index_set& ist, index_t idx );
      ~reference  () {}
      /// for b[i] = x;
      reference&  operator= (bool x);
      /// for b[i] = b[j];
      reference&  operator= (const reference& j);
      /// Flips a bit
      bool        operator~ () const;
      /// for x = b[i];
                  operator bool () const;
      /// for b[i].flip();
      reference& flip();

    private:
      index_set*  m_pst;
      index_t     m_idx;
    };
    /// Subscripting: Element access
    reference     operator[](index_t idx);
  };

  /// Size of set_value_t should be enough to contain bitset<DEFAULT_HI-DEFAULT_LO+1>
  _GLUCAT_CTAssert(sizeof(set_value_t) >= sizeof(std::bitset<DEFAULT_HI-DEFAULT_LO+1>),
           Default_index_set_too_big_for_value)

  // non-members
  /// Symmetric set difference: exclusive or
  template<const index_t LO, const index_t HI>
  const index_set<LO,HI>
  operator^ (const index_set<LO,HI>& lhs,
             const index_set<LO,HI>& rhs);

  /// Set union: or
  template<const index_t LO, const index_t HI>
  const index_set<LO,HI>
  operator| (const index_set<LO,HI>& lhs,
             const index_set<LO,HI>& rhs);

  /// "lexicographic compare" eg. {3,4,5} is less than {3,7,8}
  // -1 if a<b, +1 if a>b, 0 if a==b
  template<const index_t LO, const index_t HI>
  int
  compare(const index_set<LO,HI>& a, const index_set<LO,HI>& b);

  /// Less than operator used for comparisons, map, etc.
  template<const index_t LO, const index_t HI>
  bool
  operator< (const index_set<LO,HI>& a, const index_set<LO,HI>& b);

  /// Write out index set
  template<const index_t LO, const index_t HI>
  std::ostream&
  operator<< (std::ostream& os, const index_set<LO,HI>& ist);

  /// Read in index set
  template<const index_t LO, const index_t HI>
  std::istream&
  operator>> (std::istream& s, index_set<LO,HI>& ist);

  // Functions which support Clifford algebra operations
  /// Square of generator {j}
  int     sign_of_square(index_t j);

  /// maximum positive index, or 0 if none
  template<const index_t LO, const index_t HI>
  index_t
  max_pos(const index_set<LO,HI>& ist);

  /// minimum negative index, or 0 if none
  template<const index_t LO, const index_t HI>
  index_t
  min_neg(const index_set<LO,HI>& ist);

  /// Set containing a range of indices from range_min to range_max
  template<const index_t LO, const index_t HI>
  const index_set<LO,HI>
  index_range(const index_t range_min, const index_t range_max);
}
#endif // _GLUCAT_INDEX_SET_H
