#ifndef _GLUCAT_INDEX_SET_H
#define _GLUCAT_INDEX_SET_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    index_set.h : Declare a class for a set of non-zero integer indices
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2012 by Paul C. Leopardi
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

#include "glucat/global.h"
#include "glucat/errors.h"

#include <boost/static_assert.hpp>

#include <bitset>
#include <utility>

namespace glucat
{
  template<const index_t LO, const index_t HI>
  class index_set; // forward

  /// Symmetric set difference: exclusive or
  template<const index_t LO, const index_t HI>
  const index_set<LO,HI>
  operator^ (const index_set<LO,HI>& lhs,
             const index_set<LO,HI>& rhs);

  /// Set intersection: and
  template<const index_t LO, const index_t HI>
  const index_set<LO,HI>
  operator& (const index_set<LO,HI>& lhs,
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

  /// Index set class based on std::bitset<> in Gnu standard C++ library
  template<const index_t LO, const index_t HI>
  class index_set :
  private std::bitset<HI-LO>
  {
  private:
    BOOST_STATIC_ASSERT((LO <= 0) && (0 <= HI) && (LO < HI) && \
                       (-LO    <  _GLUCAT_BITS_PER_ULONG)   && \
                       ( HI    <  _GLUCAT_BITS_PER_ULONG)   && \
                       ( HI-LO <= _GLUCAT_BITS_PER_ULONG));
    typedef std::bitset<HI-LO>          bitset_t;
    typedef error<index_set>            error_t;
  public:
    typedef index_set                   index_set_t;
    typedef std::pair<index_t,index_t>  index_pair_t;

    static const index_t v_lo = LO;
    static const index_t v_hi = HI;

    static const std::string  classname();
    /// Default constructor creates an empty set
    index_set    () = default;
    /// Constructor from bitset_t
    index_set    (const bitset_t bst);
    /// Constructor from index
    index_set    (const index_t idx);
    /// Constructor from set value of an index set folded within the given frame
    index_set    (const set_value_t folded_val, const index_set_t frm, const bool prechecked = false);
    /// Constructor from range of indices from range.first to range.second
    index_set    (const index_pair_t& range, const bool prechecked = false);
    /// Constructor from string
    index_set    (const std::string& str);

    /// Equality
    bool         operator==  (const index_set_t rhs) const;
    /// Inequality
    bool         operator!=  (const index_set_t rhs) const;
    /// Set complement: not
    index_set_t  operator~   () const;
    /// Symmetric set difference: exclusive or
    index_set_t& operator^=  (const index_set_t rhs);
    /// Set intersection: and
    index_set_t& operator&=  (const index_set_t rhs);
    /// Set union: or
    index_set_t& operator|=  (const index_set_t rhs);
    /// Subscripting: Test idx for membership: test value of bit idx
    bool         operator[]  (const index_t idx) const;
    /// Test idx for membership: test value of bit idx
    bool         test(const index_t idx) const;
    /// Include all indices except 0: set all bits except 0
    index_set_t& set();
    /// Include idx: Set bit at idx if idx != 0
    index_set_t& set(const index_t idx);
    /// Set membership of idx to val if idx != 0: Set bit at idx to val if idx != 0
    index_set_t& set(const index_t idx, const int val);
    /// Make set empty: Set all bits to 0
    index_set_t& reset();
    /// Exclude idx:  Set bit at idx to 0
    index_set_t& reset(const index_t idx);
    /// Set complement, except 0: flip all bits, except 0
    index_set_t& flip();
    /// Complement membership of idx if idx != 0: flip bit at idx if idx != 0
    index_set_t& flip(const index_t idx);
    /// Cardinality: Number of indices included in set
    index_t      count() const;
    /// Number of negative indices included in set
    index_t      count_neg() const;
    /// Number of positive indices included in set
    index_t      count_pos() const;
    /// Minimum member
    index_t      min() const;
    /// Maximum member
    index_t      max() const;

  // Functions which support Clifford algebra operations
    /// Less than operator used for comparisons, map, etc.
    bool                  operator<     (const index_set_t rhs) const;
    /// Determine if the index set is contiguous, ie. has no gaps
    bool                  is_contiguous () const;
    /// Fold this index set within itself as a frame
    const index_set_t     fold          () const;
    /// Fold this index set within the given frame
    const index_set_t     fold          (const index_set_t frm, const bool prechecked = false) const;
    /// Unfold this index set within the given frame
    const index_set_t     unfold        (const index_set_t frm, const bool prechecked = false) const;
    /// The set value of the fold of this index set within the given frame
    set_value_t           value_of_fold (const index_set_t frm) const;
    /// Sign of geometric product of two Clifford basis elements
    int                   sign_of_mult  (const index_set_t ist) const;
    /// Sign of geometric square of a Clifford basis element
    int                   sign_of_square()                      const;

    /// Hash function
    size_t                hash_fn       ()                      const;

  // Friends
    friend const index_set_t operator^<> (const index_set_t& lhs, const index_set_t& rhs);
    friend const index_set_t operator&<> (const index_set_t& lhs, const index_set_t& rhs);
    friend const index_set_t operator|<> (const index_set_t& lhs, const index_set_t& rhs);
    friend int compare<>                 (const index_set_t& lhs, const index_set_t& rhs);

  // Member reference:
    class reference;
    friend class reference;

    /// Index set member reference
    class reference {
      friend class index_set;

    public:
      /// Default constructor is deleted
      reference() = delete;
      reference   (index_set_t& ist, index_t idx);
      ~reference  () = default;
      /// for b[i] = x;
      reference&  operator= (const bool x);
      /// for b[i] = b[j];
      reference&  operator= (const reference& j);
      /// Flips a bit
      bool        operator~ () const;
      /// for x = b[i];
                  operator bool () const;
      /// for b[i].flip();
      reference&  flip();

    private:
      index_set_t* m_pst;
      index_t      m_idx;
    };
    /// Subscripting: Element access
    reference     operator[](index_t idx);
  private:
    /// Lexicographic ordering of two sets: *this < rhs
    bool          lex_less_than (const index_set_t rhs) const;
  };

  /// Size of set_value_t should be enough to contain bitset<DEFAULT_HI-DEFAULT_LO>
  _GLUCAT_CTAssert(sizeof(set_value_t) >= sizeof(std::bitset<DEFAULT_HI-DEFAULT_LO>),
           Default_index_set_too_big_for_value)

  // non-members

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

  /// Minimum negative index, or 0 if none
  template<const index_t LO, const index_t HI>
  index_t
  min_neg(const index_set<LO,HI>& ist);

  /// Maximum positive index, or 0 if none
  template<const index_t LO, const index_t HI>
  index_t
  max_pos(const index_set<LO,HI>& ist);
}
#endif // _GLUCAT_INDEX_SET_H
