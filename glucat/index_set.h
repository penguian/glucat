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
  auto
  operator^ (const index_set<LO,HI>& lhs,
             const index_set<LO,HI>& rhs) -> const index_set<LO,HI>;

  /// Set intersection: and
  template<const index_t LO, const index_t HI>
  auto
  operator& (const index_set<LO,HI>& lhs,
             const index_set<LO,HI>& rhs) -> const index_set<LO,HI>;

  /// Set union: or
  template<const index_t LO, const index_t HI>
  auto
  operator| (const index_set<LO,HI>& lhs,
             const index_set<LO,HI>& rhs) -> const index_set<LO,HI>;

  /// "lexicographic compare" eg. {3,4,5} is less than {3,7,8}
  // -1 if a<b, +1 if a>b, 0 if a==b
  template<const index_t LO, const index_t HI>
  auto
  compare(const index_set<LO,HI>& a, const index_set<LO,HI>& b) -> int;

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
    using bitset_t = std::bitset<HI - LO>;
    using error_t = error<index_set>;
  public:
    using index_set_t = index_set;
    using index_pair_t = std::pair<index_t, index_t>;

    static const index_t v_lo = LO;
    static const index_t v_hi = HI;

    static auto classname() -> const std::string;
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

    // Default copy constructor and assignment
    index_set    (const index_set&) = default;
    auto operator= (const index_set&) -> index_set& = default;

    // Default move constructor and assignment
    index_set    (index_set&&) = default;
    auto operator= (index_set&&) -> index_set& = default;

    /// Equality
    auto operator==  (const index_set_t rhs) const -> bool;
    /// Inequality
    auto operator!=  (const index_set_t rhs) const -> bool;
    /// Set complement: not
    auto operator~   () const -> index_set_t;
    /// Symmetric set difference: exclusive or
    auto operator^=  (const index_set_t rhs) -> index_set_t&;
    /// Set intersection: and
    auto operator&=  (const index_set_t rhs) -> index_set_t&;
    /// Set union: or
    auto operator|=  (const index_set_t rhs) -> index_set_t&;
    /// Subscripting: Test idx for membership: test value of bit idx
    auto operator[]  (const index_t idx) const -> bool;
    /// Test idx for membership: test value of bit idx
    auto test(const index_t idx) const -> bool;
    /// Include all indices except 0: set all bits except 0
    auto set() -> index_set_t&;
    /// Include idx: Set bit at idx if idx != 0
    auto set(const index_t idx) -> index_set_t&;
    /// Set membership of idx to val if idx != 0: Set bit at idx to val if idx != 0
    auto set(const index_t idx, const int val) -> index_set_t&;
    /// Make set empty: Set all bits to 0
    auto reset() -> index_set_t&;
    /// Exclude idx:  Set bit at idx to 0
    auto reset(const index_t idx) -> index_set_t&;
    /// Set complement, except 0: flip all bits, except 0
    auto flip() -> index_set_t&;
    /// Complement membership of idx if idx != 0: flip bit at idx if idx != 0
    auto flip(const index_t idx) -> index_set_t&;
    /// Cardinality: Number of indices included in set
    auto count() const -> index_t;
    /// Number of negative indices included in set
    auto count_neg() const -> index_t;
    /// Number of positive indices included in set
    auto count_pos() const -> index_t;
    /// Minimum member
    auto min() const -> index_t;
    /// Maximum member
    auto max() const -> index_t;

  // Functions which support Clifford algebra operations
    /// Less than operator used for comparisons, map, etc.
    auto operator<     (const index_set_t rhs) const -> bool;
    /// Determine if the index set is contiguous, ie. has no gaps
    auto is_contiguous () const -> bool;
    /// Fold this index set within itself as a frame
    auto fold          () const -> const index_set_t;
    /// Fold this index set within the given frame
    auto fold          (const index_set_t frm, const bool prechecked = false) const -> const index_set_t;
    /// Unfold this index set within the given frame
    auto unfold        (const index_set_t frm, const bool prechecked = false) const -> const index_set_t;
    /// The set value of the fold of this index set within the given frame
    auto value_of_fold (const index_set_t frm) const -> set_value_t;
    /// Sign of geometric product of two Clifford basis elements
    auto sign_of_mult  (const index_set_t ist) const -> int;
    /// Sign of geometric square of a Clifford basis element
    auto sign_of_square()                      const -> int;

    /// Hash function
    auto hash_fn       ()                      const -> size_t;

  // Friends
    friend auto operator^<> (const index_set_t& lhs, const index_set_t& rhs) -> const index_set_t;
    friend auto operator&<> (const index_set_t& lhs, const index_set_t& rhs) -> const index_set_t;
    friend auto operator|<> (const index_set_t& lhs, const index_set_t& rhs) -> const index_set_t;
    friend auto compare<>                 (const index_set_t& lhs, const index_set_t& rhs) -> int;

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
      /// for b[i] == c[j];
      auto operator== (const reference& c_j) const -> bool;
      /// for b[i] = x;
      auto operator= (const bool x) -> reference&;
      /// for b[i] = c[j];
      auto operator= (const reference& c_j) -> reference&;
      /// Flips a bit
      auto operator~ () const -> bool;
      /// for x = b[i];
           operator bool () const;
      /// for b[i].flip();
      auto flip() -> reference&;

    private:
      index_set_t* m_pst;
      index_t      m_idx;
    };
    /// Subscripting: Element access
    auto operator[](index_t idx) -> reference;
  private:
    /// Lexicographic ordering of two sets: *this < rhs
    auto lex_less_than (const index_set_t rhs) const -> bool;
  };

  /// Size of set_value_t should be enough to contain bitset<DEFAULT_HI-DEFAULT_LO>
  _GLUCAT_CTAssert(sizeof(set_value_t) >= sizeof(std::bitset<DEFAULT_HI-DEFAULT_LO>),
           Default_index_set_too_big_for_value)

  // non-members

  /// Write out index set
  template<const index_t LO, const index_t HI>
  auto
  operator<< (std::ostream& os, const index_set<LO,HI>& ist) -> std::ostream&;

  /// Read in index set
  template<const index_t LO, const index_t HI>
  auto
  operator>> (std::istream& s, index_set<LO,HI>& ist) -> std::istream&;

  // Functions which support Clifford algebra operations
  /// Square of generator {j}
  auto sign_of_square(index_t j) -> int;

  /// Minimum negative index, or 0 if none
  template<const index_t LO, const index_t HI>
  auto
  min_neg(const index_set<LO,HI>& ist) -> index_t;

  /// Maximum positive index, or 0 if none
  template<const index_t LO, const index_t HI>
  auto
  max_pos(const index_set<LO,HI>& ist) -> index_t;
}
#endif // _GLUCAT_INDEX_SET_H
