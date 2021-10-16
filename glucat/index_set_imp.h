#ifndef _GLUCAT_INDEX_SET_IMP_H
#define _GLUCAT_INDEX_SET_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    index_set_imp.h : Implement a class for a set of non-zero integer indices
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

#include "glucat/index_set.h"

#include <sstream>

namespace glucat
{
  // References for algorithms:
  // [JA]: Joerg Arndt,  "Algorithms for programmers", http://www.jjj.de/fxt/fxtbook.pdf
  //      Chapter 1, Bit wizardry, http://www.jjj.de/bitwizardry/bitwizardrypage.html
  // [L]: Pertti Lounesto, "Clifford algebras and spinors", Cambridge UP, 1997.

  template<const index_t LO, const index_t HI>
  inline
  const std::string
  index_set<LO,HI>::
  classname()
  { return "index_set"; }

  /// Constructor from index value
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set(const index_t idx)
  { this->set(idx); }

  /// Constructor from bitset_t
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set(const bitset_t bst):
  bitset_t(bst)
  { }

  /// Constructor from set value of an index set folded within the given frame
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set(const set_value_t folded_val, const index_set_t frm, const bool prechecked)
  {
    if (!prechecked && folded_val >= set_value_t(1 << frm.count()))
        throw error_t("index_set(val,frm): cannot create: value gives an index set outside of frame");
    const index_set_t folded_frame = frm.fold();
    const index_t min_index = folded_frame.min();
    const index_t skip = min_index > 0 ? 1 : 0;
    const index_set_t folded_set = index_set_t(bitset_t(folded_val) << (min_index - skip - LO));
    *this = folded_set.unfold(frm);
  }

  /// Constructor from range of indices from range.first to range.second
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set(const index_pair_t& range, const bool prechecked)
  {
    if (!prechecked && (range.first < LO || range.second > HI))
        throw error_t("index_set(range): cannot create: range is too large");
    const index_t begin_bit = (range.first < 0)
                            ? range.first-LO
                            : range.first-LO-1;
    const index_t end_bit = (range.second < 0)
                            ? range.second-LO+1
                            : range.second-LO;
    unsigned long mask = ( (end_bit == _GLUCAT_BITS_PER_ULONG)
                           ? -1UL
                           : (1UL << end_bit)-1UL)
                         & ~((1UL << begin_bit)-1UL);
    *this = bitset_t(mask);
  }

  /// Constructor from string
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set(const std::string& str)
  {
    std::istringstream ss(str);
    ss >> *this;
    if (!ss)
      throw error_t("index_set_t(str): could not parse string");
    // Peek to see if the end of the string has been reached.
    ss.peek();
    if (!ss.eof())
      throw error_t("index_set_t(str): could not parse entire string");
  }

  /// Equality
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  operator== (const index_set_t rhs) const
  {
    const auto* pthis = static_cast<const bitset_t*>(this);
    return *pthis == static_cast<bitset_t>(rhs);
  }

  /// Inequality
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  operator!= (const index_set_t rhs) const
  {
    const auto* pthis = static_cast<const bitset_t*>(this);
    return *pthis != static_cast<bitset_t>(rhs);
  }

  /// Set complement: not
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>
  index_set<LO,HI>::
  operator~ () const
  { return bitset_t::operator~(); }

  /// Symmetric set difference: exclusive or
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  operator^= (const index_set_t rhs)
  {
    bitset_t* pthis = this;
    *pthis ^= static_cast<bitset_t>(rhs);
    return *this;
  }

  /// Symmetric set difference: exclusive or
  template<const index_t LO, const index_t HI>
  inline
  const
  index_set<LO,HI>
  operator^ (const index_set<LO,HI>& lhs,
             const index_set<LO,HI>& rhs)
  {
    using index_set_t = index_set<LO, HI>;
    using bitset_t = typename index_set_t::bitset_t;
    return static_cast<bitset_t>(lhs) ^ static_cast<bitset_t>(rhs);
  }

  /// Set intersection: and
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  operator&= (const index_set_t rhs)
  {
    bitset_t* pthis = this;
    *pthis &= static_cast<bitset_t>(rhs);
    return *this;
  }

  /// Set intersection: and
  template<const index_t LO, const index_t HI>
  inline
  const
  index_set<LO,HI>
  operator& (const index_set<LO,HI>& lhs,
             const index_set<LO,HI>& rhs)
  {
    using index_set_t = index_set<LO, HI>;
    using bitset_t = typename index_set_t::bitset_t;
    return static_cast<bitset_t>(lhs) & static_cast<bitset_t>(rhs);
  }

  /// Set union: or
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  operator|= (const index_set_t rhs)
  {
    bitset_t* pthis = this;
    *pthis |= static_cast<bitset_t>(rhs);
    return *this;
  }

  /// Set union: or
  template<const index_t LO, const index_t HI>
  inline
  const
  index_set<LO,HI>
  operator| (const index_set<LO,HI>& lhs,
             const index_set<LO,HI>& rhs)
  {
    using index_set_t = index_set<LO, HI>;
    using bitset_t = typename index_set_t::bitset_t;
    return static_cast<bitset_t>(lhs) | static_cast<bitset_t>(rhs);
  }

  /// Subscripting: Element access
  template<const index_t LO, const index_t HI>
  inline
  typename index_set<LO,HI>::reference
  index_set<LO,HI>::
  operator[] (const index_t idx)
  { return reference(*this, idx); }

  /// Subscripting: Test idx for membership: test value of bit idx
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  operator[] (const index_t idx) const
  { return this->test(idx); }

  /// Test idx for membership: test value of bit idx
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  test(const index_t idx) const
  {
    // Reference: [JA], 1.2.1
    return (idx < 0)
           ?   bool(bitset_t::to_ulong() & (1UL << (idx - LO)))
           : (idx > 0)
             ? bool(bitset_t::to_ulong() & (1UL << (idx - LO - 1)))
             : false;
  }

  /// Include all indices except 0: set all bits except 0
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  set()
  {
    bitset_t::set();
    return *this;
  }

  /// Include idx: Set bit at idx if idx != 0
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  set(index_t idx)
  {
    if (idx > 0)
      bitset_t::set(idx-LO-1);
    else if (idx < 0)
      bitset_t::set(idx-LO);
    return *this;
  }

  /// Set membership of idx to val if idx != 0: Set bit at idx to val if idx != 0
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  set(const index_t idx, const int val)
  {
    if (idx > 0)
      bitset_t::set(idx-LO-1, val);
    else if (idx < 0)
      bitset_t::set(idx-LO, val);
    return *this;
  }

  /// Make set empty: Set all bits to 0
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  reset()
  {
    bitset_t::reset();
    return *this;
  }

  /// Exclude idx:  Set bit at idx to 0
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  reset(const index_t idx)
  {
    if (idx > 0)
      bitset_t::reset(idx-LO-1);
    else if (idx < 0)
      bitset_t::reset(idx-LO);
    return *this;
  }

  /// Set complement, except 0: flip all bits, except 0
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  flip()
  {
    bitset_t::flip();
    return *this;
  }

  /// Complement membership of idx if idx != 0: flip bit at idx if idx != 0
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  flip(const index_t idx)
  {
    if (idx > 0)
      bitset_t::flip(idx-LO-1);
    else if (idx < 0)
      bitset_t::flip(idx-LO);
    return *this;
  }

  /// Cardinality: Number of indices included in set
  template<const index_t LO, const index_t HI>
  inline
  index_t
  index_set<LO,HI>::
  count() const
  {
    unsigned long val = bitset_t::to_ulong();
    // Reference: [JA], 1.3
    if (val == 0)
      return 0;
    else
    {
      index_t result = 1;
      while (val &= val-1)
        ++result;
      return result;
    }
  }

  /// Number of negative indices included in set
  template<const index_t LO, const index_t HI>
  inline
  index_t
  index_set<LO,HI>::
  count_neg() const
  {
    static const index_set_t lo_mask = bitset_t((1UL << -LO) - 1UL);
    const index_set_t neg_part = *this & lo_mask;
    return neg_part.count();
  }

  /// Number of positive indices included in set
  template<const index_t LO, const index_t HI>
  inline
  index_t
  index_set<LO,HI>::
  count_pos() const
  {
    const auto* pthis = static_cast<const bitset_t*>(this);
    const index_set_t pos_part = *pthis >> -LO;
    return pos_part.count();
  }

#if (_GLUCAT_BITS_PER_ULONG == 64)
  /// Minimum member, or 0 if none: default for 64-bit cpus
  template<const index_t LO, const index_t HI>
  inline
  index_t
  index_set<LO,HI>::
  min() const
  {
    // Reference: [JA], 1.3
    unsigned long val = bitset_t::to_ulong();
    if (val == 0)
      return 0;
    else
    {
      val -= val & (val-1); // isolate lowest bit

      index_t idx = 0;
      const index_t nbits = HI - LO;

      if (nbits > 8)
      {
        if (val & 0xffffffff00000000ul)
          idx += 32;
        if (val & 0xffff0000ffff0000ul)
          idx += 16;
        if (val & 0xff00ff00ff00ff00ul)
          idx +=  8;
      }
        if (val & 0xf0f0f0f0f0f0f0f0ul)
          idx +=  4;
        if (val & 0xccccccccccccccccul)
          idx +=  2;
        if (val & 0xaaaaaaaaaaaaaaaaul)
          idx +=  1;

      return idx + ((idx < -LO) ? LO : LO+1);
    }
  }
#elif (_GLUCAT_BITS_PER_ULONG == 32)
  /// Minimum member, or 0 if none: default for 32-bit cpus
  template<const index_t LO, const index_t HI>
  inline
  index_t
  index_set<LO,HI>::
  min() const
  {
    // Reference: [JA], 1.3
    unsigned long val = bitset_t::to_ulong();
    if (val == 0)
      return 0;
    else
    {
      val -= val & (val-1); // isolate lowest bit

      index_t idx = 0;
      const index_t nbits = HI - LO;
      if (nbits > 8)
      {
        if (val & 0xffff0000ul)
          idx += 16;
        if (val & 0xff00ff00ul)
          idx +=  8;
      }
        if (val & 0xf0f0f0f0ul)
          idx +=  4;
        if (val & 0xccccccccul)
          idx +=  2;
        if (val & 0xaaaaaaaaul)
          idx +=  1;

      return idx + ((idx < -LO) ? LO : LO+1);
    }
  }
#else
  /// Minimum member, or 0 if none
  template<const index_t LO, const index_t HI>
  index_t
  index_set<LO,HI>::
  min() const
  {
    for (index_t
        idx = LO;
        idx != 0;
        ++idx)
      if (this->test(idx))
        return idx;
    for (index_t
        idx = 1;
        idx <= HI;
        ++idx)
      if (this->test(idx))
        return idx;
    return 0;
  }
#endif

#if (_GLUCAT_BITS_PER_ULONG == 64)
  /// Maximum member, or 0 if none: default for 64-bit cpus
  template<const index_t LO, const index_t HI>
  inline
  index_t
  index_set<LO,HI>::
  max() const
  {
    // Reference: [JA], 1.6
    unsigned long val = bitset_t::to_ulong();
    if (val == 0)
      return 0;
    else
    {
      index_t idx = 0;
      const index_t nbits = HI - LO;
      if (nbits > 8)
      {
        if (val & 0xffffffff00000000ul)
          { val >>= 32; idx += 32; }
        if (val & 0x00000000ffff0000ul)
          { val >>= 16; idx += 16; }
        if (val & 0x000000000000ff00ul)
          { val >>=  8; idx +=  8; }
      }
        if (val & 0x00000000000000f0ul)
          { val >>=  4; idx +=  4; }
        if (val & 0x000000000000000cul)
          { val >>=  2; idx +=  2; }
        if (val & 0x0000000000000002ul)
          {             idx +=  1; }
      return idx + ((idx < -LO) ? LO : LO+1);
    }
  }
#elif (_GLUCAT_BITS_PER_ULONG == 32)
  /// Maximum member, or 0 if none: default for 32-bit cpus
  template<const index_t LO, const index_t HI>
  inline
  index_t
  index_set<LO,HI>::
  max() const
  {
    // Reference: [JA], 1.6
    unsigned long val = bitset_t::to_ulong();
    if (val == 0)
      return 0;
    else
    {
      index_t idx = 0;
      const index_t nbits = HI - LO;
      if (nbits > 8)
      {
        if (val & 0xffff0000ul)
          { val >>= 16; idx += 16; }
        if (val & 0x0000ff00ul)
          { val >>=  8; idx +=  8; }
      }
        if (val & 0x000000f0ul)
          { val >>=  4; idx +=  4; }
        if (val & 0x0000000cul)
          { val >>=  2; idx +=  2; }
        if (val & 0x00000002ul)
          {             idx +=  1; }
      return idx + ((idx < -LO) ? LO : LO+1);
    }
  }
#else
  /// Maximum member, or 0 if none
  template<const index_t LO, const index_t HI>
  index_t
  index_set<LO,HI>::
  max() const
  {
    for (index_t
        idx = HI;
        idx != 0;
        --idx)
      if (this->test(idx))
        return idx;
    for (index_t
        idx = -1;
        idx >= LO;
        --idx)
      if (this->test(idx))
        return idx;
    return 0;
  }
#endif

  /// Lexicographic ordering of two sets: -1 if a<b, +1 if a>b, 0 if a==b
  //  eg. {3,4,5} is less than {3,7,8}
  template<const index_t LO, const index_t HI>
  inline
  int
  compare(const index_set<LO,HI>& a, const index_set<LO,HI>& b)
  {
    return (a == b)
           ? 0
           : a.lex_less_than(b)
             ? -1
             :  1;
  }

  /// Lexicographic ordering of two sets: *this < rhs
  //  eg. {3,4,5} is less than {3,7,8}
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  lex_less_than(const index_set_t rhs) const
  { return bitset_t::to_ulong() < rhs.bitset_t::to_ulong(); }

  /// Less than operator used for comparisons, map, etc.
  // Order by count, then order lexicographically within the equivalence class of count.
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  operator< (const index_set_t rhs) const
  {
    const index_t this_grade = this->count();
    const index_t rhs_grade  = rhs.count();
    return (this_grade < rhs_grade)
           ? true
           : (this_grade > rhs_grade)
             ? false
             : this->lex_less_than(rhs);
  }

  /// Write out index set
  template<const index_t LO, const index_t HI>
  std::ostream&
  operator<< (std::ostream& os, const index_set<LO,HI>& ist)
  {
    index_t i;
    os << '{';
    for (i= LO;
        (i <= HI) && !(ist[i]);
        ++i)
    { }
    if (i <= HI)
      os << i;
    for (++i;
        i <= HI;
        ++i)
      if (ist[i])
        os << ',' << i;
    os << '}';
    return os;
  }

  /// Read in index set
  template<const index_t LO, const index_t HI>
  std::istream&
  operator>> (std::istream& s, index_set<LO,HI>& ist)
  {
    // Parsing variables.
    int c = 0;
    index_t i = 0;
    index_set<LO,HI> local_ist;
    // Parsing control variables.
    bool parse_index_list = true;
    bool expect_closing_brace = false;
    bool expect_index = false;
    // Parse an optional opening brace.
    c = s.peek();
    // If there is a failure or end of file, this ends parsing.
    if (!s.good())
      parse_index_list = false;
    else
    { // Check for an opening brace.
      expect_closing_brace = (c == int('{'));
      if (expect_closing_brace)
      { // Consume the opening brace.
        s.get();
        // The next character may be a closing brace,
        // indicating the empty index set.
        c = s.peek();
        if (s.good() && (c == int('}')))
        { // A closing brace has been parsed and is no longer expected.
          expect_closing_brace = false;
          // Consume the closing brace.
          s.get();
          // This ends parsing.
          parse_index_list = false;
        }
      }
    }
    if (s.good() && parse_index_list)
    { // Parse an optional index list.
      // The index list starts with a first index.
      for (s >> i;
          !s.fail();
          s >> i)
      { // An index has been parsed. Check to see if it is in range.
        if ((i < LO) || (i > HI))
        { // An index out of range is a failure.
          s.clear(std::istream::failbit);
          break;
        }
        // Add the index to the index set local_ist.
        local_ist.set(i);
        // Immediately after parsing an index, an index is no longer expected.
        expect_index = false;
        // Reading the index may have resulted in an end of file condition.
        // If so, this ends the index list.
        if (s.eof())
          break;
        // The index list continues with a comma, and
        // may be ended by a closing brace, if it was begun with an opening brace.
        // Parse a possible comma or closing brace.
        c = s.peek();
        if (!s.good())
          break;
        // First, test for a closing brace, if expected.
        if (expect_closing_brace && (c == int('}')))
        { // Consume the closing brace.
          s.get();
          // Immediately after parsing the closing brace, it is no longer expected.
          expect_closing_brace = false;
          // A closing brace ends the index list.
          break;
        }
        // Now test for a comma.
        if (c == int(','))
        { // Consume the comma.
          s.get();
          // A index is expected after the comma.
          expect_index = true;
        }
        else
        { // Any other character here is a failure.
          s.clear(std::istream::failbit);
          break;
        }
      }
    }
    // If an index or a closing brace is still expected, this is a failure.
    if (expect_index || expect_closing_brace)
      s.clear(std::istream::failbit);
    // End of file is not a failure.
    if (s)
    { // The index set has been successfully parsed.
      ist = local_ist;
    }
    return s;
  }

  /// Determine if the index set is contiguous, ie. has no gaps when 0 is included
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  is_contiguous () const
  {
    const index_t min_index = this->min();
    const index_t max_index = this->max();
    return (min_index < 0 && max_index > 0)
         ?  max_index - min_index == this->count()
         : (min_index == 1 || max_index == -1) &&
           (max_index - min_index == this->count() - 1);
  }

  /// Fold this index set within itself as a frame
  template<const index_t LO, const index_t HI>
  inline
  const
  index_set<LO,HI>
  index_set<LO,HI>::
  fold() const
  { return this->fold(*this, true); }

  /// Fold this index set within the given frame
  template<const index_t LO, const index_t HI>
  const
  index_set<LO,HI>
  index_set<LO,HI>::
  fold(const index_set_t frm, const bool prechecked) const
  {
    if (!prechecked && ((*this | frm) != frm))
      throw error_t("fold(frm): cannot fold from outside of frame");
    const index_t frm_min = frm.min();
    const index_t frm_max = frm.max();
    index_set_t result;
    index_t fold_idx = -1;
    index_t unfold_idx;
    for (unfold_idx = -1;
        unfold_idx >= frm_min;
        --unfold_idx)
      if (frm.test(unfold_idx))
        // result.set(fold_idx--, this->test(unfold_idx));
      {
        if (this->test(unfold_idx))
          result.set(fold_idx);
        --fold_idx;
      }
    fold_idx = 1;
    for (unfold_idx = 1;
        unfold_idx <= frm_max;
        ++unfold_idx)
      if (frm.test(unfold_idx))
        // result.set(fold_idx++, this->test(unfold_idx));
      {
        if (this->test(unfold_idx))
          result.set(fold_idx);
        ++fold_idx;
      }
    return result;
  }

  /// Unfold this index set within the given frame
  template<const index_t LO, const index_t HI>
  const
  index_set<LO,HI>
  index_set<LO,HI>::
  unfold(const index_set_t frm, const bool prechecked) const
  {
    const char* msg =
      "unfold(frm): cannot unfold into a smaller frame";
    const index_t frm_min = frm.min();
    const index_t frm_max = frm.max();
    index_set_t result;
    index_t fold_idx = -1;
    index_t unfold_idx;
    for (unfold_idx = -1;
        unfold_idx >= frm_min;
        --unfold_idx)
      if (frm.test(unfold_idx))
        if (this->test(fold_idx--))
          result.set(unfold_idx);
    if (!prechecked && ((fold_idx+1) > this->min()))
      throw error_t(msg);
    fold_idx = 1;
    for (unfold_idx = 1;
        unfold_idx <= frm_max;
        ++unfold_idx)
      if (frm.test(unfold_idx))
        if (this->test(fold_idx++))
          result.set(unfold_idx);
    if (!prechecked && ((fold_idx-1) < this->max()))
      throw error_t(msg);
    return result;
  }

  /// The set value of the fold of this index set within the given frame
  template<const index_t LO, const index_t HI>
  inline
  set_value_t
  index_set<LO,HI>::
  value_of_fold(const index_set_t frm) const
  {
    const index_t min_index = frm.fold().min();
    if (min_index == 0)
      return 0;
    else
    {
      const index_set_t folded_set = this->fold(frm);
      const index_t skip = min_index > 0 ? 1 : 0;
      return folded_set.bitset_t::to_ulong() >> (min_index-LO-skip);
    }
  }

  /// Inverse reversed Gray code
  inline
  static
  unsigned long inverse_reversed_gray(unsigned long x)
  {
    // Reference: [JA]
#if (_GLUCAT_BITS_PER_ULONG >= 64)
    x ^= x << 32; // for 64-bit words
#endif
    x ^= x << 16; // reversed_gray ** 16
    x ^= x <<  8; // reversed_gray **  8
    x ^= x <<  4; // reversed_gray **  4
    x ^= x <<  2; // reversed_gray **  2
    x ^= x <<  1; // reversed_gray **  1
    return x;
  }

  /// Inverse Gray code
  inline
  static
  unsigned long inverse_gray(unsigned long x)
  {
    // Reference: [JA]
#if (_GLUCAT_BITS_PER_ULONG >= 64)
    x ^= x >> 32; // for 64-bit words
#endif
    x ^= x >> 16; // gray ** 16
    x ^= x >>  8; // gray **  8
    x ^= x >>  4; // gray **  4
    x ^= x >>  2; // gray **  2
    x ^= x >>  1; // gray **  1
    return x;
  }

  /// Sign of geometric product of two Clifford basis elements
  template<const index_t LO, const index_t HI>
  int
  index_set<LO,HI>::
  sign_of_mult(const index_set_t rhs) const
  {
    // Implemented using Walsh functions and Gray codes.
    // Reference: [L] Chapter 21, 21.3
    // Reference: [JA]
    const unsigned long uthis = this->bitset_t::to_ulong();
    const unsigned long urhs  =   rhs.bitset_t::to_ulong();
    const index_t nbits = HI - LO;
    unsigned long negative = 0;
    if (nbits > 8)
    {
      // Set h to be the inverse reversed Gray code of rhs.
      // This sets each bit of h to be the cumulative ^ of
      // the same and lower bits of rhs.
      const unsigned long h = inverse_reversed_gray(urhs);
      // Set k to be the inverse Gray code of *this & h.
      // This sets the low bit of k to be parity(*this & h).
      const unsigned long k = inverse_gray(uthis & h);
      // Set q to be the inverse Gray code of the positive part of *this & rhs.
      const unsigned long q = inverse_gray((uthis & urhs) >> -LO);
      negative = k ^ q;
    }
    else
    {
      unsigned long h = 0;
      index_t j;
      for (j = 0;
          j < -LO;
          ++j)
      {
        h ^= urhs >> j;
        negative ^= h & (uthis >> j);
      }
      for (j = -LO;
          j < nbits;
          ++j)
      {
        negative ^= h & (uthis >> j);
        h ^= urhs >> j;
      }
    }
    return 1 - int((negative & 1) << 1);
  }

  /// Sign of geometric square of a Clifford basis element
  template<const index_t LO, const index_t HI>
  inline
  int
  index_set<LO,HI>::
  sign_of_square() const
  {
    int result = 1 - int((this->count_neg() % 2) << 1);
    switch (this->count() % 4)
    {
      case 2:
      case 3:
        result *= -1;
        break;
      default:
        break;
    }
    return result;
  }

  /// Hash function
  template<const index_t LO, const index_t HI>
  inline
  size_t
  index_set<LO,HI>::
  hash_fn() const
  {
    static const unsigned long lo_mask = (1UL << -LO) - 1UL;
    const unsigned long uthis = bitset_t::to_ulong();
    const unsigned long neg_part = uthis & lo_mask;
    const unsigned long pos_part = uthis >> -LO;
    return size_t(neg_part ^ pos_part);
  }

  /// Square of generator index j
  inline
  int
  sign_of_square(index_t j)
  { return (j < 0) ? -1 : 1; }

  /// Minimum negative index, or 0 if none
  template<const index_t LO, const index_t HI>
  inline
  index_t
  min_neg(const index_set<LO,HI>& ist)
  { return std::min(ist.min(), 0); }

  /// Maximum positive index, or 0 if none
  template<const index_t LO, const index_t HI>
  inline
  index_t
  max_pos(const index_set<LO,HI>& ist)
  { return std::max(ist.max(), 0); }

// index_set reference

  /// index_set reference
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>::reference::
  reference( index_set_t& ist, index_t idx ) :
    m_pst(&ist),
    m_idx(idx)
  { }

  /// for b[i] = x;
  template<const index_t LO, const index_t HI>
  inline
  typename index_set<LO,HI>::reference&
  index_set<LO,HI>::reference::
  operator= (bool x)
  {
    if ( x )
      m_pst->set(m_idx);
    else
      m_pst->reset(m_idx);
    return *this;
  }

  /// for b[i] = b[j];
  template<const index_t LO, const index_t HI>
  inline
  typename index_set<LO,HI>::reference&
  index_set<LO,HI>::reference::
  operator= (const reference& j)
  {
    if ( (j.m_pst)[j.m_idx] )
      m_pst->set(m_idx);
    else
      m_pst->reset(m_idx);
    return *this;
  }

  /// flips the bit
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::reference::
  operator~ () const
  { return !(m_pst->test(m_idx)); }

  /// for x = b[i];
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>::reference::
  operator bool () const
  { return m_pst->test(m_idx); }

  /// for b[i].flip();
  template<const index_t LO, const index_t HI>
  inline
  typename index_set<LO,HI>::reference&
  index_set<LO,HI>::reference::
  flip()
  {
    m_pst->flip(m_idx);
    return *this;
  }
}
#endif // _GLUCAT_INDEX_SET_IMP_H
