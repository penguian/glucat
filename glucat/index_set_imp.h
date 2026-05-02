#ifndef _GLUCAT_INDEX_SET_IMP_H
#define _GLUCAT_INDEX_SET_IMP_H
/**************************************************************************
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

#include <string>
#include <sstream>

namespace glucat
{
  // References for algorithms:
  // [JA]: Joerg Arndt,  "Algorithms for programmers", http://www.jjj.de/fxt/fxtbook.pdf
  //      Chapter 1, Bit wizardry, http://www.jjj.de/bitwizardry/bitwizardrypage.html
  // [L]: Pertti Lounesto, "Clifford algebras and spinors", Cambridge UP, 1997.

  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  classname() -> std::string_view
  { return "index_set"; }

  /*
   * @brief Constructor from index value
   * @details
   * @tparam LO
   * @tparam HI
   * @param idx Value
   */
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set(const index_t idx)
  { this->set(idx); }

  /*
   * @brief Constructor from bitset_t
   * @details
   * @tparam LO
   * @tparam HI
   * @param bst Value
   */
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set(const bitset_t bst):
  bitset_t(bst)
  { }

  /*
   * @brief Constructor from set value of an index set folded within the given frame
   * @details
   * @tparam LO
   * @tparam HI
   * @param folded_val Value
   * @param frm Value
   * @param prechecked Already checked?
   */
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set(const set_value_t folded_val, const index_set_t frm, const bool prechecked)
  {
    if (!prechecked && folded_val >= (set_value_t(1) << frm.count()))
        throw error_t("index_set(val,frm): cannot create: value gives an index set outside of frame");
    const index_set_t folded_frame = frm.fold();
    const index_t min_index = folded_frame.min();
    const index_t skip = min_index > 0 ? 1 : 0;
    const index_set_t folded_set = index_set_t(bitset_t(folded_val) << (min_index - skip - LO));
    *this = folded_set.unfold(frm);
  }

  /*
   * @brief Constructor from range of indices from range.first to range.second
   * @details
   * @tparam LO
   * @tparam HI
   * @param range Value
   * @param prechecked Already checked?
   */
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

  /*
   * @brief Constructor from string
   * @details
   * @tparam LO
   * @tparam HI
   * @param str Value
   */
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

  /*
   * @brief Equality
   * @details
   * @tparam LO
   * @tparam HI
   * @param rhs Right hand side
   * @return True if equal
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  operator== (const index_set_t& rhs) const -> bool
  {
    const auto* pthis = static_cast<const bitset_t*>(this);
    return *pthis == static_cast<const bitset_t&>(rhs);
  }

  /*
   * @brief Inequality
   * @details
   * @tparam LO
   * @tparam HI
   * @param rhs Right hand side
   * @return True if not equal
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  operator!= (const index_set_t& rhs) const -> bool
  {
    const auto* pthis = static_cast<const bitset_t*>(this);
    return *pthis != static_cast<const bitset_t&>(rhs);
  }

  /*
   * @brief Set complement: not
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  operator~ () const -> index_set_t
  { return bitset_t::operator~(); }

  /*
   * @brief Symmetric set difference: exclusive or
   * @details
   * @tparam LO
   * @tparam HI
   * @param rhs Right hand side
   * @return Reference to this
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  operator^= (const index_set_t rhs) -> index_set_t&
  {
    bitset_t* pthis = this;
    *pthis ^= static_cast<bitset_t>(rhs);
    return *this;
  }

  /*
   * @brief Symmetric set difference: exclusive or
   * @details
   *
   * Usage example:
   * Location: glucat/framed_multi_imp.h:2433
   *
   * @code
   *
   * return term_t(lhs.first ^ rhs.first, crd_of_mult(lhs, rhs));
   * @endcode
   *
   * @par Example:
   * @code
   * index_set<>("{1}") ^ index_set<>("{2}"); // Returns {1,2}
   * index_set<>("{1,2}") ^ index_set<>("{2}"); // Returns {1}
   * @endcode
   *
   * @tparam LO
   * @tparam HI
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Outer product
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  operator^ (const index_set<LO,HI>& lhs,
             const index_set<LO,HI>& rhs) -> index_set<LO,HI>
  {
    using index_set_t = index_set<LO, HI>;
    using bitset_t = typename index_set_t::bitset_t;
    return static_cast<bitset_t>(lhs) ^ static_cast<bitset_t>(rhs);
  }

  /*
   * @brief Set intersection: and
   * @details
   * @tparam LO
   * @tparam HI
   * @param rhs Right hand side
   * @return Reference to this
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  operator&= (const index_set_t rhs) -> index_set_t&
  {
    bitset_t* pthis = this;
    *pthis &= static_cast<bitset_t>(rhs);
    return *this;
  }

  /*
   * @brief Set intersection: and
   * @details
   * @par Example:
   * @code
   * index_set<>("{1}") & index_set<>("{2}"); // Returns {}
   * index_set<>("{1,2}") & index_set<>("{2}"); // Returns {2}
   * @endcode
   *
   * @tparam LO
   * @tparam HI
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Inner product
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  operator& (const index_set<LO,HI>& lhs,
             const index_set<LO,HI>& rhs) -> index_set<LO,HI>
  {
    using index_set_t = index_set<LO, HI>;
    using bitset_t = typename index_set_t::bitset_t;
    return static_cast<bitset_t>(lhs) & static_cast<bitset_t>(rhs);
  }

  /*
   * @brief Set union: or
   * @details
   * @tparam LO
   * @tparam HI
   * @param rhs Right hand side
   * @return Reference to this
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  operator|= (const index_set_t rhs) -> index_set_t&
  {
    bitset_t* pthis = this;
    *pthis |= static_cast<bitset_t>(rhs);
    return *this;
  }

  /*
   * @brief Set union: or
   * @details
   * @par Example:
   * @code
   * index_set<>("{1}") | index_set<>("{2}"); // Returns {1,2}
   * index_set<>("{1,2}") | index_set<>("{2}"); // Returns {1,2}
   * @endcode
   *
   * @tparam LO
   * @tparam HI
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Bitwise OR
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  operator| (const index_set<LO,HI>& lhs,
             const index_set<LO,HI>& rhs) -> index_set<LO,HI>
  {
    using index_set_t = index_set<LO, HI>;
    using bitset_t = typename index_set_t::bitset_t;
    return static_cast<bitset_t>(lhs) | static_cast<bitset_t>(rhs);
  }

  /*
   * @brief Subscripting: Element access
   * @details
   * @tparam LO
   * @tparam HI
   * @param idx Value
   * @return Element reference
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  operator[] (const index_t idx) -> reference
  { return reference(*this, idx); }

  /*
   * @brief Subscripting: Test idx for membership: test value of bit idx
   * @details
   * @tparam LO
   * @tparam HI
   * @param idx Value
   * @return Element reference
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  operator[] (const index_t idx) const -> bool
  { return this->test(idx); }

  /*
   * @brief Test idx for membership: test value of bit idx
   * @details
   * @tparam LO
   * @tparam HI
   * @param idx Value
   * @return True if successful or condition met
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  test(const index_t idx) const -> bool
  {
    // Reference: [JA], 1.2.1
    return (idx < 0)
           ?   bool(bitset_t::to_ulong() & (1UL << (idx - LO)))
           : (idx > 0)
             ? bool(bitset_t::to_ulong() & (1UL << (idx - LO - 1)))
             : false;
  }

  /*
   * @brief Include all indices except 0: set all bits except 0
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  set() -> index_set_t&
  {
    bitset_t::set();
    return *this;
  }

  /*
   * @brief Include idx: Set bit at idx if idx != 0
   * @details
   * @tparam LO
   * @tparam HI
   * @param idx Value
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  set(index_t idx) -> index_set_t&
  {
    if (idx > 0)
      bitset_t::set(idx-LO-1);
    else if (idx < 0)
      bitset_t::set(idx-LO);
    return *this;
  }

  /*
   * @brief Set membership of idx to val if idx != 0: Set bit at idx to val if idx != 0
   * @details
   * @tparam LO
   * @tparam HI
   * @param idx Value
   * @param val Value
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  set(const index_t idx, const int val) -> index_set_t&
  {
    if (idx > 0)
      bitset_t::set(idx-LO-1, val);
    else if (idx < 0)
      bitset_t::set(idx-LO, val);
    return *this;
  }

  /*
   * @brief Make set empty: Set all bits to 0
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  reset() -> index_set_t&
  {
    bitset_t::reset();
    return *this;
  }

  /*
   * @brief Exclude idx:  Set bit at idx to 0
   * @details
   * @tparam LO
   * @tparam HI
   * @param idx Value
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  reset(const index_t idx) -> index_set_t&
  {
    if (idx > 0)
      bitset_t::reset(idx-LO-1);
    else if (idx < 0)
      bitset_t::reset(idx-LO);
    return *this;
  }

  /*
   * @brief Set complement, except 0: flip all bits, except 0
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  flip() -> index_set<LO,HI>&
  {
    bitset_t::flip();
    return *this;
  }

  /*
   * @brief Complement membership of idx if idx != 0: flip bit at idx if idx != 0
   * @details
   * @tparam LO
   * @tparam HI
   * @param idx Value
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  flip(const index_t idx) -> index_set_t&
  {
    if (idx > 0)
      bitset_t::flip(idx-LO-1);
    else if (idx < 0)
      bitset_t::flip(idx-LO);
    return *this;
  }

  /*
   * @brief Cardinality: Number of indices included in set
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  count() const -> index_t
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

  /*
   * @brief Number of negative indices included in set
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  count_neg() const -> index_t
  {
    static const index_set_t lo_mask = bitset_t((1UL << -LO) - 1UL);
    const index_set_t neg_part = *this & lo_mask;
    return neg_part.count();
  }

  /*
   * @brief Number of positive indices included in set
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  count_pos() const -> index_t
  {
    const auto* pthis = static_cast<const bitset_t*>(this);
    const index_set_t pos_part = *pthis >> -LO;
    return pos_part.count();
  }

#if (_GLUCAT_BITS_PER_ULONG == 64)
  /*
   * @brief Minimum member, or 0 if none: default for 64-bit cpus
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  min() const -> index_t
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
  /*
   * @brief Minimum member, or 0 if none: default for 32-bit cpus
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
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
  /*
   * @brief Minimum member, or 0 if none
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  auto
  index_set<LO,HI>::
  min() const -> index_t
  {
    for (auto
        idx = LO;
        idx != 0;
        ++idx)
      if (this->test(idx))
        return idx;
    for (auto
        idx = index_t(1);
        idx <= HI;
        ++idx)
      if (this->test(idx))
        return idx;
    return 0;
  }
#endif

#if (_GLUCAT_BITS_PER_ULONG == 64)
  /*
   * @brief Maximum member, or 0 if none: default for 64-bit cpus
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  max() const -> index_t
  {
    // Reference: [JA], 1.6
    auto val = bitset_t::to_ulong();
    if (val == 0)
      return 0;
    else
    {
      auto idx = index_t(0);
      const auto nbits = HI - LO;
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
  /*
   * @brief Maximum member, or 0 if none: default for 32-bit cpus
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  max() const -> index_t
  {
    // Reference: [JA], 1.6
    auto val = bitset_t::to_ulong();
    if (val == 0)
      return 0;
    else
    {
      auto idx = index_t(0);
      const auto nbits = HI - LO;
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
  /*
   * @brief Maximum member, or 0 if none
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  auto
  index_set<LO,HI>::
  max() const -> index_t
  {
    for (auto
        idx = HI;
        idx != 0;
        --idx)
      if (this->test(idx))
        return idx;
    for (auto
        idx = index_t(-1);
        idx >= LO;
        --idx)
      if (this->test(idx))
        return idx;
    return 0;
  }
#endif

  /*
   * @brief Lexicographic ordering of two sets: -1 if a<b, +1 if a>b, 0 if a==b
   * @details
   * @param a Value
   * @param b Value
   * @return Result
   */
  //  eg. {3,4,5} is less than {3,7,8}
  template<const index_t LO, const index_t HI>
  inline
  auto
  compare(const index_set<LO,HI>& a, const index_set<LO,HI>& b) -> int
  {
    return (a == b)
           ? 0
           : a.lex_less_than(b)
             ? -1
             :  1;
  }

  /*
   * @brief Lexicographic ordering of two sets: *this < rhs
   * @details
   * @param rhs Right hand side
   * @return Result
   */
  //  eg. {3,4,5} is less than {3,7,8}
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  lex_less_than(const index_set_t& rhs) const -> bool
  { return bitset_t::to_ulong() < rhs.bitset_t::to_ulong(); }

  /*
   * @brief Less than operator used for comparisons, map, etc.
   * @details
   * @param rhs Right hand side
   */
  // Order by count, then order lexicographically within the equivalence class of count.
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  operator< (const index_set_t& rhs) const -> bool
  {
    const auto this_grade = this->count();
    const auto rhs_grade  = rhs.count();
    return (this_grade < rhs_grade)
           ? true
           : (this_grade > rhs_grade)
             ? false
             : this->lex_less_than(rhs);
  }

  /*
   * @brief Write out index set
   * @details
   * @tparam LO
   * @tparam HI
   * @param os Output stream
   * @param ist Value
   * @return Output stream
   */
  template<const index_t LO, const index_t HI>
  auto
  operator<< (std::ostream& os, const index_set<LO,HI>& ist) -> std::ostream&
  {
    index_t i;
    os << '{';
    for (i = LO;
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

  /*
   * @brief Read in index set
   * @details
   * @tparam LO
   * @tparam HI
   * @param s Value
   * @param ist Value
   * @return Input stream
   */
  template<const index_t LO, const index_t HI>
  auto
  operator>> (std::istream& s, index_set<LO,HI>& ist) -> std::istream&
  {
    // Parsing variables.
    auto i = index_t(0);
    using index_set_t = index_set<LO,HI>;
    auto local_ist = index_set_t();
    // Parsing control variables.
    auto parse_index_list = true;
    auto expect_closing_brace = false;
    auto expect_index = false;
    // Parse an optional opening brace.
    auto c = s.peek();
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

  /*
   * @brief Determine if the index set is contiguous, ie. has no gaps when 0 is included
   * @details
   * @tparam LO
   * @tparam HI
   * @return True if is contiguous
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  is_contiguous () const -> bool
  {
    const auto min_index = this->min();
    const auto max_index = this->max();
    return (min_index < 0 && max_index > 0)
         ?  max_index - min_index == this->count()
         : (min_index == 1 || max_index == -1) &&
           (max_index - min_index == this->count() - 1);
  }

  /*
   * @brief Fold this index set within itself as a frame
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  fold() const -> index_set<LO,HI>
  { return this->fold(*this, true); }

  /*
   * @brief Fold this index set within the given frame
   * @details
   * @tparam LO
   * @tparam HI
   * @param frm Value
   * @param prechecked Already checked?
   * @return True if successful or condition met
   */
  template<const index_t LO, const index_t HI>
  auto
  index_set<LO,HI>::
  fold(const index_set_t frm, const bool prechecked) const -> index_set<LO,HI>
  {
    if (!prechecked && ((*this | frm) != frm))
      throw error_t("fold(frm): cannot fold from outside of frame");
    const auto frm_min = frm.min();
    const auto frm_max = frm.max();
    auto result = index_set_t();
    auto fold_idx = index_t(-1);
    for (auto
        unfold_idx = fold_idx;
        unfold_idx >= frm_min;
        --unfold_idx)
      if (frm.test(unfold_idx))
        // result.set(fold_idx--, this->test(unfold_idx));
      {
        if (this->test(unfold_idx))
          result.set(fold_idx);
        --fold_idx;
      }
    fold_idx = index_t(1);
    for (auto
        unfold_idx = fold_idx;
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

  /*
   * @brief Unfold this index set within the given frame
   * @details
   * @tparam LO
   * @tparam HI
   * @param frm Value
   * @param prechecked Already checked?
   * @return True if successful or condition met
   */
  template<const index_t LO, const index_t HI>
  auto
  index_set<LO,HI>::
  unfold(const index_set_t frm, const bool prechecked) const -> index_set_t
  {
    const char* msg =
      "unfold(frm): cannot unfold into a smaller frame";
    const auto frm_min = frm.min();
    const auto frm_max = frm.max();
    auto result = index_set_t();
    auto fold_idx = index_t(-1);
    for (auto
        unfold_idx = fold_idx;
        unfold_idx >= frm_min;
        --unfold_idx)
      if (frm.test(unfold_idx))
        if (this->test(fold_idx--))
          result.set(unfold_idx);
    if (!prechecked && ((fold_idx+1) > this->min()))
      throw error_t(msg);
    fold_idx = index_t(1);
    for (auto
        unfold_idx = fold_idx;
        unfold_idx <= frm_max;
        ++unfold_idx)
      if (frm.test(unfold_idx))
        if (this->test(fold_idx++))
          result.set(unfold_idx);
    if (!prechecked && ((fold_idx-1) < this->max()))
      throw error_t(msg);
    return result;
  }

  /*
   * @brief The set value of the fold of this index set within the given frame
   * @details
   * @tparam LO
   * @tparam HI
   * @param frm Value
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  value_of_fold(const index_set_t frm) const -> set_value_t
  {
    const auto min_index = frm.fold().min();
    if (min_index == 0)
      return 0;
    else
    {
      const auto folded_set = this->fold(frm);
      const auto skip = min_index > 0 ? index_t(1) : index_t(0);
      return folded_set.bitset_t::to_ulong() >> (min_index-LO-skip);
    }
  }

  /*
   * @brief Inverse reversed Gray code
   * @details
   * @param x Value
   * @return Inverse
   */
  inline
  static
  auto inverse_reversed_gray(unsigned long x) -> unsigned long
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

  /*
   * @brief Inverse Gray code
   * @details
   * @param x Value
   * @return Inverse
   */
  inline
  static
  auto inverse_gray(unsigned long x) -> unsigned long
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

  /*
   * @brief Sign of geometric product of two Clifford basis elements
   * @details
   *
   * Usage example:
   * Location: glucat/framed_multi_imp.h:2414
   *
   * @code
   *
   * { return lhs.first.sign_of_mult(rhs.first) * lhs.second * rhs.second; }
   * @endcode
   *
   * @tparam LO
   * @tparam HI
   * @param rhs Right hand side
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  auto
  index_set<LO,HI>::
  sign_of_mult(const index_set_t& rhs) const -> int
  {
    // Implemented using Walsh functions and Gray codes.
    // Reference: [L] Chapter 21, 21.3
    // Reference: [JA]
    const auto uthis = this->bitset_t::to_ulong();
    const auto urhs  =   rhs.bitset_t::to_ulong();
    const auto nbits = HI - LO;
    auto negative = 0UL;
    if (nbits > 8)
    {
      // Set h to be the inverse reversed Gray code of rhs.
      // This sets each bit of h to be the cumulative ^ of
      // the same and lower bits of rhs.
      const auto h = inverse_reversed_gray(urhs);
      // Set k to be the inverse Gray code of *this & h.
      // This sets the low bit of k to be parity(*this & h).
      const auto k = inverse_gray(uthis & h);
      // Set q to be the inverse Gray code of the positive part of *this & rhs.
      const auto q = inverse_gray((uthis & urhs) >> -LO);
      negative = k ^ q;
    }
    else
    {
      auto h = 0UL;
      for (auto
          j = index_t(0);
          j < -LO;
          ++j)
      {
        h ^= urhs >> j;
        negative ^= h & (uthis >> j);
      }
      for (auto
          j = index_t(-LO);
          j < nbits;
          ++j)
      {
        negative ^= h & (uthis >> j);
        h ^= urhs >> j;
      }
    }
    return 1 - int((negative & 1) << 1);
  }

  /*
   * @brief Sign of geometric square of a Clifford basis element
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  sign_of_square() const -> int
  {
    auto result = 1 - int((this->count_neg() % 2) << 1);
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

  /*
   * @brief Hash function
   * @details
   * @tparam LO
   * @tparam HI
   * @return Size
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::
  hash_fn() const -> size_t
  {
    static const auto lo_mask = (1UL << -LO) - 1UL;
    const auto uthis = bitset_t::to_ulong();
    const auto neg_part = uthis & lo_mask;
    const auto pos_part = uthis >> -LO;
    return size_t(neg_part ^ pos_part);
  }

  /*
   * @brief Square of generator index j
   * @details
   * @param j Column index
   * @return Result
   */
  inline
  auto
  sign_of_square(index_t j) -> int
  { return (j < 0) ? -1 : 1; }

  /*
   * @brief Minimum negative index, or 0 if none
   * @details
   * @tparam LO
   * @tparam HI
   * @param ist Value
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  min_neg(const index_set<LO,HI>& ist) -> index_t
  { return std::min(ist.min(), 0); }

  /*
   * @brief Maximum positive index, or 0 if none
   * @details
   * @tparam LO
   * @tparam HI
   * @param ist Value
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  max_pos(const index_set<LO,HI>& ist) -> index_t
  { return std::max(ist.max(), 0); }

// index_set reference

  /*
   * @brief index_set reference
   * @details
   * @tparam LO
   * @tparam HI
   * @param ist Value
   * @param idx Value
   */
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>::reference::
  reference( index_set_t& ist, index_t idx ) :
    m_pst(&ist),
    m_idx(idx)
  { }

  /*
   * @brief for b[i] == c[j];
   * @details
   * @tparam LO
   * @tparam HI
   * @param c_j Value
   * @return True if equal
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::reference::
  operator== (const reference& c_j) const -> bool
  { return m_pst == c_j.m_pst && m_idx == c_j.m_idx; }

  /*
   * @brief for b[i] = x;
   * @details
   * @tparam LO
   * @tparam HI
   * @param x Value
   * @return Reference to this
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::reference::
  operator= (bool x) -> reference&
  {
    if ( x )
      m_pst->set(m_idx);
    else
      m_pst->reset(m_idx);
    return *this;
  }

  /*
   * @brief for b[i] = c[j];
   * @details
   * @tparam LO
   * @tparam HI
   * @param c_j Value
   * @return Reference to this
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::reference::
  operator= (const reference& c_j) -> reference&
  {
    if (&c_j != this && c_j != *this)
    {
      if ( (*c_j.m_pst)[c_j.m_idx] )
        m_pst->set(m_idx);
      else
        m_pst->reset(m_idx);
    }
    return *this;
  }

  /*
   * @brief flips the bit
   * @details
   * @tparam LO
   * @tparam HI
   * @return True if successful or condition met
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::reference::
  operator~ () const -> bool
  { return !(m_pst->test(m_idx)); }

  /*
   * @brief for x = b[i];
   * @details
   * @tparam LO
   * @tparam HI
   * @return True if successful or condition met
   */
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>::reference::
  operator bool () const
  { return m_pst->test(m_idx); }

  /*
   * @brief for b[i].flip();
   * @details
   * @tparam LO
   * @tparam HI
   * @return Result
   */
  template<const index_t LO, const index_t HI>
  inline
  auto
  index_set<LO,HI>::reference::
  flip() -> reference&
  {
    m_pst->flip(m_idx);
    return *this;
  }
}
#ifdef GLUCAT_DOCTEST
#include <iostream>
#include <sstream>

TEST_CASE("index_set<LO,HI>") {
  using is_t = glucat::index_set<-32, 32>;

  SUBCASE("Metadata") {
    CHECK(is_t::classname() == "index_set");
  }

  SUBCASE("Constructor and string representation") {
    is_t s1(1);
    std::ostringstream oss1;
    oss1 << s1;
    CHECK(oss1.str() == "{1}");

    is_t s2("{1,2}");
    std::ostringstream oss2;
    oss2 << s2;
    CHECK(oss2.str() == "{1,2}");

    is_t s3("");
    std::ostringstream oss3;
    oss3 << s3;
    CHECK(oss3.str() == "{}");
  }

  SUBCASE("Comparisons") {
    CHECK(is_t(1) == is_t("{1}"));
    CHECK(is_t("{1}") != is_t("{2}"));
    CHECK(is_t("{1}") < is_t("{2}"));
    CHECK(is_t("{1}") <= is_t("{2}"));
    CHECK_FALSE(is_t("{1}") > is_t("{2}"));
    CHECK_FALSE(is_t("{1}") >= is_t("{2}"));
  }

  SUBCASE("Set operations") {
    is_t s1("{1}");
    is_t s2("{2}");
    CHECK((s1 ^ s2) == is_t("{1,2}"));
    CHECK((is_t("{1,2}") ^ s2) == s1);
    CHECK((is_t("{1,2}") & s2) == s2);
    CHECK((is_t("{1}") & s2) == is_t("{}"));
    CHECK((s1 | s2) == is_t("{1,2}"));
  }

  SUBCASE("Cardinality and bounds") {
    is_t s("{-1,1,2}");
    CHECK(s.count() == 3);
    CHECK(s.count_neg() == 1);
    CHECK(s.count_pos() == 2);
    CHECK(s.min() == -1);
    CHECK(s.max() == 2);
  }

  SUBCASE("Sign functions") {
    is_t s1("{1,2}");
    is_t s2("{-1}");
    CHECK(s1.sign_of_mult(s2) == 1);
    CHECK(s1.sign_of_square() == -1);
  }

  SUBCASE("Adversarial and edge cases") {
    // sign_of_square for single index
    CHECK(glucat::sign_of_square(1) == 1);
    CHECK(glucat::sign_of_square(-1) == -1);

    // reference operators
    is_t s;
    s[1] = true;
    CHECK(s.test(1));
    s[1] = false;
    CHECK_FALSE(s.test(1));
    s[1].flip();
    CHECK(s.test(1));

    is_t s2;
    s2[2] = s[1];
    CHECK(s2.test(2));
    s2[2] = s[3]; // s[3] is false
    CHECK_FALSE(s2.test(2));

    // operator~ and operator bool (via reference)
    is_t s3("{1}");
    CHECK(s3[1]);
    CHECK_FALSE(s3[2]);
    CHECK((~s3).test(1) == false);
    CHECK((~s3).test(2) == true);

    // reference operators
    is_t s_ref;
    s_ref[1] = true;
    CHECK(s_ref[1]);
    CHECK_FALSE(~s_ref[1]);
    CHECK(s_ref[1] == s_ref[1]); // same index set and index
    s_ref[1] = s_ref[1]; // self assignment

    // Comparisons with more varied sets
    CHECK(is_t("{1,2}") < is_t("{1,2,3}"));
    CHECK(is_t("{1,3}") > is_t("{1,2}"));
    CHECK(is_t("{-1,1}") != is_t("{1}"));
  }

  SUBCASE("Exceptions") {
    CHECK_THROWS(is_t("{invalid}"));
    CHECK_THROWS(is_t("{1,invalid}"));
    CHECK_THROWS(is_t("{1,,2}"));
    // Out of frame: is_t is <-32, 32>, so index 33 should throw if we try to set it via string or other means
    // that check bounds.
    // The constructor index_set(index_t val) uses prechecked=false by default, which should throw if out of bounds.
    CHECK_THROWS(is_t(33));
    CHECK_THROWS(is_t(-33));
    CHECK_THROWS(is_t("{1} garbage"));
  }

  SUBCASE("Static Factory and Bit-Wizardry") {
    // Static factory methods (verified at compile-time, runtime check for consistency)
    CHECK(is_t::from_index<1>() == is_t(1));
    CHECK(is_t::from_range<-1, 1>() == is_t(typename is_t::index_pair_t(-1, 1)));

    // Bit-wizardry coverage for 64-bit paths (requires LO <= -32, HI >= 32)
    using large_is_t = glucat::index_set<-32, 32>;
    large_is_t s_large;
    s_large.set(32);
    CHECK(s_large.max() == 32);
    CHECK(s_large.min() == 32);
    s_large.set(-32);
    CHECK(s_large.min() == -32);
    CHECK(s_large.max() == 32);

    large_is_t s_mid;
    s_mid.set(16);
    CHECK(s_mid.max() == 16);
    s_mid.set(8);
    CHECK(s_mid.min() == 8);
  }
}
#endif

#endif // _GLUCAT_INDEX_SET_IMP_H
