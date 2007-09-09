#ifndef _GLUCAT_INDEX_SET_IMP_H
#define _GLUCAT_INDEX_SET_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    index_set_imp.h : Implement a class for a set of non-zero integer indices
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2007 by Paul C. Leopardi
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
  // [JA]: Joerg Arndt,  "Algorithms for programmers", http://www.jjj.de/fxt/fxtbook.pdf
  //      Chapter 1, Bit wizardry, http://www.jjj.de/bitwizardry/bitwizardrypage.html
  // [L]: Pertti Lounesto, "Clifford algebras and spinors", Cambridge UP, 1997.

  template<const index_t LO, const index_t HI>
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
    if (!prechecked)
      if (folded_val > set_value_t(1 << frm.count()))
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
  index_set(const std::pair<index_t,index_t>& range, const bool prechecked)
  {
    if (!prechecked)
      if ((range.first < LO) || (range.second > HI))
        throw error_t("index_set(range): cannot create: range is too large");
    for (index_t
        idx =  range.first;
        idx <= range.second;
        ++idx)
      if (idx != 0)
        this->set(idx);
  }

  /// Constructor from string
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set(const std::string& str)
  {
    std::istringstream ss(str);
    ss >> *this;
  }

  /// Equality
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  operator== (const index_set_t rhs) const
  {
    const bitset_t* pthis = this;
    const bitset_t* pthat = &rhs;
    return *pthis == *pthat;
  }

  /// Inequality
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  operator!= (const index_set_t rhs) const
  {
    const bitset_t* pthis = this;
    const bitset_t* pthat = &rhs;
    return *pthis != *pthat;
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
    const bitset_t* pthat = &rhs;
    *pthis ^= *pthat;
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
    index_set<LO,HI> result = lhs;
    return result ^= rhs;
  }

  /// Set intersection: and
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  operator&= (const index_set_t rhs)
  {
    bitset_t* pthis = this;
    const bitset_t* pthat = &rhs;
    *pthis &= *pthat;
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
    index_set<LO,HI> result = lhs;
    return result &= rhs;
  }

  /// Set union: or
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  operator|= (const index_set_t rhs)
  {
    bitset_t* pthis = this;
    const bitset_t* pthat = &rhs;
    *pthis |= *pthat;
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
    index_set<LO,HI> result = lhs;
    return result |= rhs;
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
    return idx < 0 
         ? bitset_t::test(idx-LO) 
         : idx > 0 
         ? bitset_t::test(idx-LO-1)
         : false; 
  }

  /// Test idx for membership: test value of bit idx
  template< >
  inline
  bool
  index_set<DEFAULT_LO,DEFAULT_HI>::
  test(const index_t idx) const
  {
    // Reference: [JA], 1.2.1
    return idx < 0 
         ? bool(bitset_t::to_ulong() & (1UL << (idx-DEFAULT_LO)))
         : idx > 0 
         ? bool(bitset_t::to_ulong() & (1UL << (idx-DEFAULT_LO-1)))
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
  { return bitset_t::count(); }

  /// Cardinality: Number of indices included in set
  template< >
  inline
  index_t
  index_set<DEFAULT_LO,DEFAULT_HI>::
  count() const
  { 
    // Reference: [JA], 1.3
    set_value_t val = bitset_t::to_ulong();
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
    static const index_set_t hi_mask = bitset_t((1UL <<  HI) - 1UL);
    const bitset_t* pthis = this;
    const index_set_t pos_part = (*pthis >> -LO) & hi_mask;
    return pos_part.count();
  }

  /// Number of positive indices included in set
  template< >
  inline
  index_t
  index_set<DEFAULT_LO,DEFAULT_HI>::
  count_pos() const
  {
    const bitset_t* pthis = this;
    const index_set_t pos_part = *pthis >> -DEFAULT_LO;
    return pos_part.count();
  }

  /// Minimum member, or 0 if none
  template<const index_t LO, const index_t HI>
  inline
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

#if (_GLUCAT_BITS_PER_ULONG == 64)
  /// Minimum member, or 0 if none: default for 64-bit cpus
  template< >
  inline
  index_t
  index_set<-32,32>::
  min() const
  {
    // Reference: [JA], 1.3
    unsigned long val = bitset_t::to_ulong();
    if (val == 0)
      return 0;
    else
    {
      val -= val & (val-1); // isolate lowest bit

      index_t idx = -32;
      if (val & 0xffffffff00000000)
        idx  =  1;
      if (val & 0xffff0000ffff0000)
        idx += 16;
      if (val & 0xff00ff00ff00ff00)
        idx +=  8;
      if (val & 0xf0f0f0f0f0f0f0f0)
        idx +=  4;
      if (val & 0xcccccccccccccccc)
        idx +=  2;
      if (val & 0xaaaaaaaaaaaaaaaa)
        idx +=  1;

      return idx;
    }
  }
#endif

#if (_GLUCAT_BITS_PER_ULONG == 32)
  /// Minimum member, or 0 if none: default for 32-bit cpus
  template< >
  inline
  index_t
  index_set<-16,16>::
  min() const
  {
    // Reference: [JA], 1.3
    unsigned long val = bitset_t::to_ulong();
    if (val == 0)
      return 0;
    else
    {
      val -= val & (val-1); // isolate lowest bit

      index_t idx = -16;
      if (val & 0xffff0000)
        idx  =  1;
      if (val & 0xff00ff00)
        idx +=  8;
      if (val & 0xf0f0f0f0)
        idx +=  4;
      if (val & 0xcccccccc)
        idx +=  2;
      if (val & 0xaaaaaaaa)
        idx +=  1;

      return idx;
    }
  }
#endif

  /// Maximum member, or 0 if none
  template<const index_t LO, const index_t HI>
  inline
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

#if (_GLUCAT_BITS_PER_ULONG == 64)
  /// Maximum member, or 0 if none: default for 64-bit cpus
  template< >
  inline
  index_t
  index_set<-32,32>::
  max() const
  {
    // Reference: [JA], 1.6
    unsigned long val = bitset_t::to_ulong();
    if (val == 0)
      return 0;
    else
    {
      index_t idx = -32;
      if (val & 0xffffffff00000000)
        { val >>= 32; idx  =  1; }
      if (val & 0x00000000ffff0000)
        { val >>= 16; idx += 16; }
      if (val & 0x000000000000ff00)
        { val >>=  8; idx +=  8; }
      if (val & 0x00000000000000f0)
        { val >>=  4; idx +=  4; }
      if (val & 0x000000000000000c)
        { val >>=  2; idx +=  2; }
      if (val & 0x0000000000000002)
        {             idx +=  1; }
      return idx;
    }
  }
#endif

#if (_GLUCAT_BITS_PER_ULONG == 32)
  /// Maximum member, or 0 if none: default for 32-bit cpus
  template< >
  inline
  index_t
  index_set<-16,16>::
  max() const
  {
    // Reference: [JA], 1.6
    unsigned long val = bitset_t::to_ulong();
    if (val == 0)
      return 0;
    else
    {
      index_t idx = -16;
      if (val & 0xffff0000)
        { val >>= 16; idx  =  1; }
      if (val & 0x0000ff00)
        { val >>=  8; idx +=  8; }
      if (val & 0x000000f0)
        { val >>=  4; idx +=  4; }
      if (val & 0x0000000c)
        { val >>=  2; idx +=  2; }
      if (val & 0x00000002)
        {             idx +=  1; }
      return idx;
    }
  }
#endif

  /// Lexicographic ordering of two sets: -1 if a<b, +1 if a>b, 0 if a==b
  //  eg. {3,4,5} is less than {3,7,8}
  template<const index_t LO, const index_t HI>
  inline
  int
  compare(const index_set<LO,HI>& a, const index_set<LO,HI>& b)
  {
    for (index_t
        i = LO;
        i <= HI;
        i++)
      if (a[i] != b[i])
        return( a[i] < b[i] ? -1 : +1 );
    return 0;    // all elements are equal => a == b
  }

  /// Lexicographic ordering of two sets: *this < rhs
  //  eg. {3,4,5} is less than {3,7,8}
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  lex_less_than(const index_set_t rhs) const
  {
    for (index_t
        idx = HI;
        idx >= 1;
        --idx)
      if (this->test(idx) != rhs.test(idx))
        return this->test(idx) < rhs.test(idx);
    for (index_t
        idx = -1;
        idx >= LO;
        --idx)
      if (this->test(idx) != rhs.test(idx))
        return this->test(idx) < rhs.test(idx);
    return false;
  }

  /// Less than operator used for comparisons, map, etc.
  // Order by count, then order lexicographically within the equivalence class of count.
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  operator< (const index_set_t rhs) const
  {
    index_t this_grade = this->count();
    index_t rhs_grade  = rhs.count();
    return this_grade < rhs_grade 
         ? true 
         : this_grade > rhs_grade 
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
    char c = 0;
    index_t i = 0;
    index_set<LO,HI> local;
    bool bracketed;
    s >> c;
    if (s.fail() || s.bad())
      return s;
    bracketed = (c == '{');
    if (!bracketed)
      s.putback(c);
    for (s >> i;
        !s.fail() && !s.bad();
        s >> i)
    {
      if ((i < LO) || (i > HI))
      {
        s.clear(std::ios_base::badbit); // set state to error
        break;
      }
      local.set(i);
      s >> c;
      if (!s.fail() && (c != ','))
        s.clear(std::ios_base::failbit); // set state to fail
    }
    if (bracketed && (c != '}'))
      s.clear(std::ios_base::badbit); // set state to error
    else if (!s.bad())
    {
      s.clear();
      ist = local;
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
      { 
        if (this->test(unfold_idx))
          result.set(fold_idx);
        fold_idx--;
      } 
    fold_idx = 1;
    for (unfold_idx = 1;
        unfold_idx <= frm_max;
        ++unfold_idx)
      if (frm.test(unfold_idx))
      { 
        if (this->test(unfold_idx))
          result.set(fold_idx);
        fold_idx++;
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

  /// Sign of product of index sets, using Walsh functions and Gray codes.
  // Reference: [L] Chapter 21, 21.3
  template<const index_t LO, const index_t HI>
  int
  index_set<LO,HI>::
  sign_of_mult(const index_set_t ist) const
  {
    bool h = false;
    bool negative = false;
    for (index_t
        j = LO;
        j != 0;
        ++j)
    {
      h ^= ist.test(j);
      if (this->test(j))
        negative ^= h;
    }
    for (index_t
        j = 1;
        j <= HI;
        ++j)
    {
      const bool b = ist.test(j);
      h ^= b;
      if (this->test(j))
      {
        negative ^= h;
        negative ^= b;
      }
    }
    return negative ? -1 : 1;
  }

  /// Sign of product of index sets, using Walsh functions and Gray codes.
  // Reference: [L] Chapter 21, 21.3
  template< >
  int
  index_set<DEFAULT_LO,DEFAULT_HI>::
  sign_of_mult(const index_set_t ist) const
  {
    // Reference: [JA]
    unsigned long uthis =    bitset_t::to_ulong();
    unsigned long uist = ist.bitset_t::to_ulong();
    unsigned long h = 0;
    unsigned long negative = 0;
    for (index_t
        j = 0;
        j != -DEFAULT_LO;
        ++j)
    {
      h ^= uist >> j;
      negative ^= h & (uthis >> j);
    }
    for (index_t
        j = -DEFAULT_LO;
        j != DEFAULT_HI-DEFAULT_LO;
        ++j)
    {
      const unsigned long b = uist >> j;
      h ^= b;
      negative ^= (h ^ b) & (uthis >> j);
    }
    return (negative & 1) ? -1 : 1;
  }

  /// Hash function
  template<const index_t LO, const index_t HI>
  inline
  size_t
  index_set<LO,HI>::
  hash_fn() const
  {
    static const index_set_t lo_mask = bitset_t((1UL << -LO) - 1);
    const index_set_t neg_part = *this & lo_mask;
    static const index_set_t hi_mask = bitset_t((1UL <<  HI) - 1);
    const index_set_t pos_part = (*this >> -LO) & hi_mask;
    return size_t(neg_part.bitset_t::to_ulong() ^
                  pos_part.bitset_t::to_ulong());
  }

  /// Hash function
  template< >
  inline
  size_t
  index_set<DEFAULT_LO,DEFAULT_HI>::
  hash_fn() const
  {
    static const unsigned long lo_mask = (1UL << -DEFAULT_LO) - 1UL;
    const unsigned long uthis = bitset_t::to_ulong();
    const unsigned long neg_part = uthis & lo_mask;
    const unsigned long pos_part = uthis >> -DEFAULT_LO;
    return size_t(neg_part ^ pos_part);
  }

  /// Square of generator index j
  inline
  int
  sign_of_square(index_t j)
  { return j < 0 ? -1 : 1; }

  /// Maximum positive index, or 0 if none
  template<const index_t LO, const index_t HI>
  inline
  index_t
  max_pos(const index_set<LO,HI>& ist)
  { return std::max(ist.max(), 0); }

  /// Minimum negative index, or 0 if none
  template<const index_t LO, const index_t HI>
  inline
  index_t
  min_neg(const index_set<LO,HI>& ist)
  { return std::min(ist.min(), 0); }

  /// Set containing a range of indices from range_min to range_max
  template<const index_t LO, const index_t HI>
  const index_set<LO,HI>
  index_range(const index_t range_min, const index_t range_max)
  {
    typedef index_set<LO,HI> index_set_t;
    const index_t safe_min = std::max(range_min, LO);
    const index_t safe_max = std::min(range_max, HI);
    index_set_t result;
    for (index_t
        idx = safe_min;
        idx <= safe_max;
        ++idx)
      if (idx != 0)
        result.set(idx);
    return result;
  }

// index_set reference

  /// index_set reference
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::reference::
  reference( index_set<LO,HI>& ist, index_t idx ) :
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
  operator= (const typename index_set<LO,HI>::reference& j)
  {
    if ( (j.m_pst)[j.m_idx] )
      m_pst->set(m_idx);
    else
      m_pst->reset(m_idx);
    return *this;
  }

  /// flips the bit
  template<const index_t LO, const index_t HI>
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
  typename index_set<LO,HI>::reference&
  index_set<LO,HI>::reference::
  flip()
  {
    m_pst->flip(m_idx);
    return *this;
  }
}
#endif // _GLUCAT_INDEX_SET_IMP_H
