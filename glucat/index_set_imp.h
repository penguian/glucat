#ifndef _GLUCAT_INDEX_SET_IMP_H
#define _GLUCAT_INDEX_SET_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    index_set_imp.h : Implement a class for a set of non-zero integer indices
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
  template<const index_t LO, const index_t HI>
  const std::string
  index_set<LO,HI>::
  classname()
  { return "index_set"; }

  /// Constructor from index value
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set(const index_t& idx)
  {
    this->set(idx);
  }

  /// Constructor from bitset_t
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set(const bitset_t& bst)
  {
    const index_set_t* pbst = static_cast<const index_set_t*>(&bst);
    *this = *pbst;
  }

  /// Constructor from set value of an index set folded within the given frame
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set(const set_value_t& folded_val, const index_set_t& frm, const bool prechecked)
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
  operator== (const index_set& rhs) const
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
  operator!= (const index_set& rhs) const
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
  {
    index_set_t result;
    bitset_t* presult = &result;
    const bitset_t* pthis = this;
    *presult = ~(*pthis);
    return result;
  }

  /// Symmetric set difference: exclusive or
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  operator^= (const index_set<LO,HI>& rhs)
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
  operator&= (const index_set<LO,HI>& rhs)
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
  operator|= (const index_set<LO,HI>& rhs)
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

  /// Subscripting: Test idx for membership: test value of bit idx
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  operator[] (index_t idx) const
  { return idx > 0 ? bitset_t::test(idx-LO-1) :
           idx < 0 ? bitset_t::test(idx-LO) : false; }

  /// Subscripting: Element access
  template<const index_t LO, const index_t HI>
  inline
  typename index_set<LO,HI>::reference
  index_set<LO,HI>::
  operator[] (index_t idx)
  { return reference(*this, idx); }

  /// Test idx for membership: test value of bit idx
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  test(index_t idx) const
  { return idx > 0 ? bitset_t::test(idx-LO-1) :
           idx < 0 ? bitset_t::test(idx-LO) : false; }

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
  set(index_t idx, int val)
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
  reset(index_t idx)
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
  flip(index_t idx)
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

  /// Number of negative indices included in set
  template<const index_t LO, const index_t HI>
  inline
  index_t
  index_set<LO,HI>::
  count_neg() const
  {
    index_t result = 0;
    for (size_t 
        idx = 0; 
        idx != -LO; 
        ++idx)
      if (bitset_t::test(idx))
        ++result;
    return result;
  }

  /// Number of positive indices included in set
  template<const index_t LO, const index_t HI>
  inline
  index_t
  index_set<LO,HI>::
  count_pos() const
  {
    index_t result = 0;
    for (size_t 
        idx = HI-LO-1; 
        idx != -LO-1; 
        --idx)
      if (bitset_t::test(idx))
        ++result;
    return result;
  }

  /// Minimum member, or 0 if none
  template<const index_t LO, const index_t HI>
  inline
  index_t
  index_set<LO,HI>::
  min() const
  {
    for (size_t 
        idx = 0; 
        idx != -LO; 
        ++idx)
      if (bitset_t::test(idx))
        return idx+LO;
    for (size_t 
        idx = -LO; 
        idx != HI-LO; 
        ++idx)
      if (bitset_t::test(idx))
        return idx+LO+1;
    return 0;
  }

  /// Maximum member, or 0 if none
  template<const index_t LO, const index_t HI>
  inline
  index_t
  index_set<LO,HI>::
  max() const
  {
    for (size_t 
        idx = HI-LO-1; 
        idx != -LO-1; 
        --idx)
      if (bitset_t::test(idx))
        return idx+LO+1;
    for (size_t 
        idx = -LO; 
        idx != 0; 
        --idx)
      if (bitset_t::test(idx))
        return idx+LO;
    if (bitset_t::test(0))
      return LO;
    return 0;
  }

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
      if(a[i] != b[i])
        return( (a[i] < b[i]) ? -1 : +1 );
    return 0;    // all elements are equal => a == b
  }

  /// Lexicographic ordering of two sets: *this < rhs
  //  eg. {3,4,5} is less than {3,7,8}
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  lex_less_than(const index_set<LO,HI>& rhs) const
  {
    const bitset_t* prhs = &rhs;
    for (size_t 
        idx = 0; 
        idx != HI - LO; 
        ++idx)
      if(bitset_t::test(idx) != prhs->test(idx))
        return bitset_t::test(idx) > prhs->test(idx);
    return false;
  }

  /// Less than operator used for comparisons, map, etc.
  // Order by count, then order lexicographically within the equivalence class of count.
  template<const index_t LO, const index_t HI>
  inline
  bool
  index_set<LO,HI>::
  operator< (const index_set<LO,HI>& rhs) const
  {
    index_t this_grade = this->count();
    index_t rhs_grade  = rhs.count();
    return this_grade < rhs_grade ? true :
           this_grade > rhs_grade ? false:
           this->lex_less_than(rhs);
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
    if(i <= HI)
      os << i;
    for (++i; 
        i <= HI; 
        ++i)
      if(ist[i])
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
  fold(const index_set& frm, const bool prechecked) const
  {
    if (!prechecked && ((*this | frm) != frm))
      throw error_t("fold(frm): cannot fold from outside of frame");
    index_set_t result;
    index_t fold_idx = -1;
    index_t unfold_idx;
    for (unfold_idx = -1; 
        unfold_idx >= LO; 
        --unfold_idx)
      if (frm[unfold_idx])
        result.set(fold_idx--, this->test(unfold_idx));
    fold_idx = 1;
    for (unfold_idx = 1; 
        unfold_idx <= HI; 
        ++unfold_idx)
      if (frm[unfold_idx])
        result.set(fold_idx++, this->test(unfold_idx));
    return result;
  }

  /// Unfold this index set within the given frame
  template<const index_t LO, const index_t HI>
  const
  index_set<LO,HI>
  index_set<LO,HI>::
  unfold(const index_set& frm, const bool prechecked) const
  {
    const char* msg =
      "unfold(frm): cannot unfold into a smaller frame";
    index_set_t result;
    index_t fold_idx = -1;
    index_t unfold_idx;
    for (unfold_idx = -1; 
        unfold_idx >= LO; 
        --unfold_idx)
      if (frm[unfold_idx])
        result.set(unfold_idx, this->test(fold_idx--));
    if (!prechecked && ((fold_idx+1) > this->min()))
      throw error_t(msg);
    fold_idx = 1;
    for (unfold_idx = 1; 
        unfold_idx <= HI; 
        ++unfold_idx)
      if (frm[unfold_idx])
        result.set(unfold_idx, this->test(fold_idx++));
    if (!prechecked && ((fold_idx-1) < this->max()))
      throw error_t(msg);
    return result;
  }

  /// The set value of the fold of this index set within the given frame
  template<const index_t LO, const index_t HI>
  set_value_t
  index_set<LO,HI>::
  value_of_fold(const index_set& frm) const
  {
    const index_set_t folded_set = this->fold(frm);
    const index_t min_index = frm.fold().min();
    if (min_index == 0)
      return 0;
    const index_t skip = min_index > 0 ? 1 : 0;
    const bitset_t* pfolded = &folded_set;
    return ((*pfolded) >> (min_index - skip - LO)).to_ulong();
  }

  template<const index_t LO, const index_t HI>
  inline
  int
  index_set<LO,HI>::
  sign_of_mult(const index_set<LO,HI>& ist) const
  {
    int result = 1;
    index_t i = 0;
    index_t j = 0;

    i = this->count();
    for (j = LO; 
        j <= HI; 
        ++j)
      if (this->test(j))
      {
        i--;
        if (ist[j])
          result *= (i % 2 == 0) ? sign_of_square(j) :
                                  -sign_of_square(j);
      }
      else if (ist[j])
        if (i % 2 == 1)
          result *= -1;
    return result;
  }

  /// Hash function
  template< const index_t LO, const index_t HI>
  inline
  size_t
  index_set<LO,HI>::
  hash_fn() const
  {
    static const index_set_t lo_mask = bitset_t((1UL << sizeof(size_t)) - 1);
    const index_set_t neg_part = *this & lo_mask;
    const index_set_t pos_part = *this >> -LO;
    const bitset_t* pneg_part = &neg_part;
    const bitset_t* ppos_part = &pos_part;
    return size_t((*pneg_part).to_ulong() ^
                  (*ppos_part).to_ulong());
  }

  /// Square of generator index j
  inline
  int
  sign_of_square(index_t j)
  { return j < 0 ? -1 : 1; }

  /// maximum positive index, or 0 if none
  template<const index_t LO, const index_t HI>
  inline
  index_t
  max_pos(const index_set<LO,HI>& ist)
  { return std::max(ist.max(), index_t(0)); }

  /// minimum negative index, or 0 if none
  template<const index_t LO, const index_t HI>
  inline
  index_t
  min_neg(const index_set<LO,HI>& ist)
  { return std::min(ist.min(), index_t(0)); }

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
        result |= index_set_t(idx);
    return result;
  }

// index_set reference

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
