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
  const char*
  index_set<LO,HI>::
  classname()
  { return "index_set"; }

  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set<LO,HI>(const index_t& idx)
  {
    set(idx);
    reset(0);
  }

  /// Constructor from set value of an index set folded within the given frame
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set<LO,HI>(const set_value_t& folded_val, const index_set& frm, const bool prechecked)
  {
    if (!prechecked)
      if (folded_val > set_value_t(1 << frm.count()))
        throw error_t("index_set(val,frm): cannot create: value gives an index set outside of frame");
    const index_set_t folded_frame = frm.fold();
    const index_t min_index = folded_frame.min();
    const index_t max_index = folded_frame.max();
    const index_t skip = min_index > 0 ? 0 : 1;
    index_set<LO,HI> folded_set;
    for (index_t idx = -1; idx >= min_index; --idx)
      if ((folded_val >> (idx - min_index)) & 1)
        folded_set.set(idx);
    for (index_t idx = 1; idx <= max_index; ++idx)
      if ((folded_val >> (idx - skip - min_index)) & 1)
        folded_set.set(idx);
    folded_set.reset(0);
    *this = folded_set.unfold(frm);
  }

  /// Constructor from string
  template<const index_t LO, const index_t HI>
  index_set<LO,HI>::
  index_set<LO,HI>(const string& str)
  {
    istringstream ss(str);
    ss >> *this;
  }

  /// Equality
  template<const index_t LO, const index_t HI>
  inline
  const
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
  const bool
  index_set<LO,HI>::
  operator!= (const index_set& rhs) const
  {
    const bitset_t* pthis = this;
    const bitset_t* pthat = &rhs;
    return *pthis != *pthat;
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
    reset(0);
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
    index_set<LO,HI> result(lhs);
    return result ^= rhs;
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
    reset(0);
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
    index_set<LO,HI> result(lhs);
    return result |= rhs;
  }

  /// Subscripting: Test idx for membership: test value of bit idx
  template<const index_t LO, const index_t HI>
  inline
  const bool
  index_set<LO,HI>::
  operator[] (index_t idx) const
  { return test(idx); }

  /// Subscripting: Element access
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>::reference
  index_set<LO,HI>::
  operator[] (index_t idx)
  {
    return reference(*this, idx);
  }

  /// Test idx for membership: test value of bit idx
  template<const index_t LO, const index_t HI>
  inline
  const bool
  index_set<LO,HI>::
  test(index_t idx) const
  { return bitset_t::test(idx-LO); }

  /// Include all indices except 0: set all bits except 0
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  set()
  {
    bitset_t::set();
    reset(0);
    return *this;
  }

  /// Include idx: Set bit at idx if idx != 0
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  set(index_t idx)
  {
    if(idx != 0)
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
    if(idx != 0)
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
    reset(0);
    return *this;
  }

  /// Complement membership of idx if idx != 0: flip bit at idx if idx != 0
  template<const index_t LO, const index_t HI>
  inline
  index_set<LO,HI>&
  index_set<LO,HI>::
  flip(index_t idx)
  {
    if(idx != 0)
      bitset_t::flip(idx-LO);
  }

  /// Cardinality: Number of indices included in set
  template<const index_t LO, const index_t HI>
  inline
  const index_t
  index_set<LO,HI>::
  count() const
  { return bitset_t::count(); }

  /// Minimum member
  template<const index_t LO, const index_t HI>
  inline
  const index_t
  index_set<LO,HI>::
  min() const
  {
    for (index_t result = LO; result <= HI; ++result)
      if (test(result))
        return result;
    return 0;
  }

  /// Maximum member
  template<const index_t LO, const index_t HI>
  inline
  const index_t
  index_set<LO,HI>::
  max() const
  {
    for (index_t result = HI; result >= LO; --result)
      if (test(result))
        return result;
    return 0;
  }

  /// Lexicographic ordering of two sets: -1 if a<b, +1 if a>b, 0 if a==b
  //  eg. {3,4,5} is less than {3,7,8}
  template<const index_t LO, const index_t HI>
  inline
  const int
  compare(const index_set<LO,HI>& a, const index_set<LO,HI>& b)
  {
    for(index_t i = LO; i <= HI; i++)
      if(a[i] != b[i])
        return( (a[i] < b[i]) ? -1 : +1 );
    return 0;    // all elements are equal => a == b
  }

  /// Less than operator used for comparisons, map, etc.
  // Order by count, then order lexicographically within the equivalence class of count.
  template<const index_t LO, const index_t HI>
  inline
  const bool
  operator< (const index_set<LO,HI>& a, const index_set<LO,HI>& b)
  {
    index_t a_grade = a.count();
    index_t b_grade = b.count();
    return a_grade < b_grade ? true :
           a_grade > b_grade ? false:
           compare( a, b ) > 0 ? true : false ;
  }

  /// Write out index set
  template<const index_t LO, const index_t HI>
  ostream&
  operator<< (ostream& os, const index_set<LO,HI>& ist)
  {
    index_t i;
    os << '{';
    for(i= LO; (i <= HI) & !(ist[i]); ++i)
    { }
    if(i <= HI)
      os << i;
    for(++i ; i <= HI; ++i)
      if(ist[i])
        os << ',' << i;
    os << '}';
    return os;
  }

  /// Read in index set
  template<const index_t LO, const index_t HI>
  istream&
  operator>> (istream& s, index_set<LO,HI>& ist)
  {
    char c = 0;
    index_t i = 0;
    index_set<LO,HI> local;
    bool bracketed;
    s >> c;
    bracketed = (c == '{');
    if (!bracketed)
      s.putback(c);
    for (s >> i >> c; c == ','; s >> i >> c)
    {
      if ((i < LO) || (i > HI))
      {
        s.clear(std::ios_base::failbit); // set state to error
        break;
      }
      local.set(i);
    }
    if (bracketed && (c != '}'))
      s.clear(ios_base::failbit); // set state to error
    else if (!s.fail())
      if ((i < LO) || (i > HI))
        s.clear(ios_base::failbit); // set state to error
      else
      {
        local.set(i);
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
  { return fold(*this, true); }

  /// Fold this index set within the given frame
  template<const index_t LO, const index_t HI>
  const
  index_set<LO,HI>
  index_set<LO,HI>::
  fold(const index_set& frm, const bool prechecked) const
  {
    if (!prechecked && ((*this | frm) != frm))
      throw error_t("fold(frm): cannot fold from outside of frame");
    index_set<LO,HI> result;
    index_t fold_idx = -1;
    index_t unfold_idx;
    for (unfold_idx = -1; unfold_idx >= LO; --unfold_idx)
      if (frm[unfold_idx])
      {
        result[fold_idx] = (*this)[unfold_idx];
        --fold_idx;
      }
    fold_idx = 1;
    for (unfold_idx = 1; unfold_idx <= HI; ++unfold_idx)
      if (frm[unfold_idx])
      {
        result[fold_idx] = (*this)[unfold_idx];
        ++fold_idx;
      }
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
    for (unfold_idx = -1; unfold_idx >= LO; --unfold_idx)
      if (frm[unfold_idx])
      {
        result[unfold_idx] = (*this)[fold_idx];
        --fold_idx;
      }
    if (!prechecked && ((fold_idx+1) > this->min()))
      throw error_t(msg);
    fold_idx = 1;
    for (unfold_idx = 1; unfold_idx <= HI; ++unfold_idx)
      if (frm[unfold_idx])
      {
        result[unfold_idx] = (*this)[fold_idx];
        ++fold_idx;
      }
    if (!prechecked && ((fold_idx-1) < this->max()))
      throw error_t(msg);
    return result;
  }

  /// The set value of the fold of this index set within the given frame
  template<const index_t LO, const index_t HI>
  const set_value_t
  index_set<LO,HI>::
  value_of_fold(const index_set& frm) const
  {
    if (frm.count() > BITS_PER_SET_VALUE)
      return 0;
    const index_set<LO,HI> folded_set = fold(frm);
    const index_t min_index = frm.fold().min();
    if (min_index == HI+1)
      return 0;
    const index_t skip = min_index > 0 ? 0 : 1;
    set_value_t result = 0;
    for (index_t idx = -1; idx >= min_index; --idx)
      if (folded_set[idx])
        result |= 1 << (idx - min_index);
    for (index_t idx = 1; idx <= HI; ++idx)
      if (folded_set[idx])
        result |= 1 << (idx - skip - min_index);
    return result;
  }

  template<const index_t LO, const index_t HI>
  inline
  const int
  index_set<LO,HI>::
  sign_of_mult(const index_set<LO,HI>& ist) const
  {
    int result = 1;
    index_t i = 0;
    index_t j = 0;

    i = count();
    for(j = LO; j <= HI; ++j)
      if (test(j))
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

  /// Square of generator index j
  inline
  const int
  sign_of_square(index_t j)
  { return j < 0 ? -1 : 1; }

  /// maximum positive index, or 0 if none
  template<const index_t LO, const index_t HI>
  inline
  const index_t
  max_pos(const index_set<LO,HI>& ist)
  { return max(ist.max(), index_t(0)); }

  /// minimum negative index, or 0 if none
  template<const index_t LO, const index_t HI>
  inline
  const index_t
  min_neg(const index_set<LO,HI>& ist)
  { return min(ist.min(), index_t(0)); }

  /// Set containing a range of indices from range_min to range_max
  template<const index_t LO, const index_t HI>
  const index_set<LO,HI>
  index_range(const index_t range_min, const index_t range_max)
  {
    typedef index_set<LO,HI> index_set_t;
    const index_t safe_min = max(range_min, LO);
    const index_t safe_max = min(range_max, HI);
    index_set_t result;
    for (index_t idx = safe_min; idx <= safe_max; ++idx)
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
  index_set<LO,HI>::reference&
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
  index_set<LO,HI>::reference&
  index_set<LO,HI>::reference::
  operator= (const index_set<LO,HI>::reference& j)
  {
    if ( (j.m_pst)[j.m_idx] )
      m_pst->set(m_idx);
    else
      m_pst->reset(m_idx);
    return *this;
  }

  /// flips the bit
  template<const index_t LO, const index_t HI>
  const bool
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
  index_set<LO,HI>::reference&
  index_set<LO,HI>::reference::
  flip()
  {
    m_pst->flip(m_idx);
    return *this;
  }
}
#endif // _GLUCAT_INDEX_SET_IMP_H
