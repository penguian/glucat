#ifndef _GLUCAT_FRAMED_MULTI_IMP_H
#define _GLUCAT_FRAMED_MULTI_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    framed_multi_imp.h : Implement the coordinate map representation of a
    Clifford algebra element
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
  /// Class name used in messages
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const std::string
  framed_multi<Scalar_T,LO,HI>::
  classname()
  { return "framed_multi"; }

  /// Default constructor
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi()
  { }

  /// Construct a multivector, within a given frame, from a given multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const multivector_t& val,
               const index_set_t frm, const bool prechecked)
  {
    if (!prechecked && (val.frame() | frm) != frm)
      throw error_t("multivector_t(val,frm): cannot initialize with value outside of frame");
    *this = val;
  }

  /// Construct a multivector from an index set and a scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const index_set_t ist, const Scalar_T& crd)
  {
    if (crd != Scalar_T(0))
      this->insert(term_t(ist, crd));
  }

  /// Construct a multivector, within a given frame, from an index set and a scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const index_set_t ist, const Scalar_T& crd,
               const index_set_t frm, const bool prechecked)
  {
    if (!prechecked && (ist | frm) != frm)
      throw error_t("multivector_t(ist,crd,frm): cannot initialize with value outside of frame");
    if (crd != Scalar_T(0))
      this->insert(term_t(ist, crd));
  }

  /// Construct a multivector from a scalar (within a frame, if given)
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const Scalar_T& scr, const index_set_t frm)
  {
    if (scr != Scalar_T(0))
      this->insert(term_t(index_set_t(), scr));
  }

  /// Construct a multivector from an int (within a frame, if given)
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const int scr, const index_set_t frm)
  {
    if (scr != 0)
      this->insert(term_t(index_set_t(), Scalar_T(scr)));
  }

  /// Construct a multivector, within a given frame, from a given vector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const vector_t& vec,
               const index_set_t frm, const bool prechecked)
  {
    if (!prechecked && index_t(vec.size()) != frm.count())
      throw error_t("multivector_t(vec,frm): cannot initialize with vector not matching frame");
    typename vector_t::const_iterator vec_it = vec.begin();
    const index_t begin_index = frm.min();
    const index_t end_index = frm.max()+1;
    for (index_t
        idx = begin_index;
        idx != end_index;
        ++idx)
      if (frm[idx])
      {
        *this += term_t(index_set_t(idx), *vec_it);
        ++vec_it;
      }
  }

  /// Construct a multivector from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const std::string& str)
  {
    std::istringstream ss(str);
    ss >> *this;
  }

  /// Construct a multivector, within a given frame, from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const std::string& str, const index_set_t frm, const bool prechecked)
  {
    if (prechecked)
      *this = multivector_t(str);
    else
      *this = multivector_t(multivector_t(str), frm, false);
  }

  /// Construct a multivector from a matrix_multi_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const matrix_multi_t& val)
  {
    if (val.m_matrix.size1() >= Tune_P::inv_fast_dim_threshold)
      try
      {
        *this = val.fast_framed_multi();
        return;
      }
      catch (const glucat_error& e)
      { }
    const index_set_t frm = val.frame();
    const set_value_t end_set_value = 1 << frm.count();
    for (set_value_t
        stv = 0;
        stv != end_set_value;
        stv++)
    {
      const index_set_t ist = index_set_t(stv, frm, true);
      const Scalar_T& crd =
        matrix::inner<Scalar_T>(val.basis_element(ist), val.m_matrix);
      if (crd != Scalar_T(0))
        this->insert(term_t(ist, crd));
    }
  }

  /// Test for equality of multivectors
  template< typename Scalar_T, const index_t LO, const index_t HI >
  bool
  framed_multi<Scalar_T,LO,HI>::
  operator==  (const multivector_t& rhs) const
  {
    if (this->size() != rhs.size())
      return false;
    else if (compare_types< map_t, sorted_map_t >::are_same)
    {
      const_iterator this_it = this->begin();
      const_iterator rhs_it = rhs.begin();
      for (;
          (this_it != this->end()) || (rhs_it != rhs.end());
          this_it++, rhs_it++)
        if (*this_it != *rhs_it)
          return false;
      return (this_it == this->end()) && (rhs_it == rhs.end());
    }
    else
    {
      for (const_iterator
          this_it = this->begin();
          this_it != this->end();
          this_it++)
      {
        const const_iterator& rhs_it = rhs.find(this_it->first);
        if (rhs_it == rhs.end())
          return false;
        else if (rhs_it->second != this_it->second)
          return false;
      }
      return true;
    }
  }

  /// Test for equality of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  bool
  framed_multi<Scalar_T,LO,HI>::
  operator==  (const Scalar_T& scr) const
  { 
    switch (this->size())
    {
    case 0:
      return scr == Scalar_T(0);
    case 1:
      {
        const const_iterator& this_it = this->begin();
        return this_it->first == index_set_t() && this_it->second == scr;
      }
    default:
      return false;
    }
  }


  /// Geometric sum of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator+= (const Scalar_T& scr)
  {
    *this += term_t(index_set_t(),scr);
    return *this;
  }

  /// Geometric sum
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator+= (const multivector_t& rhs)
  { // simply add terms
    for (const_iterator
        rhs_it = rhs.begin();
        rhs_it != rhs.end();
        ++rhs_it)
      *this += *rhs_it;
    return *this;
  }

  /// Geometric difference
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator-= (const multivector_t& rhs)
  {
    for (const_iterator
        rhs_it = rhs.begin();
        rhs_it != rhs.end();
        ++rhs_it)
      *this += term_t(rhs_it->first, -(rhs_it->second));
    return *this;
  }

  /// Unary -
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  operator- () const
  { return *this * Scalar_T(-1); }

  /// Product of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator*= (const Scalar_T& scr)
  { // multiply coordinates of all terms by scalar
    if (numeric_traits<Scalar_T>::isNaN_or_isInf(scr))
      return *this = numeric_traits<Scalar_T>::NaN();
    if (scr == Scalar_T(0))
      if (this->isnan())
        *this = numeric_traits<Scalar_T>::NaN();
      else
        this->clear();
    else
      for (iterator
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
        this_it->second *= scr;
    return *this;
  }

  /// Geometric product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator*= (const multivector_t& rhs)
  {
    const index_set_t our_frame = this->frame() | rhs.frame();
    const index_t frm_count = our_frame.count();
    const set_value_t algebra_dim = 1 << frm_count;
    const bool direct_mult = double(this->size())*rhs.size() <= double(algebra_dim)
#ifdef _GLUCAT_USE_GNU_CXX_HASH_MAP
                           || frm_count < Tune_P::mult_matrix_threshold
#endif
                           ;
    if (direct_mult)
    { // If we have a sparse multiply, store the result directly
      multivector_t  result;
      for (const_iterator
          this_it = this->begin();
          this_it != this->end();
          this_it++)
        for (const_iterator
            rhs_it = rhs.begin();
            rhs_it != rhs.end();
            rhs_it++)
          result += (*this_it) * (*rhs_it);
      return *this = result;
    }
#ifndef _GLUCAT_USE_GNU_CXX_HASH_MAP
    else if (frm_count < Tune_P::mult_matrix_threshold)
    { // Fastest dense algorithm in low dimensions stores result in array
      std::vector< Scalar_T > result_array( algebra_dim, Scalar_T(0) );
      for (const_iterator
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
        for (const_iterator
            rhs_it = rhs.begin();
            rhs_it != rhs.end();
            ++rhs_it)
        {
          const term_t& term = (*this_it) * (*rhs_it);
          const set_value_t stv = term.first.value_of_fold(our_frame);
          result_array[stv] += term.second;
        }
      this->clear();
      for (set_value_t
          stv = 0;
          stv < algebra_dim;
          ++stv)
        if (result_array[stv] != Scalar_T(0))
        {
          const index_set_t result_ist = index_set_t(stv, our_frame);
          this->insert(term_t(result_ist, result_array[stv]));
        }
      return *this;
    }
#endif
    else
      // Past a certain threshold, the matrix algorithm is fastest
      return *this = matrix_multi_t(*this, our_frame, true) *
                     matrix_multi_t(rhs, our_frame, true);
  }

  /// Left contraction
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator%= (const multivector_t& rhs)
  {
    // Reference: Leo Dorst, "Honing geometric algebra for its use in the computer sciences",
    // in Geometric Computing with Clifford Algebras, ed. G. Sommer,
    // Springer 2001, Chapter 6, pp. 127-152.
    // http://staff.science.uva.nl/~leo/clifford/index.html

    multivector_t result;
    for (const_iterator
        this_it = this->begin();
        this_it != this->end();
        this_it++)
      for (const_iterator
          rhs_it = rhs.begin();
          rhs_it != rhs.end();
          rhs_it++)
        if ((this_it->first | rhs_it->first) == rhs_it->first)
          result += (*this_it) * (*rhs_it);
    return *this = result;
  }

  /// Inner product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator&= (const multivector_t& rhs)
  { // Arvind Raja's original reference:
    // "old clical, innerproduct(p,q:pterm):pterm in file compmod.pas"
    multivector_t result;
    for (const_iterator
        this_it = this->begin();
        this_it != this->end();
        this_it++)
      for (const_iterator
          rhs_it = rhs.begin();
          rhs_it != rhs.end();
          rhs_it++)
      {
        const index_t this_grade = this_it->first.count();
        const index_t rhs_grade = rhs_it->first.count();
        if (this_grade > 0 && rhs_grade > 0)
        {
          const index_t grade = std::abs(this_grade - rhs_grade);
          const term_t& term = (*this_it) * (*rhs_it);
          if (term.first.count() == grade)
            result += term;
        }
      }
    return *this = result;
  }

  /// Outer product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator^= (const multivector_t& rhs)
  { // Arvind Raja's original reference:
    // "old clical, outerproduct(p,q:pterm):pterm in file compmod.pas"
    multivector_t result;
    for (const_iterator
        this_it = this->begin();
        this_it != this->end();
        this_it++)
      for (const_iterator
          rhs_it = rhs.begin();
          rhs_it != rhs.end();
          rhs_it++)
      {
        const index_t grade = this_it->first.count() + rhs_it->first.count();
        if (grade <= HI-LO)
        {
          const term_t& term = (*this_it) * (*rhs_it);
          if (term.first.count() == grade)
            result += term;
        }
      }
    return *this = result;
  }

  /// Quotient of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator/= (const Scalar_T& scr)
  { // Divide coordinates of all terms by scr
    if (numeric_traits<Scalar_T>::isNaN(scr))
      return *this = numeric_traits<Scalar_T>::NaN();
    if (numeric_traits<Scalar_T>::isInf(scr))
      if (this->isnan())
        *this = numeric_traits<Scalar_T>::NaN();
      else
        this->clear();
    else
      for (iterator
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
        this_it->second /= scr;
    return *this;
  }

  /// Geometric quotient
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator/= (const multivector_t& rhs)
  {
    const index_set_t our_frame = this->frame() | rhs.frame();
    matrix_multi_t result = matrix_multi_t(*this, our_frame, true);
    result /= matrix_multi_t(rhs, our_frame, true);
    return *this = result;
  }

  /// Clifford multiplicative inverse
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  inv() const
  {
    matrix_multi_t result = matrix_multi_t(Scalar_T(1), this->frame());
    return result /= matrix_multi_t(*this);
  }

  /// Subscripting: map from index set to scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  Scalar_T
  framed_multi<Scalar_T,LO,HI>::
  operator[] (const index_set_t ist) const
  {
      const const_iterator& this_it = this->find(ist);
      if (this_it == this->end())
        return Scalar_T(0);
      else
        return this_it->second;
  }

  /// Main involution, each {i} is replaced by -{i} in each term
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  involute() const
  {
    multivector_t result = *this;
    for (iterator
        result_it = result.begin();
        result_it != result.end();
        ++result_it)
    { // for a k-vector u, involute(u) == (-1)^k * u
      if ((result_it->first.count() % 2) == 1)
        result_it->second *= Scalar_T(-1);
    }
    return result;
  }

  /// Reversion, order of {i} is reversed in each term
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  reverse() const
  {
    multivector_t result = *this;
    for (iterator
        result_it = result.begin();
        result_it != result.end();
        ++result_it)
      // For a k-vector u, reverse(u) = { -u, k == 2,3 (mod 4)
      //                                {  u, k == 0,1 (mod 4)
      switch (result_it->first.count() % 4)
      {
      case 2:
      case 3:
        result_it->second *= Scalar_T(-1);
        break;
      default:
        break;
      }
    return result;
  }

  /// Conjugation, conj == reverse o involute == involute o reverse
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  conj() const
  {
    multivector_t result = *this;
    for (iterator
        result_it = result.begin();
        result_it != result.end();
        ++result_it)
      // For a k-vector u, conj(u) = { -u, k == 1,2 (mod 4)
      //                             {  u, k == 0,3 (mod 4)
      switch (result_it->first.count() % 4)
      {
      case 1:
      case 2:
        result_it->second *= Scalar_T(-1);
        break;
      default:
        break;
      }
    return result;
  }

  /// Quadratic form := scalar part of rev(x)*x
  template< typename Scalar_T, const index_t LO, const index_t HI >
  Scalar_T
  framed_multi<Scalar_T,LO,HI>::
  quad() const
  {
    // scalar(conj(x)*x) = 2*quad(even(x)) - quad(x)
    // ref: old clical: quadfunction(p:pter):pterm in file compmod.pas
    Scalar_T result = Scalar_T(0);
    for (const_iterator
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
    {
      int sign = 1;
      for (index_t
          i = LO;
          i < 0;
          ++i)
        sign *= (this_it->first.test(i)) ? -1 : 1;
      result += sign * (this_it->second) * (this_it->second);
    }
    return result;
  }

  /// Norm squared := sum of norm squared of coordinates
  template< typename Scalar_T, const index_t LO, const index_t HI >
  Scalar_T
  framed_multi<Scalar_T,LO,HI>::
  norm() const
  {
    Scalar_T result = Scalar_T(0);
    for (const_iterator
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
    {
      const Scalar_T& n2 = ublas::type_traits<Scalar_T>::norm_2(this_it->second);
      result +=  n2 * n2;
    }
    return result;
  }

  /// *this to the m
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  pow(int m) const
  { return glucat::pow(*this, m); }

  /// Outer product power
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  outer_pow(int m) const
  {
    if (m < 0)
      throw error_t("outer_pow(int): negative exponent");
    multivector_t result = Scalar_T(1);
    multivector_t a = *this;
    for (;
        m != 0;
        m >>= 1, a ^= a)
      if (m & 1)
        result ^= a;
    return result;
  }

  /// Grading: part where each term is a grade-vector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  operator() (index_t grade) const
  {
    if ((grade < 0) || (grade > HI-LO))
      return Scalar_T(0);
    else
    {
      multivector_t result;
      for (const_iterator
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
        if (this_it->first.count() == grade)
          result += *this_it;
      return result;
    }
  }

  /// Even part, sum of the even grade terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  even() const
  { // even part of x, sum of the pure(count) with even count
    multivector_t result;
    for (const_iterator
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
      if ((this_it->first.count() % 2) == 0)
        result.insert(*this_it);
    return result;
  }

  /// Odd part, sum of the odd grade terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  odd() const
  { // even part of x, sum of the pure(count) with even count
    multivector_t result;
    for (const_iterator
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
      if ((this_it->first.count() % 2) == 1)
        result.insert(*this_it);
    return result;
  }

  /// Vector part of multivector, as a vector_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const typename framed_multi<Scalar_T,LO,HI>::vector_t
  framed_multi<Scalar_T,LO,HI>::
  vector_part() const
  {
    const index_set_t frm = this->frame();
    vector_t result;
    result.reserve(frm.count());
    const index_t frm_end = frm.max()+1;
    for (index_t
        idx  = frm.min();
        idx != frm_end;
        ++idx)
      // Frame may contain indices which do not correspond to a grade 1 term but
      // frame cannot omit any index corresponding to a grade 1 term
      if (frm[idx])
        result.push_back((*this)[index_set_t(idx)]);
    return result;
  }

  /// Write multivector to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  void
  framed_multi<Scalar_T,LO,HI>::
  write(const std::string& msg) const
  { std::cout << msg << std::endl << "  " << (*this) << std::endl; }

  /// Write multivector to file
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  void
  framed_multi<Scalar_T,LO,HI>::
  write(std::ofstream& ofile, const std::string& msg) const
  {
    if (!ofile)
      throw error_t("write(ofile,msg): cannot write to output file");
    ofile << msg << std::endl << "  " << (*this) << std::endl;
  }

  /// Sorted range for use with output
  template< typename Map_T,typename Sorted_Map_T >
  class sorted_range
  {
  public:
    typedef Map_T map_t;
    typedef Sorted_Map_T sorted_map_t;
    typedef typename Sorted_Map_T::const_iterator sorted_iterator;

    sorted_range (Sorted_Map_T &sorted_val, const Map_T& val)
    {
      for (typename map_t::const_iterator
          val_it = val.begin();
          val_it != val.end();
          ++val_it)
        sorted_val.insert(*val_it);
      sorted_begin = sorted_val.begin();
      sorted_end   = sorted_val.end();
    }
    sorted_iterator sorted_begin;
    sorted_iterator sorted_end;
  };

  template< typename Sorted_Map_T >
  class sorted_range< Sorted_Map_T, Sorted_Map_T >
  {
  public:
    typedef Sorted_Map_T map_t;
    typedef Sorted_Map_T sorted_map_t;
    typedef typename Sorted_Map_T::const_iterator sorted_iterator;

    sorted_range (Sorted_Map_T &sorted_val, const Sorted_Map_T& val)
    : sorted_begin( val.begin() ),
      sorted_end( val.end() )
    { }
    sorted_iterator sorted_begin;
    sorted_iterator sorted_end;
  };

  /// Write multivector to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::ostream&
  operator<< (std::ostream& os, const framed_multi<Scalar_T,LO,HI>& val)
  {
    if (val.empty())
      os << 0;
    else if (val.isnan())
      os << std::numeric_limits<Scalar_T>::quiet_NaN();
    else
    {
      typedef framed_multi<Scalar_T,LO,HI>  multivector_t;
      typedef typename multivector_t::map_t map_t;
      typedef typename multivector_t::sorted_map_t sorted_map_t;
      typedef typename sorted_map_t::const_iterator sorted_iterator;
      sorted_map_t sorted_val;
      sorted_range< map_t, sorted_map_t > sorted_val_range(sorted_val, val);
      sorted_iterator sorted_it = sorted_val_range.sorted_begin;
      os << *sorted_it;
      for (++sorted_it;
          sorted_it != sorted_val_range.sorted_end;
          ++sorted_it)
      {
        Scalar_T scr = sorted_it->second;
        if ((scr != std::real(scr)) || (std::real(scr) >= 0.0))
          os << '+';
        os << *sorted_it;
      }
    }
    return os;
  }

  /// Write term to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::ostream&
  operator<< (std::ostream& os, const std::pair< const index_set<LO,HI>, Scalar_T >& term)
  {
    if (term.first.count() == 0)
      os << term.second;
    else if (term.second == Scalar_T(-1))
    {
      os << '-';
      os << term.first;
    }
    else if (term.second != Scalar_T(1))
    {
      os << term.second;
      os << term.first;
    }
    else
      os << term.first;
    return os;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::istream&
  operator>> (std::istream& s, framed_multi<Scalar_T,LO,HI> & val)
  { // Input looks like 1.0 -2.0{1,2} +3.2{3,4}
    typedef framed_multi<Scalar_T,LO,HI> multivector_t;
    multivector_t local;
    int c = 0;
    while (s)
    {  // Default coordinate value is Scalar_T(0)
      Scalar_T     coordinate = Scalar_T(0);
      // Default index set is empty
      index_set<LO,HI> ist;
      ist.reset();
      c = s.peek();
      // Look for a leading '-'
      const bool negative = (c == int('-'));
      if (negative)
      { // consume the '-'
        c = s.get();
        if (s)
          c = s.peek();
      }
      // Look for a leading '{'
      if (c == int('{') && s)
        // Index set without explicit coordinate.
        // Default is Scalar_T(1)
        coordinate = Scalar_T(1);
      else if (s)
        // Try to read a coordinate value
        s >> coordinate;
      // coordinate is now Scalar_T(0), Scalar_T(1) or a Scalar_T value
      // Expect an index set or a '+' or a '-'
      if (s)
      {
        c = s.peek();
        if (c == int('{') && s)
          // Try to read index set.
          s >> ist;
      }
      if (!s.bad() && coordinate != Scalar_T(0))
      { // Insert whatever we have
        coordinate = negative ? -coordinate : coordinate;
        typedef typename multivector_t::term_t term_t;
        local += term_t(ist, coordinate);
      }
      if (s)
      {
        c = s.peek();
        if (c == int('+'))
          // consume the '+'
          c = s.get();
      }
    }
    if (!s.bad())
      // If s.bad() then we have a corrupt input
      // otherwise we are fine and can copy the resulting framed_multi
      val = local;
    return s;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator+= (const term_t& term)
  { // Insert a term into a multivector_t, add terms with same index set.
    // Do not insert terms with 0 coordinate
    if (term.second != Scalar_T(0))
    {
      const iterator& this_it = this->find(term.first);
      if (this_it == this->end())
        this->insert(term);
      else if (this_it->second + term.second == Scalar_T(0))
        // Erase term if resulting coordinate is 0
        this->erase(this_it);
      else
        this_it->second += term.second;
    }
    return *this;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const index_set<LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  frame() const
  {
    index_set_t result;
    for (const_iterator
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
      result |= this_it->first;
    return result;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  Scalar_T
  framed_multi<Scalar_T,LO,HI>::
  max_abs() const
  {
    Scalar_T result = Scalar_T(0);
    for (const_iterator
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
    {
      const Scalar_T& abs_term = std::abs(this_it->second);
      if (abs_term > result)
        result = abs_term;
    }
    return result;
  }

  /// Check if a multivector contains any IEEE NaN values
  template< typename Scalar_T, const index_t LO, const index_t HI >
  bool
  framed_multi<Scalar_T,LO,HI>::
  isnan() const
  {
    if (std::numeric_limits<Scalar_T>::has_quiet_NaN)
      for (const_iterator
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
          if (numeric_traits<Scalar_T>::isNaN(this_it->second))
            return true;
    return false;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  truncated(const Scalar_T& limit) const
  {
    if (this->isnan())
      return *this;
    Scalar_T top = max_abs();
    multivector_t result;
    if (top != Scalar_T(0))
      for (const_iterator
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
        if (std::abs(this_it->second / top) > limit)
          result.insert(term_t(this_it->first, this_it->second));
    return result;
  }

  /// Subalgebra isomorphism: fold each term within the given frame
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  fold(const index_set_t frm) const
  {
    if (frm.is_contiguous())
      return *this;
    else
    {
      multivector_t result;
      for (const_iterator
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
        result.insert(term_t(this_it->first.fold(frm), this_it->second));
      return result;
    }
  }

  /// Subalgebra isomorphism: unfold each term within the given frame
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  unfold(const index_set_t frm) const
  {
    if (frm.is_contiguous())
      return *this;
    else
    {
      multivector_t result;
      for (const_iterator
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
        result.insert(term_t(this_it->first.unfold(frm, true), this_it->second));
      return result;
    }
  }

  /// Subalgebra isomorphism: R_{p,q} to R_{p-4,q+4}
  // Reference: [L] 16.4 Periodicity of 8, p216
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>&
  framed_multi<Scalar_T,LO,HI>::
  centre_pm4_qp4(index_t& p, index_t& q)
  {
    // We add 4 to q by subtracting 4 from p
    if (q+4 > -LO)
      throw error_t("centre_pm4_qp4(p,q): LO is too high to represent this value");
    if (this->frame().max() > p-4)
    {
      const index_set_t pm3210 =
        index_set_t(p-3)  | index_set_t(p-2)  | index_set_t(p-1)  | index_set_t(p);
      const index_set_t qm4321 =
        index_set_t(-q-4) | index_set_t(-q-3) | index_set_t(-q-2) | index_set_t(-q-1);
      const term_t& tqm4321 = term_t(qm4321, Scalar_T(1));
      multivector_t result;
      for (const_iterator
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
      {
        index_set_t ist = this_it->first;
        if (ist.max() > p-4)
        {
          var_term_t term;
          for (index_t
              n = 0;
              n != 4;
              ++n)
            if (ist[n+p-3])
              term *= term_t(index_set_t(n-q-4), Scalar_T(1)) * tqm4321;
          // Mask out {p-3}..{p}
          result.insert(term_t(ist & ~pm3210, this_it->second) *
                        term_t(term.first, term.second));
        }
        else
          result.insert(*this_it);
      }
      *this = result;
    }
    p -=4; q += 4;
    return *this;
  }

  /// Subalgebra isomorphism: R_{p,q} to R_{p+4,q-4}
  // Reference: [L] 16.4 Periodicity of 8, p216
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>&
  framed_multi<Scalar_T,LO,HI>::
  centre_pp4_qm4(index_t& p, index_t& q)
  {
    // We add 4 to p by subtracting 4 from q
    if (p+4 > HI)
      throw error_t("centre_pp4_qm4(p,q): HI is too low to represent this value");
    if (this->frame().min() < -q+4)
    {
      const index_set_t qp0123 =
        index_set_t(-q)  | index_set_t(-q+1) | index_set_t(-q+2) | index_set_t(-q+3);
      const index_set_t pp1234 =
        index_set_t(p+1) | index_set_t(p+2)  | index_set_t(p+3)  | index_set_t(p+4);
      const term_t& tpp1234 = term_t(pp1234, Scalar_T(1));
      multivector_t result;
      for (const_iterator
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
      {
        index_set_t ist = this_it->first;
        if (ist.min() < -q+4)
        {
          var_term_t term;
          for (index_t
              n = 0;
              n != 4;
              ++n)
            if (ist[n-q])
              term *= term_t(index_set_t(n+p+1), Scalar_T(1)) * tpp1234;
          // Mask out {-q}..{-q+3}
          result.insert(term_t(term.first, term.second) *
                        term_t(ist & ~qp0123, this_it->second));
        }
        else
          result.insert(*this_it);
      }
      *this = result;
    }
    p +=4; q -= 4;
    return *this;
  }

  /// Subalgebra isomorphism: R_{p,q} to R_{q+1,p-1}
  // Reference: [P] Proposition 15.20, p 131
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>&
  framed_multi<Scalar_T,LO,HI>::
  centre_qp1_pm1(index_t& p, index_t& q)
  {
    if (q+1 > HI)
      throw error_t("centre_qp1_pm1(p,q): HI is too low to represent this value");
    if (p-1 > -LO)
      throw error_t("centre_qp1_pm1(p,q): LO is too high to represent this value");
    const index_set_t qp1 = index_set_t(q+1);
    const term_t& tqp1 = term_t(qp1, Scalar_T(1));
    multivector_t result;
    for (const_iterator
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
    {
      const index_set_t ist = this_it->first;
      var_term_t term = var_term_t(index_set_t(), this_it->second);
      for (index_t
          n = -q;
          n != p;
          ++n)
        if (n != 0 && ist[n])
          term *= term_t(index_set_t(-n) | qp1, Scalar_T(1));
      if (p != 0 && ist[p])
        term *= tqp1;
      result.insert(term_t(term.first, term.second));
    }
    index_t orig_p = p;
    p = q+1;
    q = orig_p-1;
    return *this = result;
  }

  /// Divide multivector into quotient with terms divisible by index set, and remainder
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const std::pair< const framed_multi<Scalar_T,LO,HI>, const framed_multi<Scalar_T,LO,HI> >
  framed_multi<Scalar_T,LO,HI>::
  divide(const index_set_t ist) const
  {
    multivector_t quo;
    multivector_t rem;
    for (const_iterator
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
      if ((this_it->first | ist) == this_it->first)
        quo.insert(term_t(this_it->first ^ ist, this_it->second));
      else
        rem.insert(*this_it);
    return framed_pair_t(quo, rem);
  }

  /// Generalized FFT from framed_multi_t to matrix_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const typename framed_multi<Scalar_T,LO,HI>::matrix_t
  framed_multi<Scalar_T,LO,HI>::
  fast(const index_t level, const bool odd) const
  {
    // Assume val is already folded and centred
    if (level == 0)
      return matrix::unit<matrix_t>(1) * scalar(*this);

    typedef typename matrix_multi_t::basis_matrix_t basis_matrix_t;
    const basis_matrix_t& I = matrix::unit<basis_matrix_t>(2);
    basis_matrix_t J(2,2);
    J(0,1)  = Scalar_T(-1);
    J(1,0)  = Scalar_T( 1);
    basis_matrix_t K(2,2);
    K(0,1)  = Scalar_T( 1);
    K(1,0)  = Scalar_T( 1);
    basis_matrix_t JK(2,2);
    JK(0,0) = Scalar_T(-1);
    JK(1,1) = Scalar_T( 1);

    const index_set_t ist_mn = index_set_t(-level);
    const index_set_t ist_pn = index_set_t(level);
    if (level == 1)
    {
      if (odd)
        return J * (*this)[ist_mn] + K  * (*this)[ist_pn];
      else
        return I * scalar(*this)   + JK * (*this)[ist_mn ^ ist_pn];
    }
    else
    {
      const framed_pair_t& pair_mn = this->divide(ist_mn);
      const framed_multi_t& quo_mn = pair_mn.first;
      const framed_multi_t& rem_mn = pair_mn.second;
      const framed_pair_t& pair_quo_mnpn = quo_mn.divide(ist_pn);
      const framed_multi_t& val_mnpn = pair_quo_mnpn.first;
      const framed_multi_t& val_mn   = pair_quo_mnpn.second;
      const framed_pair_t& pair_rem_mnpn = rem_mn.divide(ist_pn);
      const framed_multi_t& val_pn   = pair_rem_mnpn.first;
      const framed_multi_t& val_1    = pair_rem_mnpn.second;
      using matrix::kron;
      if (odd)
        return - kron(JK, val_1.fast   (level-1, 1))
               + kron(I,  val_mnpn.fast(level-1, 1))
               + kron(J,  val_mn.fast  (level-1, 0))
               + kron(K,  val_pn.fast  (level-1, 0));
      else
        return   kron(I,  val_1.fast   (level-1, 0))
               + kron(JK, val_mnpn.fast(level-1, 0))
               + kron(K,  val_mn.fast  (level-1, 1))
               - kron(J,  val_pn.fast  (level-1, 1));
    }
  }

  /// Use generalized FFT to construct a matrix_multi_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const typename framed_multi<Scalar_T,LO,HI>::matrix_multi_t
  framed_multi<Scalar_T,LO,HI>::
  fast_matrix_multi(const index_set_t frm) const
  {
    // Fold val
    multivector_t val = this->fold(frm);
    index_t p = frm.count_pos();
    index_t q = frm.count_neg();
    const index_t bott_offset = gen::offset_to_super[pos_mod(p - q, 8)];
    p += std::max(bott_offset,index_t(0));
    q -= std::min(bott_offset,index_t(0));
    if (p > HI)
      throw error_t("fast_matrix_multi(frm): HI is too low to represent this value");
    if (q > -LO)
      throw error_t("fast_matrix_multi(frm): LO is too high to represent this value");
    // Centre val
    while (p - q > 4)
      val.centre_pm4_qp4(p, q);
    while (p - q < -3)
      val.centre_pp4_qm4(p, q);
    if (p - q > 1)
      val.centre_qp1_pm1(p, q);
    const index_t level = (p + q)/2;

    // Do the fast transform
    const multivector_t& ev_val = val.even();
    const multivector_t& od_val = val.odd();
    return matrix_multi_t(ev_val.fast(level, 0) + od_val.fast(level, 1), frm);
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  fast_framed_multi() const
  { return *this; }

  /// Class name used in messages
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const std::string
  framed_multi<Scalar_T,LO,HI>::var_term_t::
  classname()
  { return "var_term"; }

  /// Default constructor
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::var_term_t::
  var_term() :
  var_pair_t(index_set_t(), Scalar_T(1))
  { }

  /// Construct a variable term from an index set and a scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::var_term_t::
  var_term(const index_set_t ist, const Scalar_T& crd) :
  var_pair_t(ist, crd)
  { }

  /// Product of variable term and term
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  typename framed_multi<Scalar_T,LO,HI>::var_term_t&
           framed_multi<Scalar_T,LO,HI>::var_term_t::
  operator*= (const term_t& rhs)
  {
    this->second *= rhs.second * this->first.sign_of_mult(rhs.first);
    this->first  ^= rhs.first;
    return *this;
  }

  /// Product of terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const std::pair<const index_set<LO,HI>, Scalar_T>
  operator* (const std::pair<const index_set<LO,HI>, Scalar_T>& lhs,
             const std::pair<const index_set<LO,HI>, Scalar_T>& rhs)
  {
    typedef std::pair<const index_set<LO,HI>, Scalar_T> term_t;
    return term_t(
           lhs.first ^ rhs.first,
           lhs.second * rhs.second *
           lhs.first.sign_of_mult(rhs.first));
  }
}
#endif  // _GLUCAT_FRAMED_MULTI_IMP_H
