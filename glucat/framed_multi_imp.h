#ifndef _GLUCAT_FRAMED_MULTI_IMP_H
#define _GLUCAT_FRAMED_MULTI_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    framed_multi_imp.h : Implement the coordinate map representation of a
    Clifford algebra element
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
               const index_set_t& frm, const bool prechecked)
  {
    if (!prechecked && (val.frame() | frm) != frm)
      throw error_t("multivector_t(val,frm): cannot initialize with value outside of frame");
    *this = val;
  }

  /// Construct a multivector from an index set and a scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const index_set_t& ist, const Scalar_T& crd)
  {
    if (crd != Scalar_T(0))
      this->insert(pair_t(ist, crd));
  }

  /// Construct a multivector, within a given frame, from an index set and a scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const index_set_t& ist, const Scalar_T& crd,
               const index_set_t& frm, const bool prechecked)
  {
    if (!prechecked && (ist | frm) != frm)
      throw error_t("multivector_t(ist,crd,frm): cannot initialize with value outside of frame");
    if (crd != Scalar_T(0))
      this->insert(pair_t(ist, crd));
  }

  /// Construct a multivector from a scalar (within a frame, if given)
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const Scalar_T& scr, const index_set_t& frm)
  {
    if (scr != 0)
      this->insert(pair_t(index_set_t(), scr));
  }

  /// Construct a multivector from an int (within a frame, if given)
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const int scr, const index_set_t& frm)
  {
    if (scr != 0)
      this->insert(pair_t(index_set_t(), Scalar_T(scr)));
  }

  /// Construct a multivector, within a given frame, from a given vector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const vector_t& vec,
               const index_set_t& frm, const bool prechecked)
  {
    if (!prechecked && index_t(vec.size()) != frm.count())
      throw error_t("multivector_t(vec,frm): cannot initialize with vector not matching frame");
    typename vector_t::const_iterator scvec = vec.begin();
    const index_t begin_index = frm.min();
    const index_t end_index = frm.max()+1;
    for (index_t idx = begin_index; idx != end_index; ++idx)
      if (frm[idx])
      {
        *this += pair_t(index_set_t(idx), *scvec);
        ++scvec;
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
  framed_multi(const std::string& str, const index_set_t& frm, const bool prechecked)
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
    const index_set_t frm = val.frame();
    const set_value_t end_set_value = 1 << frm.count();
    for (set_value_t stv = 0; stv != end_set_value; stv++)
    {
      const index_set_t ist(stv, frm, true);
      const Scalar_T crd =
        matrix::inner<Scalar_T>(
          basis_element<Scalar_T,LO,HI>(ist, frm), val.m_matrix);
      if (crd != Scalar_T(0))
        this->insert(pair_t(ist, crd));
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
    else
    {
      const_iterator scthis = this->begin();
      const_iterator scthat = rhs.begin();
      for (;
           (scthis != this->end()) || (scthat != rhs.end());
           scthis++, scthat++)
        if (*scthis != *scthat)
          return false;
      return (scthis == this->end()) && (scthat == rhs.end());
    }
  }

  // Test for equality of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  bool
  framed_multi<Scalar_T,LO,HI>::
  operator==  (const Scalar_T& scr) const
  { return frame().count() == 0 && scalar(*this) == scr; }


  /// Geometric sum of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator+= (const Scalar_T& scr)
  {
    *this += pair_t(index_set_t(),scr);
    return *this;
  }

  /// Geometric sum
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator+= (const multivector_t& rhs)
  { // simply add terms
    for(const_iterator 
        scan = rhs.begin(); scan != rhs.end();  ++scan)
      *this += *scan;
    return *this;
  }

  /// Geometric difference
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator-= (const multivector_t& rhs)
  {
    for(const_iterator 
        scan = rhs.begin(); scan != rhs.end();  ++scan)
      *this += pair_t(scan->first, -(scan->second));
    return *this;
  }

  /// Unary -
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  operator- () const
  { return *this * Scalar_T(-1.0); }

  /// Product of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator*= (const Scalar_T& scr)
  { // multiply coordinates of all terms by scalar
    if (scr == Scalar_T(0))
      this->clear();
    else
      for (iterator scan = this->begin(); scan != this->end(); ++scan)
        scan->second *= scr;
    return *this;
  }

  /// Geometric product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator*= (const multivector_t& rhs)
  {
    const index_set_t our_frame = frame() | rhs.frame();
    const index_t frm_count = our_frame.count();
    const set_value_t array_size = 1 << frm_count;
    if (this->size()*rhs.size() <= array_size)
    { // If we have a sparse multiply, store the result directly
      multivector_t  result;
      for (const_iterator scthis = this->begin(); scthis != this->end(); scthis++)
        for (const_iterator scthat = rhs.begin(); scthat != rhs.end(); scthat++)
          result += (*scthis) * (*scthat);
      return *this = result;
    }
    if (frm_count >= Tune_P::mult_matrix_threshold)
      // Past a certain threshold, the matrix algorithm is fastest
      return *this = matrix_multi_t(*this, our_frame, true) *
                     matrix_multi_t(rhs, our_frame, true);
    else
    { // Fastest dense algorithm in low dimensions stores result in array
      std::vector< Scalar_T > result_array( array_size, Scalar_T(0) );
      for (const_iterator scthis = this->begin(); scthis != this->end(); ++scthis)
        for (const_iterator scthat = rhs.begin(); scthat != rhs.end(); ++scthat)
        {
          const pair_t result_pair = (*scthis) * (*scthat);
          const set_value_t stv = result_pair.first.value_of_fold(our_frame);
          result_array[stv] += result_pair.second;
        }
      *this = Scalar_T(0);
      for (set_value_t stv = 0; stv < array_size; ++stv)
        if (result_array[stv] != Scalar_T(0))
        {
          const index_set_t result_ist(stv, our_frame);
          this->insert(pair_t(result_ist, result_array[stv]));
        }
      return *this;
    }
  }

  /// Left contraction
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator%= (const multivector_t& rhs)
  {
    // Reference: Leo Dorst, "Honing geometric algebra for its use in the computer sciences",
    // http://carol.wins.uva.nl/~leo/clifford/index.html
    // also in Geometric Computing with Clifford Algebras, ed. G. Sommer,
    // Springer 2001, Chapter 6, pp. 127-152
    multivector_t result;
    for (const_iterator scthis = this->begin(); scthis != this->end(); scthis++)
      for (const_iterator scthat = rhs.begin(); scthat != rhs.end(); scthat++)
      {
        const pair_t term = (*scthis) * (*scthat);
        if ((scthis->first | scthat->first) == scthat->first)
          result += term;
      }
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
    for (const_iterator scthis = this->begin(); scthis != this->end(); scthis++)
      for (const_iterator scthat = rhs.begin(); scthat != rhs.end(); scthat++)
      {
        const index_t grade = std::abs(scthis->first.count() - scthat->first.count());
        const pair_t& term = (*scthis) * (*scthat);
        if (term.first.count() == grade)
          result += term;
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
    multivector_t  result;
    for (const_iterator scthis = this->begin(); scthis != this->end(); scthis++)
      for (const_iterator scthat = rhs.begin(); scthat != rhs.end(); scthat++)
      {
        const index_t grade = scthis->first.count() + scthat->first.count();
        if (grade <= HI-LO)
        {
          const pair_t& term = (*scthis) * (*scthat);
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
    for (iterator scan = this->begin(); scan != this->end(); ++scan)
      scan->second /= scr;
    return *this;
  }

  /// Geometric quotient
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator/= (const multivector_t& rhs)
  {
    index_set_t our_frame = frame() | rhs.frame();
    matrix_multi_t result(*this, our_frame, true);
    return *this = result /= matrix_multi_t(rhs, our_frame, true);
  }

  /// Clifford multiplicative inverse
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  inv() const
  {
    matrix_multi_t result(1, frame());
    return result /= matrix_multi_t(*this);
  }

  /// Subscripting: map from index set to scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  Scalar_T
  framed_multi<Scalar_T,LO,HI>::
  operator[] (const index_set_t& ist) const
  {
      const_iterator scan = this->find(ist);
      if (scan == this->end())
        return 0;
      else
        return scan->second;
  }

  /// Main involution, each {i} is replaced by -{i} in each term
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  involute() const
  {
    multivector_t result(*this);
    for (iterator scan = result.begin(); scan != result.end(); ++scan)
    { // for a k-vector u, involute(u) == (-1)^k * u
      if ((scan->first.count() % 2) == 1)
        scan->second *= Scalar_T(-1.0);
    }
    return result;
  }

  /// Reversion, order of {i} is reversed in each term
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  reverse() const
  {
    multivector_t result(*this);
    for (iterator scan = result.begin(); scan != result.end(); ++scan)
      // For a k-vector u, reverse(u) = { -u, k == 2,3 (mod 4)
      //                                {  u, k == 0,1 (mod 4)
      switch (scan->first.count() % 4)
      {
      case 2:
      case 3:
        scan->second *= Scalar_T(-1.0);
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
    multivector_t result(*this);
    for (iterator scan = result.begin(); scan != result.end(); ++scan)
      // For a k-vector u, conj(u) = { -u, k == 1,2 (mod 4)
      //                             {  u, k == 0,3 (mod 4)
      switch (scan->first.count() % 4)
      {
      case 1:
      case 2:
        scan->second *= Scalar_T(-1.0);
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
    for (const_iterator scan = this->begin(); scan != this->end(); ++scan)
    {
      int sign = 1;
      for(index_t i = LO; i < 0; ++i)
        sign *= (scan->first.test(i)) ? -1 : 1;
      result += sign * (scan->second) * (scan->second);
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
    for (const_iterator scan = this->begin(); scan != this->end(); ++scan)
    {
      const Scalar_T n2 = ublas::type_traits<Scalar_T>::norm_2(scan->second);
      result +=  n2 * n2;
    }
    return result;
  }

  /// *this to the m
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  pow(int m) const
  {
    multivector_t a;
    if (m < 0)
    {
      m = -m;
      a = (*this).inv();
    }
    else
      a = *this;
    multivector_t result = 1;
    for (; m != 0; m >>= 1, a *= a)
      if (m & 1)
        result *= a;
    return result;
  }

  /// Outer product power
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const
  framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  outer_pow(int m) const
  {
    if (m < 0)
      throw error_t("outer_pow(int): negative exponent");
    multivector_t result = 1;
    multivector_t a = *this;
    for (; m != 0; m >>= 1, a ^= a)
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
      return 0;
    else
    {
      multivector_t result;
      const_iterator scan = this->begin();
      while (scan != this->end() && scan->first.count() < grade)
        ++scan;
      for (; scan != this->end() && scan->first.count() == grade; ++scan)
        result += *scan;
      return result;
    }
  }

  /// Even part
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  even() const
  { // even part of x, sum of the pure(count) with even count
    multivector_t result;
    for(index_t even_grade = 0; even_grade <= HI-LO; even_grade += 2)
      result += (*this)(even_grade);
    return result;
  }

  /// Vector part of multivector, as a vector_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const typename framed_multi<Scalar_T,LO,HI>::vector_t
  framed_multi<Scalar_T,LO,HI>::
  vector_part() const
  {
    vector_t result;
    index_set_t frm = frame();
    const_iterator scan = this->begin();
    while (scan != this->end() && scan->first.count() < 1)
      ++scan;
    index_t idx = frm.min();
    const index_t end_index = frm.max()+1;
    for (; scan != this->end() && scan->first.count() == 1; ++scan, ++idx)
    {
      // Within framed_multi_t, terms with grade 1 are in order of increasing index
      // Frame may contain indices which do not correspond to a grade 1 term but
      // frame cannot omit any index corresponding to a grade 1 term
      const index_t scidx = scan->first.min();
      for(; idx != scidx; ++idx)
        if (frm[idx])
          result.push_back(Scalar_T(0));
      result.push_back(scan->second);
    }
    for(; idx != end_index; ++idx)
      if (frm[idx])
        result.push_back(Scalar_T(0));
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

  /// Write multivector to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::ostream&
  operator<< (std::ostream& os, const framed_multi<Scalar_T,LO,HI>& val)
  {
    typedef framed_multi<Scalar_T,LO,HI>  multivector_t;
    if(val.empty())
      os << 0;
    else
    {
      typename multivector_t::const_iterator scan = val.begin();
      os << *scan;
      for(++scan; scan != val.end(); ++scan)
      {
        Scalar_T scr = scan->second;
        if ((scr != std::real(scr)) || (std::real(scr) >= 0.0))
          os << '+';
        os << *scan;
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
    else if (term.second == Scalar_T(-1.0))
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
    framed_multi<Scalar_T,LO,HI>  local;
    char c = 0;
    while (s)
    {  // Default coordinate value is Scalar_T(0)
      Scalar_T     coordinate = Scalar_T(0);
      // Default index set is empty
      index_set<LO,HI> ist;
      ist.reset();
      c = s.peek();
      // Look for a leading '-'
      const bool negative = (c == '-');
      if (negative)
      { // consume the '-'
        s >> c;
        if (s)
          c = s.peek();
      }
      // Look for a leading '{'
      if (c == '{' && s)
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
        if (c == '{' && s)
          // Try to read index set.
          s >> ist;
      }
      if (!s.bad() && coordinate != Scalar_T(0))
      { // Insert whatever we have
        coordinate = negative ? -coordinate : coordinate;
        local += framed_multi<Scalar_T,LO,HI> (ist, coordinate);
      }
      if(s)
      {
        c = s.peek();
        if (c == '+')
          // consume the '+'
          s >> c;
      }
    }
    if (!s.bad())
      // If s.bad() then we have a corrupt input
      // otherwise we are fine and can copy the resulting framed_multi
      val = local;
    return s;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >   inline
  framed_multi<Scalar_T,LO,HI> &
  framed_multi<Scalar_T,LO,HI>::
  operator+= (const pair_t& term)
  { // Insert a term into a multivector_t, add terms with same index set.
    // Do not insert terms with 0 coordinate
    if (term.second != Scalar_T(0))
    {
      multivector_t::iterator scan = this->find(term.first);
      if (scan == this->end())
        this->insert(term);
      else if (scan->second + term.second == Scalar_T(0))
        // Erase term if resulting coordinate is 0
        this->erase(scan);
      else
        scan->second += term.second;
    }
    return *this;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  const index_set<LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  frame() const
  {
    index_set_t result;
    for (const_iterator scan = this->begin(); scan != this->end(); ++scan)
      result |= scan->first;
    return result;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  Scalar_T
  framed_multi<Scalar_T,LO,HI>::
  max_abs() const
  {
    Scalar_T result = Scalar_T(0);
    for (const_iterator scan = this->begin(); scan != this->end(); ++scan)
    {
      const Scalar_T scan_abs = std::abs(scan->second);
      if (scan_abs > result)
        result = scan_abs;
    }
    return result;
  }

  /// Check if a multivector contains any IEEE NaN values
  template< typename Scalar_T, const index_t LO, const index_t HI >
  bool
  framed_multi<Scalar_T,LO,HI>::
  isnan() const
  { // The distinguishing feature is that NaN != NaN
    for (const_iterator scan = this->begin(); scan != this->end(); ++scan)
        if (*scan != *scan)
          return true;
    return false;
  }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  framed_multi<Scalar_T,LO,HI>::
  truncated(const Scalar_T& limit) const
  {
    Scalar_T top = max_abs();
    multivector_t result;
    if (top != Scalar_T(0))
      for (const_iterator scan = this->begin(); scan != this->end(); ++scan)
        if (std::abs(scan->second / top) > limit)
          result.insert(pair_t(scan->first, scan->second));
    return result;
  }

  /// Product of terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  const std::pair<const index_set<LO,HI>, Scalar_T>
  operator* (const std::pair<const index_set<LO,HI>, Scalar_T>& lhs,
             const std::pair<const index_set<LO,HI>, Scalar_T>& rhs)
  {
    typedef std::pair<const index_set<LO,HI>, Scalar_T> pair_t;
    return (pair_t(
            lhs.first ^ rhs.first,
            lhs.second * rhs.second *
            lhs.first.sign_of_mult(rhs.first)));
  }
}
#endif  // _GLUCAT_FRAMED_MULTI_IMP_H
