#ifndef _GLUCAT_FRAMED_MULTI_IMP_H
#define _GLUCAT_FRAMED_MULTI_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    framed_multi_imp.h : Implement the coordinate map representation of a
    Clifford algebra element
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

#include "glucat/framed_multi.h"

#include "glucat/scalar.h"
#include "glucat/random.h"
#include "glucat/generation.h"
#include "glucat/matrix.h"

#if defined(_GLUCAT_USE_BOOST_POOL_ALLOC)
// Use the Boost pool allocator
#include <boost/pool/pool_alloc.hpp>
#endif

#include <sstream>
#include <fstream>

namespace glucat
{
  /// Class name used in messages
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  classname() -> const std::string
  { return "framed_multi"; }

#if defined(_GLUCAT_MAP_IS_HASH)
#define _GLUCAT_HASH_N(x) (x)
#define _GLUCAT_HASH_SIZE_T(x) (typename multivector_t::hash_size_t)(x)
#else
#define _GLUCAT_HASH_N(x)
#define _GLUCAT_HASH_SIZE_T(x)
#endif

  /// Default constructor
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi()
  : map_t(_GLUCAT_HASH_N(0))
  { }

  /// Private constructor using hash_size
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const hash_size_t& hash_size)
  : map_t(_GLUCAT_HASH_N(hash_size()))
  { }

  /// Construct a multivector from a multivector with a different scalar type
  template< typename Scalar_T, const index_t LO, const index_t HI >
  template< typename Other_Scalar_T >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const framed_multi<Other_Scalar_T,LO,HI>& val)
  : map_t(_GLUCAT_HASH_N(val.size()))
  {
    using other_multivector_t = framed_multi<Other_Scalar_T, LO, HI>;
    using other_const_iterator = typename other_multivector_t::const_iterator;
    const other_const_iterator val_begin = val.begin();
    const other_const_iterator val_end   = val.end();
    for (other_const_iterator val_it = val_begin; val_it != val_end; ++val_it)
      this->insert(term_t(val_it->first, numeric_traits<Scalar_T>::to_scalar_t(val_it->second)));
  }

  /// Construct a multivector, within a given frame, from a given multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  template< typename Other_Scalar_T >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const framed_multi<Other_Scalar_T,LO,HI>& val,
               const index_set_t frm, const bool prechecked)
  : map_t(_GLUCAT_HASH_N(val.size()))
  {
    using other_multivector_t = framed_multi<Other_Scalar_T, LO, HI>;
    using other_const_iterator = typename other_multivector_t::const_iterator;
    const other_const_iterator val_begin = val.begin();
    const other_const_iterator val_end   = val.end();
    for (other_const_iterator val_it = val_begin; val_it != val_end; ++val_it)
      this->insert(term_t(val_it->first, static_cast<Scalar_T>(val_it->second)));
  }

  /// Construct a multivector, within a given frame, from a given multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const multivector_t& val,
               const index_set_t frm, const bool prechecked)
  : map_t(val)
  { }

  /// Construct a multivector from an index set and a scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const index_set_t ist, const Scalar_T& crd)
  : map_t(_GLUCAT_HASH_N(1))
  {
    if (crd != Scalar_T(0))
      this->insert(term_t(ist, crd));
  }

  /// Construct a multivector, within a given frame, from an index set and a scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const index_set_t ist, const Scalar_T& crd,
               const index_set_t frm, const bool prechecked)
  : map_t(_GLUCAT_HASH_N(1))
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
  : map_t(_GLUCAT_HASH_N(1))
  {
    if (scr != Scalar_T(0))
      this->insert(term_t(index_set_t(), scr));
  }

  /// Construct a multivector from an int (within a frame, if given)
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const int scr, const index_set_t frm)
  : map_t(_GLUCAT_HASH_N(1))
  {
    if (scr != Scalar_T(0))
      this->insert(term_t(index_set_t(), Scalar_T(scr)));
  }

  /// Construct a multivector, within a given frame, from a given vector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const vector_t& vec,
               const index_set_t frm, const bool prechecked)
  : map_t(_GLUCAT_HASH_N(vec.size()))
  {
    if (!prechecked && index_t(vec.size()) != frm.count())
      throw error_t("multivector_t(vec,frm): cannot initialize with vector not matching frame");
    auto vec_it = vec.begin();
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
  : map_t(_GLUCAT_HASH_N(0))
  {
    std::istringstream ss(str);
    ss >> *this;
    if (!ss)
      throw error_t("multivector_t(str): could not parse string");
    // Peek to see if the end of the string has been reached.
    ss.peek();
    if (!ss.eof())
      throw error_t("multivector_t(str): could not parse entire string");
  }

  /// Construct a multivector, within a given frame, from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
  template< typename Scalar_T, const index_t LO, const index_t HI >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const std::string& str, const index_set_t frm, const bool prechecked)
  : map_t(_GLUCAT_HASH_N(0))
  {
    if (prechecked)
      *this = multivector_t(str);
    else
      *this = multivector_t(multivector_t(str), frm, false);
  }

  /// Construct a multivector from a matrix_multi_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  template< typename Other_Scalar_T >
  framed_multi<Scalar_T,LO,HI>::
  framed_multi(const matrix_multi<Other_Scalar_T,LO,HI>& val)
  : map_t(_GLUCAT_HASH_N(1))
  {
    if (val == Other_Scalar_T(0))
      return;

    using matrix_index_t = typename matrix_multi_t::matrix_index_t;
    const matrix_index_t dim = val.m_matrix.size1();
    using traits_t = numeric_traits<Scalar_T>;
    if (dim == 1)
    {
      this->insert(term_t(index_set_t(), traits_t::to_scalar_t(val.m_matrix(0, 0))));
      return;
    }
    if (dim >= Tune_P::inv_fast_dim_threshold)
      try
      {
        *this = val.template fast_framed_multi<Scalar_T>();
        return;
      }
      catch (const glucat_error& e)
      { }

    const index_set_t frm = val.frame();
    const set_value_t algebra_dim = 1 << frm.count();
    const Scalar_T val_norm = traits_t::to_scalar_t(val.norm());
    if (traits_t::isNaN_or_isInf(val_norm))
    {
      *this = traits_t::NaN();
      return;
    }
    const Scalar_T eps = std::numeric_limits<Scalar_T>::epsilon();
    const Scalar_T tol = traits_t::abs(val_norm * eps * eps);
#if defined(_GLUCAT_MAP_IS_HASH)
    const size_t max_size = std::min<size_t>(algebra_dim, matrix::nnz(val.m_matrix));
    *this = multivector_t(_GLUCAT_HASH_SIZE_T(max_size));
#endif
    for (set_value_t
        stv = 0;
        stv != algebra_dim;
        stv++)
    {
      const index_set_t ist = index_set_t(stv, frm, true);
      const Scalar_T crd =
        traits_t::to_scalar_t(matrix::inner<Other_Scalar_T>(val.basis_element(ist), val.m_matrix));
      const Scalar_T abs_crd = traits_t::abs(crd);
      if ((abs_crd * abs_crd) > tol)
        this->insert(term_t(ist, crd));
    }
  }

  /// Test for equality of multivectors
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator==  (const multivector_t& rhs) const -> bool
  {
    if (this->size() != rhs.size())
      return false;
    const const_iterator this_begin = this->begin();
    const const_iterator this_end   = this->end();
    const const_iterator rhs_end    = rhs.end();
#if defined(_GLUCAT_MAP_IS_ORDERED)
    const_iterator this_it = this_begin;
    const_iterator rhs_it = rhs.begin();
    for (;
        (this_it != this_end) && (rhs_it != rhs_end);
        this_it++, rhs_it++)
      if (*this_it != *rhs_it)
        return false;
    return (this_it == this_end) && (rhs_it == rhs_end);
#else
    for (const_iterator
        this_it = this_begin;
        this_it != this_end;
        this_it++)
    {
      const const_iterator& rhs_it = rhs.find(this_it->first);
      if (rhs_it == rhs_end || rhs_it->second != this_it->second)
        return false;
    }
    return true;
#endif
  }

  /// Test for equality of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator==  (const Scalar_T& scr) const -> bool
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
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator+= (const Scalar_T& scr) -> multivector_t&
  {
    *this += term_t(index_set_t(), scr);
    return *this;
  }

  /// Geometric sum
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator+= (const multivector_t& rhs) -> multivector_t&
  { // simply add terms
    for (auto
        rhs_it = rhs.begin();
        rhs_it != rhs.end();
        ++rhs_it)
      *this += *rhs_it;
    return *this;
  }

  /// Geometric difference of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator-= (const Scalar_T& scr) -> multivector_t&
  {
    *this += term_t(index_set_t(), -scr);
    return *this;
  }

  /// Geometric difference
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator-= (const multivector_t& rhs) -> multivector_t&
  {
    for (auto
        rhs_it = rhs.begin();
        rhs_it != rhs.end();
        ++rhs_it)
      *this += term_t(rhs_it->first, -(rhs_it->second));
    return *this;
  }

  /// Unary -
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator- () const -> const multivector_t
  { return *this * Scalar_T(-1); }

  /// Product of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator*= (const Scalar_T& scr) -> multivector_t&
  { // multiply coordinates of all terms by scalar
    using traits_t = numeric_traits<Scalar_T>;

    if (traits_t::isNaN_or_isInf(scr))
      return *this = traits_t::NaN();
    if (scr == Scalar_T(0))
      if (this->isnan())
        *this = traits_t::NaN();
      else
        this->clear();
    else
      for (auto
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
        this_it->second *= scr;
    return *this;
  }

  /// Geometric product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator* (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> const framed_multi<Scalar_T,LO,HI>
  {
    using multivector_t = framed_multi<Scalar_T, LO, HI>;
    using traits_t = numeric_traits<Scalar_T>;
    using index_set_t = typename multivector_t::index_set_t;
    using term_t = typename multivector_t::term_t;
    using map_t = typename multivector_t::map_t;
    using const_iterator = typename map_t::const_iterator;

    if (lhs.isnan() || rhs.isnan())
      return traits_t::NaN();

    const double lhs_size = lhs.size();
    const double rhs_size = rhs.size();
    const index_set_t our_frame = lhs.frame() | rhs.frame();
    const index_t frm_count = our_frame.count();
    const set_value_t algebra_dim = 1 << frm_count;
    const bool direct_mult = lhs_size * rhs_size <= double(algebra_dim)
#if defined(_GLUCAT_MAP_IS_HASH)
                           || frm_count < Tune_P::mult_matrix_threshold
#endif
                           ;
    if (direct_mult)
    { // If we have a sparse multiply, store the result directly
      multivector_t result =
        multivector_t(_GLUCAT_HASH_SIZE_T(size_t(std::min(lhs_size * rhs_size, double(algebra_dim)))));
      const const_iterator lhs_begin = lhs.begin();
      const const_iterator lhs_end   = lhs.end();
      const const_iterator rhs_begin = rhs.begin();
      const const_iterator rhs_end   = rhs.end();

      for (const_iterator
          lhs_it = lhs_begin;
          lhs_it != lhs_end;
          ++lhs_it)
      {
        const term_t& lhs_term = *lhs_it;
        for (const_iterator
            rhs_it = rhs_begin;
            rhs_it != rhs_end;
            ++rhs_it)
          result += lhs_term * *rhs_it;
      }
      return result;
    }
#if !defined(_GLUCAT_MAP_IS_HASH)
    else if (frm_count < Tune_P::mult_matrix_threshold)
    { // Fastest dense algorithm in low dimensions stores result in array
      typedef std::vector<Scalar_T> array_t;
      array_t result_array(algebra_dim, Scalar_T(0));

      const const_iterator lhs_begin = lhs.begin();
      const const_iterator lhs_end   = lhs.end();
      const const_iterator rhs_begin = rhs.begin();
      const const_iterator rhs_end   = rhs.end();

      for (const_iterator
          lhs_it = lhs_begin;
          lhs_it != lhs_end;
          ++lhs_it)
      {
        const term_t& lhs_term = *lhs_it;
        for (const_iterator
            rhs_it = rhs_begin;
            rhs_it != rhs_end;
            ++rhs_it)
        {
          const term_t& term = lhs_term * *rhs_it;
          const set_value_t stv = term.first.value_of_fold(our_frame);
          result_array[stv] += term.second;
        }
      }
      multivector_t result;
      for (set_value_t
          stv = 0;
          stv != algebra_dim;
          ++stv)
        if (result_array[stv] != Scalar_T(0))
          result.insert(term_t(index_set_t(stv, our_frame, true), result_array[stv]));
      return result;
    }
#endif
    else
    { // Past a certain threshold, the matrix algorithm is fastest
      using matrix_multi_t = typename multivector_t::matrix_multi_t;
      return matrix_multi_t(lhs, our_frame, true) *
             matrix_multi_t(rhs, our_frame, true);
    }
  }

  /// Geometric product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator*= (const multivector_t& rhs) -> multivector_t&
  { return *this = *this * rhs; }

  /// Outer product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator^ (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> const framed_multi<Scalar_T,LO,HI>
  { // Arvind Raja's original reference:
    // "old clical, outerproduct(p,q:pterm):pterm in file compmod.pas"

    if (lhs.empty() || rhs.empty())
      return Scalar_T(0);

    using multivector_t = framed_multi<Scalar_T, LO, HI>;
    using index_set_t = typename multivector_t::index_set_t;
    using term_t = typename multivector_t::term_t;
    using map_t = typename multivector_t::map_t;
    using const_iterator = typename map_t::const_iterator;

    multivector_t result;
    const index_set_t empty_set = index_set_t();

    const const_iterator lhs_begin = lhs.begin();
    const const_iterator lhs_end   = lhs.end();
    const const_iterator rhs_begin = rhs.begin();
    const const_iterator rhs_end   = rhs.end();
    const double lhs_size = lhs.size();
    const double rhs_size = rhs.size();

    if (lhs_size * rhs_size > double(Tune_P::products_size_threshold))
    {
      const index_set_t lhs_frame = lhs.frame();
      const index_set_t rhs_frame = rhs.frame();

      const index_set_t our_frame = lhs_frame | rhs_frame;
      const set_value_t algebra_dim = 1 << our_frame.count();
      multivector_t result =
        multivector_t(_GLUCAT_HASH_SIZE_T(size_t(std::min(lhs_size * rhs_size, double(algebra_dim)))));
      for (set_value_t
          result_stv = 0;
          result_stv != algebra_dim;
          ++result_stv)
      {
        const index_set_t result_ist = index_set_t(result_stv, our_frame, true);
        const index_set_t lhs_result_frame = lhs_frame & result_ist;
        const set_value_t lhs_result_dim = 1 << lhs_result_frame.count();
        auto result_crd = Scalar_T(0);
        for (set_value_t
            lhs_stv = 0;
            lhs_stv != lhs_result_dim;
            ++lhs_stv)
        {
          const index_set_t lhs_ist = index_set_t(lhs_stv, lhs_result_frame, true);
          const index_set_t rhs_ist = result_ist ^ lhs_ist;
          if ((rhs_ist | rhs_frame) == rhs_frame)
          {
            const const_iterator lhs_it = lhs.find(lhs_ist);
            if (lhs_it != lhs_end)
            {
              const const_iterator rhs_it = rhs.find(rhs_ist);
              if (rhs_it != rhs_end)
                result_crd += crd_of_mult(*lhs_it, *rhs_it);
            }
          }
        }
        if (result_crd != Scalar_T(0))
          result.insert(term_t(result_ist, result_crd));
      }
      return result;
    }
    else
    {
      multivector_t result;
      for (const_iterator
          rhs_it = rhs_begin;
          rhs_it != rhs_end;
          ++rhs_it)
      {
        const term_t& rhs_term = *rhs_it;
        const index_set_t rhs_ist = rhs_term.first;
        for (const_iterator
            lhs_it = lhs_begin;
            lhs_it != lhs_end;
            ++lhs_it)
        {
          const term_t& lhs_term = *lhs_it;
          const index_set_t lhs_ist = lhs_term.first;
          if ((lhs_ist & rhs_ist) == empty_set)
            result += lhs_term * rhs_term;
        }
      }
      return result;
    }
  }

  /// Outer product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator^= (const multivector_t& rhs) -> multivector_t&
  { return *this = *this ^ rhs; }

  /// Inner product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator& (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> const framed_multi<Scalar_T,LO,HI>
  { // Arvind Raja's original reference:
    // "old clical, innerproduct(p,q:pterm):pterm in file compmod.pas"

    if (lhs.empty() || rhs.empty())
      return Scalar_T(0);

    using multivector_t = framed_multi<Scalar_T, LO, HI>;
    using index_set_t = typename multivector_t::index_set_t;
    using term_t = typename multivector_t::term_t;
    using map_t = typename multivector_t::map_t;
    using const_iterator = typename map_t::const_iterator;

    const const_iterator lhs_end   = lhs.end();
    const const_iterator rhs_end   = rhs.end();
    const double lhs_size = lhs.size();
    const double rhs_size = rhs.size();

    if (lhs_size * rhs_size > double(Tune_P::products_size_threshold))
    {
      const index_set_t lhs_frame = lhs.frame();
      const index_set_t rhs_frame = rhs.frame();

      const index_set_t our_frame = lhs_frame | rhs_frame;
      const set_value_t algebra_dim = 1 << our_frame.count();
      multivector_t result =
        multivector_t(_GLUCAT_HASH_SIZE_T(size_t(std::min(lhs_size * rhs_size, double(algebra_dim)))));
      for (set_value_t
          result_stv = 0;
          result_stv != algebra_dim;
          ++result_stv)
      {
        const index_set_t result_ist = index_set_t(result_stv, our_frame, true);
        const index_set_t comp_frame = our_frame & ~result_ist;
        const set_value_t comp_dim = 1 << comp_frame.count();
        auto result_crd = Scalar_T(0);
        for (set_value_t
            comp_stv = 1;
            comp_stv != comp_dim;
            ++comp_stv)
        {
          const index_set_t comp_ist = index_set_t(comp_stv, comp_frame, true);
          const index_set_t our_ist = result_ist ^ comp_ist;
          if ((our_ist | lhs_frame) == lhs_frame)
          {
            const const_iterator lhs_it = lhs.find(our_ist);
            if (lhs_it != lhs_end)
            {
              const const_iterator rhs_it = rhs.find(comp_ist);
              if (rhs_it != rhs_end)
                result_crd += crd_of_mult(*lhs_it, *rhs_it);
            }
          }
          if (result_stv != 0)
          {
            if ((our_ist | rhs_frame) == rhs_frame)
            {
              const const_iterator rhs_it = rhs.find(our_ist);
              if (rhs_it != rhs_end)
              {
                const const_iterator lhs_it = lhs.find(comp_ist);
                if (lhs_it != lhs_end)
                  result_crd += crd_of_mult(*lhs_it, *rhs_it);
              }
            }
          }
        }
        if (result_crd != Scalar_T(0))
          result.insert(term_t(result_ist, result_crd));
      }
      return result;
    }
    else
    {
      const index_set_t empty_set = index_set_t();

      const const_iterator lhs_begin = lhs.begin();
      const const_iterator rhs_begin = rhs.begin();

      multivector_t result;
      for (const_iterator
          lhs_it = lhs_begin;
          lhs_it != lhs_end;
          ++lhs_it)
      {
        const term_t& lhs_term = *lhs_it;
        const index_set_t lhs_ist = lhs_term.first;
        if (lhs_ist != empty_set)
          for (const_iterator
              rhs_it = rhs_begin;
              rhs_it != rhs_end;
              ++rhs_it)
          {
            const term_t& rhs_term = *rhs_it;
            const index_set_t rhs_ist = rhs_term.first;
            if (rhs_ist != empty_set)
            {
              const index_set_t our_ist = lhs_ist | rhs_ist;
              if ((lhs_ist == our_ist) || (rhs_ist == our_ist))
                result += lhs_term * rhs_term;
            }
          }
      }
      return result;
    }
  }

  /// Inner product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator&= (const multivector_t& rhs) -> multivector_t&
  { return *this = *this & rhs; }

  /// Left contraction
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator% (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> const framed_multi<Scalar_T,LO,HI>
  {
    // Reference: Leo Dorst, "Honing geometric algebra for its use in the computer sciences",
    // in Geometric Computing with Clifford Algebras, ed. G. Sommer,
    // Springer 2001, Chapter 6, pp. 127-152.
    // http://staff.science.uva.nl/~leo/clifford/index.html

    if (lhs.empty() || rhs.empty())
      return Scalar_T(0);

    using multivector_t = framed_multi<Scalar_T, LO, HI>;
    using index_set_t = typename multivector_t::index_set_t;
    using term_t = typename multivector_t::term_t;
    using map_t = typename multivector_t::map_t;
    using const_iterator = typename map_t::const_iterator;

#if defined(_GLUCAT_MAP_IS_ORDERED)
    // Both lhs and rhs are sorted by increasing grade, then lexicographically,
    // and a "larger" index set cannot be a subset of a "smaller" one.

    const const_iterator lhs_begin = lhs.begin();

    typedef typename map_t::const_reverse_iterator const_reverse_iterator;
    const const_reverse_iterator rhs_rbegin = rhs.rbegin();
    const const_reverse_iterator rhs_rlower_bound =
          static_cast<const_reverse_iterator>(rhs.lower_bound(lhs_begin->first));

    multivector_t result;

    for (const_reverse_iterator
        rhs_it = rhs_rbegin;
        rhs_it != rhs_rlower_bound;
        ++rhs_it)
    {
      const term_t& rhs_term = *rhs_it;
      const index_set_t rhs_ist = rhs_term.first;
      const const_iterator lhs_upper_bound = lhs.upper_bound(rhs_ist);
      for (const_iterator
          lhs_it = lhs_begin;
          lhs_it != lhs_upper_bound;
          ++lhs_it)
      {
        const term_t& lhs_term = *lhs_it;
        const index_set_t lhs_ist = lhs_term.first;
        if ((lhs_ist | rhs_ist) == rhs_ist)
          result += lhs_term * rhs_term;
      }
    }
    return result;
#else
    const const_iterator lhs_end   = lhs.end();
    const const_iterator rhs_end   = rhs.end();
    const double lhs_size = lhs.size();
    const double rhs_size = rhs.size();

    if (lhs_size * rhs_size > double(Tune_P::products_size_threshold))
    {
      const index_set_t lhs_frame = lhs.frame();
      const index_set_t rhs_frame = rhs.frame();

      const index_set_t our_frame = lhs_frame | rhs_frame;
      const set_value_t algebra_dim = 1 << our_frame.count();
      multivector_t result =
        multivector_t(_GLUCAT_HASH_SIZE_T(size_t(std::min(lhs_size * rhs_size, double(algebra_dim)))));
      for (set_value_t
          result_stv = 0;
          result_stv != algebra_dim;
          ++result_stv)
      {
        const index_set_t result_ist = index_set_t(result_stv, our_frame, true);
        const index_set_t comp_frame = lhs_frame & ~result_ist;
        const set_value_t comp_dim = 1 << comp_frame.count();
        auto result_crd = Scalar_T(0);
        for (set_value_t
            comp_stv = 0;
            comp_stv != comp_dim;
            ++comp_stv)
        {
          const index_set_t comp_ist = index_set_t(comp_stv, comp_frame, true);
          const index_set_t rhs_ist = result_ist ^ comp_ist;
          if ((rhs_ist | rhs_frame) == rhs_frame)
          {
            const const_iterator rhs_it = rhs.find(rhs_ist);
            if (rhs_it != rhs_end)
            {
              const const_iterator lhs_it = lhs.find(comp_ist);
              if (lhs_it != lhs_end)
                result_crd += crd_of_mult(*lhs_it, *rhs_it);
            }
          }
        }
        if (result_crd != Scalar_T(0))
          result.insert(term_t(result_ist, result_crd));
      }
      return result;
    }
    else
    {
      const const_iterator rhs_begin = rhs.begin();
      const const_iterator lhs_begin = lhs.begin();

      multivector_t result;
      for (const_iterator
          rhs_it = rhs_begin;
          rhs_it != rhs_end;
          ++rhs_it)
      {
        const term_t& rhs_term = *rhs_it;
        const index_set_t rhs_ist = rhs_term.first;
        for (const_iterator
            lhs_it = lhs_begin;
            lhs_it != lhs_end;
            ++lhs_it)
        {
          const term_t& lhs_term = *lhs_it;
          const index_set_t lhs_ist = lhs_term.first;
          if ((lhs_ist | rhs_ist) == rhs_ist)
            result += lhs_term * rhs_term;
        }
      }
      return result;
    }
#endif
  }

  /// Left contraction
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator%= (const multivector_t& rhs) -> multivector_t&
  { return *this = *this % rhs; }

  /// Hestenes scalar product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  star(const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> Scalar_T
  {
    using multivector_t = framed_multi<Scalar_T, LO, HI>;
    using map_t = typename multivector_t::map_t;
    using const_iterator = typename map_t::const_iterator;
    using index_set_t = typename multivector_t::index_set_t;

    auto result = Scalar_T(0);
    const bool small_star_large = lhs.size() < rhs.size();
    const multivector_t* smallp =
      small_star_large
      ? &lhs
      : &rhs;
    const multivector_t* largep =
      small_star_large
      ? &rhs
      : &lhs;

    for (auto
         small_it = smallp->begin();
         small_it != smallp->end();
         ++small_it)
    {
      const index_set_t small_ist = small_it->first;
      const Scalar_T    large_crd = (*largep)[small_ist];
      if (large_crd != Scalar_T(0))
        result += small_ist.sign_of_square() * small_it->second * large_crd;
    }
    return result;
  }

  /// Quotient of multivector and scalar
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator/= (const Scalar_T& scr) -> multivector_t&
  { // Divide coordinates of all terms by scr
    using traits_t = numeric_traits<Scalar_T>;

    if (traits_t::isNaN(scr))
      return *this = traits_t::NaN();
    if (traits_t::isInf(scr))
      if (this->isnan())
        *this = traits_t::NaN();
      else
        this->clear();
    else
      for (auto
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
        this_it->second /= scr;
    return *this;
  }

  /// Geometric quotient
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  operator/ (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> const framed_multi<Scalar_T,LO,HI>
  {
    using multivector_t = framed_multi<Scalar_T, LO, HI>;
    using traits_t = numeric_traits<Scalar_T>;
    using index_set_t = typename multivector_t::index_set_t;
    using matrix_multi_t = typename multivector_t::matrix_multi_t;

    if (rhs == Scalar_T(0))
      return traits_t::NaN();

    const index_set_t our_frame = lhs.frame() | rhs.frame();
    return matrix_multi_t(lhs, our_frame, true) / matrix_multi_t(rhs, our_frame, true);
  }

  /// Geometric quotient
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator/= (const multivector_t& rhs) -> multivector_t&
  { return *this = *this / rhs; }

  /// Transformation via twisted adjoint action
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  operator| (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> const framed_multi<Scalar_T,LO,HI>
  {
    using multivector_t = framed_multi<Scalar_T, LO, HI>;
    using matrix_multi_t = typename multivector_t::matrix_multi_t;

    return matrix_multi_t(rhs) * matrix_multi_t(lhs) / matrix_multi_t(rhs.involute());
  }

  /// Transformation via twisted adjoint action
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator|= (const multivector_t& rhs) -> multivector_t&
  { return *this = *this | rhs; }

  /// Clifford multiplicative inverse
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  inv() const -> const multivector_t
  {
    matrix_multi_t result = matrix_multi_t(Scalar_T(1), this->frame());
    return result /= matrix_multi_t(*this);
  }

  /// Integer power of multivector: *this to the m
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  pow(int m) const -> const multivector_t
  { return glucat::pow(*this, m); }

  /// Outer product power of multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  outer_pow(int m) const -> const multivector_t
  {
    if (m < 0)
      throw error_t("outer_pow(int): negative exponent");
    multivector_t result = Scalar_T(1);
    multivector_t a = *this;
    for (;
        m != 0;
        m >>= 1, a = a ^ a)
      if (m & 1)
        result ^= a;
    return result;
  }

  /// Frame of multivector: union of index sets of terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  frame() const -> const index_set_t
  {
    index_set_t result;
    for (auto
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
      result |= this_it->first;
    return result;
  }

  /// Grade of multivector: maximum of the grades of each term
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  grade() const -> index_t
  {
    index_t result = 0;
    for (auto
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
      result = std::max( result, this_it->first.count() );
    return result;
  }

  /// Subscripting: map from index set to scalar coordinate
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator[] (const index_set_t ist) const -> Scalar_T
  {
    const const_iterator& this_it = this->find(ist);
    if (this_it == this->end())
      return Scalar_T(0);
    else
      return this_it->second;
  }

  /// Grading: part where each term is a grade-vector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator() (index_t grade) const -> const multivector_t
  {
    if ((grade < 0) || (grade > HI-LO))
      return Scalar_T(0);
    else
    {
      multivector_t result;
      for (auto
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
        if (this_it->first.count() == grade)
          result += *this_it;
      return result;
    }
  }

  /// Scalar part
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  scalar() const -> Scalar_T
  { return (*this)[index_set_t()]; }

  /// Pure part
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  pure() const -> const multivector_t
  { return *this - this->scalar(); }

  /// Even part, sum of the even grade terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  even() const -> const multivector_t
  { // even part of x, sum of the pure(count) with even count
    multivector_t result;
    for (auto
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
      if ((this_it->first.count() % 2) == 0)
        result.insert(*this_it);
    return result;
  }

  /// Odd part, sum of the odd grade terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  odd() const -> const multivector_t
  { // even part of x, sum of the pure(count) with even count
    multivector_t result;
    for (auto
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
      if ((this_it->first.count() % 2) == 1)
        result.insert(*this_it);
    return result;
  }

  /// Vector part of multivector, as a vector_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  vector_part() const -> const vector_t
  { return this->vector_part(this->frame(), true); }

  /// Vector part of multivector, as a vector_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  vector_part(const index_set_t frm, const bool prechecked) const -> const vector_t
  {
    if (!prechecked && (this->frame() | frm) != frm)
      throw error_t("vector_part(frm): value is outside of requested frame");
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

  /// Main involution, each {i} is replaced by -{i} in each term
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  involute() const -> const multivector_t
  {
    multivector_t result = *this;
    for (auto
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
  auto
  framed_multi<Scalar_T,LO,HI>::
  reverse() const -> const multivector_t
  {
    multivector_t result = *this;
    for (auto
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
  auto
  framed_multi<Scalar_T,LO,HI>::
  conj() const -> const multivector_t
  {
    multivector_t result = *this;
    for (auto
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
  auto
  framed_multi<Scalar_T,LO,HI>::
  quad() const -> Scalar_T
  {
    // scalar(conj(x)*x) = 2*quad(even(x)) - quad(x)
    // ref: old clical: quadfunction(p:pter):pterm in file compmod.pas
    auto result = Scalar_T(0);
    for (auto
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
    {
      const Scalar_T sign =
        (this_it->first.count_neg() % 2)
        ? -Scalar_T(1)
        :  Scalar_T(1);
      result += sign * (this_it->second) * (this_it->second);
    }
    return result;
  }

  /// Norm squared := sum of norm squared of coordinates
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  norm() const -> Scalar_T
  {
    using traits_t = numeric_traits<Scalar_T>;

    auto result = Scalar_T(0);
    for (auto
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
    {
      const Scalar_T abs_crd = traits_t::abs(this_it->second);
      result +=  abs_crd * abs_crd;
    }
    return result;
  }

  /// Maximum of absolute values of components of multivector: multivector infinity norm
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  max_abs() const -> Scalar_T
  {
    using traits_t = numeric_traits<Scalar_T>;

    auto result = Scalar_T(0);
    for (auto
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
    {
      const Scalar_T abs_crd = traits_t::abs(this_it->second);
      if (abs_crd > result)
        result = abs_crd;
    }
    return result;
  }

  /// Random multivector within a frame
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  random(const index_set_t frm, Scalar_T fill) -> const multivector_t
  {
    using multivector_t = framed_multi<Scalar_T, LO, HI>;
    using index_set_t = typename multivector_t::index_set_t;
    using term_t = typename multivector_t::term_t;

    using random_generator_t = random_generator<Scalar_T>;
    random_generator_t& generator = random_generator_t::generator();

    fill =
      (fill < Scalar_T(0))
      ? Scalar_T(0)
      : (fill > Scalar_T(1))
        ? Scalar_T(1)
        : fill;
    const set_value_t algebra_dim = 1 << frm.count();
    using traits_t = numeric_traits<Scalar_T>;
    const Scalar_T mean_abs = traits_t::sqrt(Scalar_T(double(algebra_dim)));
    multivector_t result;
    for (set_value_t
        stv = 0;
        stv != algebra_dim;
        ++stv)
      if (generator.uniform() < fill)
      {
        const Scalar_T& result_crd = generator.normal() / mean_abs;
        result.insert(term_t(index_set_t(stv, frm, true), result_crd));
      }
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
    using map_t = Map_T;
    using sorted_map_t = Sorted_Map_T;
    using sorted_iterator = typename Sorted_Map_T::const_iterator;

    sorted_range (Sorted_Map_T &sorted_val, const Map_T& val)
    {
      for (auto
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
    using map_t = Sorted_Map_T;
    using sorted_map_t = Sorted_Map_T;
    using sorted_iterator = typename Sorted_Map_T::const_iterator;

    sorted_range (Sorted_Map_T &sorted_val, const Sorted_Map_T& val)
    : sorted_begin( val.begin() ),
      sorted_end( val.end() )
    { }
    sorted_iterator sorted_begin;
    sorted_iterator sorted_end;
  };

  /// Write multivector to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator<< (std::ostream& os, const framed_multi<Scalar_T,LO,HI>& val) -> std::ostream&
  {
    if (val.empty())
      os << 0;
    else if (val.isnan())
      os << std::numeric_limits<Scalar_T>::quiet_NaN();
    else if (val.isinf())
    {
      const Scalar_T& inf = std::numeric_limits<Scalar_T>::infinity();
      os << (scalar(val) < 0.0 ? -inf : inf);
    }
    else
    {
      using multivector_t = framed_multi<Scalar_T, LO, HI>;
      using map_t = typename multivector_t::map_t;
      using sorted_map_t = typename multivector_t::sorted_map_t;
      using sorted_iterator = typename sorted_map_t::const_iterator;
      sorted_map_t sorted_val;
      sorted_range< map_t, sorted_map_t > sorted_val_range(sorted_val, val);
      auto sorted_it = sorted_val_range.sorted_begin;
      os << *sorted_it;
      for (++sorted_it;
          sorted_it != sorted_val_range.sorted_end;
          ++sorted_it)
      {
        const Scalar_T& scr = sorted_it->second;
        if (scr >= 0.0)
          os << '+';
        os << *sorted_it;
      }
    }
    return os;
  }

  /// Write term to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator<< (std::ostream& os, const std::pair< const index_set<LO,HI>, Scalar_T >& term) -> std::ostream&
  {
    const double second_as_double = numeric_traits<Scalar_T>::to_double(term.second);
    const bool use_double =
      (os.precision() <= std::numeric_limits<double>::digits10) ||
      (term.second == Scalar_T(second_as_double));
    if (term.first.count() == 0)
      if (use_double)
        os << second_as_double;
      else
        os << term.second;
    else if (term.second == Scalar_T(-1))
    {
      os << '-';
      os << term.first;
    }
    else if (term.second != Scalar_T(1))
    {
      if (use_double)
      {
        double tol = std::pow(10.0,-os.precision());
        if ( std::fabs(second_as_double + 1.0) < tol )
          os << '-';
        else if ( std::fabs(second_as_double - 1.0) >= tol )
          os << second_as_double;
      }
      else
        os << term.second;
      os << term.first;
    }
    else
      os << term.first;
    return os;
  }

  /// Read multivector from input
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator>> (std::istream& s, framed_multi<Scalar_T,LO,HI> & val) -> std::istream&
  { // Input looks like 1.0-2.0{1,2}+3.2{3,4}.
    using multivector_t = framed_multi<Scalar_T, LO, HI>;
    // Parsing variables.
    multivector_t local_val;
    int c = 0;
    // Parsing control variables.
    bool negative = false;
    bool expect_term = true;
    // The multivector may begin with '+' or '-'. Check for this.
    c = s.peek();
    if (s.good() && (c == int('+') || c == int('-')))
    { // A '-' here negates the following term.
      negative = (c == int('-'));
      // Consume the '+' or '-'.
      s.get();
    }
    while (s.good())
    { // Parse a term.
      // A term consists of an optional scalar, followed by an optional index set.
      // At least one of the two must be present.
      // Default coordinate is Scalar_T(1).
      auto  coordinate  = Scalar_T(1);
      // Default index set is empty.
      index_set<LO,HI> ist;
      // First, check for an opening brace.
      c = s.peek();
      if (s.good())
      { // If the character is not an opening brace,
        // a coordinate value is expected here.
        if (c != int('{'))
        { // Try to read a coordinate value.
          double coordinate_as_double;
          s >> coordinate_as_double;
          // Reading the coordinate may have resulted in an end of file condition.
          // This is not a failure.
          if (s)
            coordinate = Scalar_T(coordinate_as_double);
        }
      }
      else
      { // End of file here ends parsing while a term may still be expected.
        break;
      }
      // Coordinate is now Scalar_T(1) or a Scalar_T value.
      // Parse an optional index set.
      if (s.good())
      {
        c = s.peek();
        if (s.good() && c == int('{'))
        { // Try to read index set.
          s >> ist;
        }
      }
      // Reading the term may have resulted in an end of file condition.
      // This is not a failure.
      if (s)
      {
        // Immediately after parsing a term, another term is not expected.
        expect_term = false;
        if (coordinate != Scalar_T(0))
        {
          // Add the term to the local multivector.
          coordinate =
            negative
            ? -coordinate
            :  coordinate;
          using term_t = typename multivector_t::term_t;
          local_val += term_t(ist, coordinate);
        }
      }
      // Check if anything follows the current term.
      if (s.good())
      {
        c = s.peek();
        if (s.good())
        { // Only '+' and '-' are valid here.
          if (c == int('+') || c == int('-'))
          { // A '-' here negates the following term.
            negative = (c == int('-'));
            // Consume the '+' or '-'.
            s.get();
            // Immediately after '+' or '-',
            // expect another term.
            expect_term = true;
          }
          else
          { // Any other character here is a not failure,
            // but still ends the parsing of the multivector.
            break;
          }
        }
      }
    }
    // If a term is still expected, this is a failure.
    if (expect_term)
      s.clear(std::istream::failbit);
    // End of file is not a failure.
    if (s)
    { // The multivector has been successfully parsed.
      val = local_val;
    }
    return s;
  }

  /// Number of terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  nbr_terms () const -> unsigned long
  { return this->size(); }

  /// Insert a term into a multivector, add terms with same index set.
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  operator+= (const term_t& term) -> multivector_t&
  { // Do not insert terms with 0 coordinate
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

  /// Check if a multivector contains any infinite values
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  isinf() const -> bool
  {
    using traits_t = numeric_traits<Scalar_T>;

    if (std::numeric_limits<Scalar_T>::has_infinity)
      for (auto
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
          if (traits_t::isInf(this_it->second))
            return true;
    return false;
  }

  /// Check if a multivector contains any IEEE NaN values
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  isnan() const -> bool
  {
    using traits_t = numeric_traits<Scalar_T>;

    if (std::numeric_limits<Scalar_T>::has_quiet_NaN)
      for (auto
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
          if (traits_t::isNaN(this_it->second))
            return true;
    return false;
  }

  /// Remove all terms with relative size smaller than limit
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  truncated(const Scalar_T& limit) const -> const multivector_t
  {
    using traits_t = numeric_traits<Scalar_T>;

    const Scalar_T abs_limit = traits_t::abs(limit);
    if (this->isnan())
      return *this;
    Scalar_T top = max_abs();
    multivector_t result;
    if (top != Scalar_T(0))
      for (auto
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
        if (traits_t::abs(this_it->second / top) > abs_limit)
          result.insert(term_t(this_it->first, this_it->second));
    return result;
  }

  /// Subalgebra isomorphism: fold each term within the given frame
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  fold(const index_set_t frm) const -> multivector_t
  {
    if (frm.is_contiguous())
      return *this;
    else
    {
      multivector_t result;
      for (auto
          this_it = this->begin();
          this_it != this->end();
          ++this_it)
        result.insert(term_t(this_it->first.fold(frm), this_it->second));
      return result;
    }
  }

  /// Subalgebra isomorphism: unfold each term within the given frame
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  unfold(const index_set_t frm) const -> multivector_t
  {
    if (frm.is_contiguous())
      return *this;
    else
    {
      multivector_t result;
      for (auto
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
  auto
  framed_multi<Scalar_T,LO,HI>::
  centre_pm4_qp4(index_t& p, index_t& q) -> multivector_t&
  {
    // We add 4 to q by subtracting 4 from p
    if (q+4 > -LO)
      throw error_t("centre_pm4_qp4(p,q): LO is too high to represent this value");
    if (this->frame().max() > p-4)
    {
      using index_pair_t = typename index_set_t::index_pair_t;
      const index_set_t pm3210(index_pair_t(p-3,p), true);
      const index_set_t qm4321(index_pair_t(-q-4,-q-1), true);
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
  auto
  framed_multi<Scalar_T,LO,HI>::
  centre_pp4_qm4(index_t& p, index_t& q) -> multivector_t&
  {
    // We add 4 to p by subtracting 4 from q
    if (p+4 > HI)
      throw error_t("centre_pp4_qm4(p,q): HI is too low to represent this value");
    if (this->frame().min() < -q+4)
    {
      using index_pair_t = typename index_set_t::index_pair_t;
      const index_set_t qp0123(index_pair_t(-q,-q+3), true);
      const index_set_t pp1234(index_pair_t(p+1,p+4), true);
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
  auto
  framed_multi<Scalar_T,LO,HI>::
  centre_qp1_pm1(index_t& p, index_t& q) -> multivector_t&
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
  auto
  framed_multi<Scalar_T,LO,HI>::
  divide(const index_set_t ist) const -> const framed_pair_t
  {
    multivector_t quo;
    multivector_t rem;
    for (auto
        this_it = this->begin();
        this_it != this->end();
        ++this_it)
      if ((this_it->first | ist) == this_it->first)
        quo.insert(term_t(this_it->first ^ ist, this_it->second));
      else
        rem.insert(*this_it);
    return framed_pair_t(quo, rem);
  }

  /// Generalized FFT from multivector_t to matrix_t
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  framed_multi<Scalar_T,LO,HI>::
  fast(const index_t level, const bool odd) const -> const matrix_t
  {
    // Assume val is already folded and centred
    if (this->empty())
    {
      using matrix_index_t = typename matrix_multi_t::matrix_index_t;
      const matrix_index_t dim = 1 << level;
      matrix_t result(dim, dim);
      result.clear();
      return result;
    }
    if (level == 0)
      return matrix::unit<matrix_t>(1) * this->scalar();

    using basis_matrix_t = typename matrix_multi_t::basis_matrix_t;
    using basis_scalar_t = typename basis_matrix_t::value_type;

    const basis_matrix_t&  I = matrix::unit<basis_matrix_t>(2);
    basis_matrix_t J(2,2,2);
    J.clear();
    J(0,1)  = basis_scalar_t(-1);
    J(1,0)  = basis_scalar_t( 1);
    basis_matrix_t K = J;
    K(0,1)  = basis_scalar_t( 1);
    basis_matrix_t JK = I;
    JK(0,0) = basis_scalar_t(-1);

    const index_set_t ist_mn = index_set_t(-level);
    const index_set_t ist_pn = index_set_t(level);
    if (level == 1)
    {
      if (odd)
        return matrix_t(J) * (*this)[ist_mn] + matrix_t(K)  * (*this)[ist_pn];
      else
        return matrix_t(I) * this->scalar()  + matrix_t(JK) * (*this)[ist_mn ^ ist_pn];
    }
    else
    {
      const framed_pair_t& pair_mn = this->divide(ist_mn);
      const multivector_t& quo_mn = pair_mn.first;
      const multivector_t& rem_mn = pair_mn.second;
      const framed_pair_t& pair_quo_mnpn = quo_mn.divide(ist_pn);
      const multivector_t& val_mnpn = pair_quo_mnpn.first;
      const multivector_t& val_mn   = pair_quo_mnpn.second;
      const framed_pair_t& pair_rem_mnpn = rem_mn.divide(ist_pn);
      const multivector_t& val_pn   = pair_rem_mnpn.first;
      const multivector_t& val_1    = pair_rem_mnpn.second;
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
  template< typename Other_Scalar_T >
  auto
  framed_multi<Scalar_T,LO,HI>::
  fast_matrix_multi(const index_set_t frm) const -> const matrix_multi<Other_Scalar_T,LO,HI>
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
      return matrix_multi<Other_Scalar_T,LO,HI>(ev_val.fast(level, 0) + od_val.fast(level, 1), frm);
    }

  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  framed_multi<Scalar_T,LO,HI>::
  fast_framed_multi() const -> const multivector_t
  { return *this; }

  /// Coordinate of product of terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  static
  auto
  crd_of_mult(const std::pair<const index_set<LO,HI>, Scalar_T>& lhs,
              const std::pair<const index_set<LO,HI>, Scalar_T>& rhs) -> Scalar_T
  { return lhs.first.sign_of_mult(rhs.first) * lhs.second * rhs.second; }

  /// Product of terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  operator* (const std::pair<const index_set<LO,HI>, Scalar_T>& lhs,
             const std::pair<const index_set<LO,HI>, Scalar_T>& rhs) -> const std::pair<const index_set<LO,HI>, Scalar_T>
  {
    using term_t = std::pair<const index_set<LO, HI>, Scalar_T>;
    return term_t(lhs.first ^ rhs.first, crd_of_mult(lhs, rhs));
  }

  /// Square root of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  sqrt(const framed_multi<Scalar_T,LO,HI>& val, const framed_multi<Scalar_T,LO,HI>& i, bool prechecked) -> const framed_multi<Scalar_T,LO,HI>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    check_complex(val, i, prechecked);

    const Scalar_T realval = val.scalar();
    if (val == realval)
    {
      if (realval < Scalar_T(0))
        return i * traits_t::sqrt(-realval);
      else
        return traits_t::sqrt(realval);
    }
    using matrix_multi_t = typename framed_multi<Scalar_T, LO, HI>::matrix_multi_t;
    return sqrt(matrix_multi_t(val), matrix_multi_t(i), prechecked);
  }

  /// Exponential of multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  exp(const framed_multi<Scalar_T,LO,HI>& val) -> const framed_multi<Scalar_T,LO,HI>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    const Scalar_T s = scalar(val);
    if (val == s)
      return traits_t::exp(s);

    const double size = val.size();
    const index_t frm_count = val.frame().count();
    const set_value_t algebra_dim = 1 << frm_count;

    if( (size * size <= double(algebra_dim)) || (frm_count < Tune_P::mult_matrix_threshold))
    {
      switch (Tune_P::function_precision)
      {
      case precision_demoted:
        {
          using demoted_scalar_t = typename traits_t::demoted::type;
          using demoted_multivector_t = framed_multi<demoted_scalar_t, LO, HI>;

          const demoted_multivector_t& demoted_val = demoted_multivector_t(val);
          return clifford_exp(demoted_val);
        }
        break;
      case precision_promoted:
        {
          using promoted_scalar_t = typename traits_t::promoted::type;
          using promoted_multivector_t = framed_multi<promoted_scalar_t, LO, HI>;

          const promoted_multivector_t& promoted_val = promoted_multivector_t(val);
          return clifford_exp(promoted_val);
        }
        break;
      default:
        return clifford_exp(val);
      }
    }
    else
    {
      using matrix_multi_t = matrix_multi<Scalar_T, LO, HI>;
      return exp(matrix_multi_t(val));
    }
  }

  /// Natural logarithm of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  log(const framed_multi<Scalar_T,LO,HI>& val, const framed_multi<Scalar_T,LO,HI>& i, bool prechecked) -> const framed_multi<Scalar_T,LO,HI>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();

    check_complex(val, i, prechecked);

    const Scalar_T realval = val.scalar();
    if (val == realval)
    {
      if (realval < Scalar_T(0))
        return i * traits_t::pi() + traits_t::log(-realval);
      else
        return traits_t::log(realval);
    }
    using matrix_multi_t = typename framed_multi<Scalar_T, LO, HI>::matrix_multi_t;
    return log(matrix_multi_t(val), matrix_multi_t(i), prechecked);
  }
}
#endif  // _GLUCAT_FRAMED_MULTI_IMP_H
