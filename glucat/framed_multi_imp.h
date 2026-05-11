#ifndef _GLUCAT_FRAMED_MULTI_IMP_H
#define _GLUCAT_FRAMED_MULTI_IMP_H
/**************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    framed_multi_imp.h : Implement the coordinate map representation of a
    Clifford algebra element
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2021 by Paul C. Leopardi
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
 ***************************************************************************
 ***************************************************************************/

#include "glucat/framed_multi.h"

#include "glucat/scalar.h"
#include "glucat/random.h"
#include "glucat/generation.h"
#include "glucat/matrix.h"

#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>

namespace glucat
{
  /*
   * @brief Class name used in messages
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  classname() -> std::string_view
  { return "framed_multi"; }

#define _GLUCAT_HASH_N(x) (x)
#define _GLUCAT_HASH_SIZE_T(x) (typename multivector_t::hash_size_t)(x)

  /*
   * @brief Default constructor
   * @details
   *
   * Usage example:
   * Location: glucat/framed_multi_imp.h:177
   *
   * @code
   *
   * framed_multi()
   * @endcode
   *
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi()
  : map_t(_GLUCAT_HASH_N(0))
  { }

  /*
   * @brief Move constructor
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param other Other matrix
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi(framed_multi&& other) noexcept(std::is_nothrow_move_constructible_v<Scalar_T>)
  : map_t(std::move(other))
  { }

  /*
   * @brief Private constructor using hash_size
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param hash_size Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi(const hash_size_t& hash_size)
  : map_t(_GLUCAT_HASH_N(hash_size()))
  { }

  /*
   * @brief Construct a multivector from a multivector with a different scalar type
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param val Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  template< typename Other_Scalar_T, typename Other_Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi(const framed_multi<Other_Scalar_T,LO,HI,Other_Tune_P>& val)
  : map_t(_GLUCAT_HASH_N(val.size()))
  {
    for (auto& val_term : val)
      this->insert(term_t(val_term.first, numeric_traits<Scalar_T>::to_scalar_t(val_term.second)));
  }

  /*
   * @brief Construct a multivector, within a given frame, from a given multivector
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param val Value
   * @param frm Frame
   * @param prechecked Bool: true if i has already been checked?
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  template< typename Other_Scalar_T, typename Other_Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi(const framed_multi<Other_Scalar_T,LO,HI,Other_Tune_P>& val,
               const index_set_t frm, const bool prechecked)
  : map_t(_GLUCAT_HASH_N(val.size()))
  {
    if (!prechecked && (val.frame() | frm) != frm)
      throw error_t("multivector_t(val,frm): cannot initialize with value outside of frame");
    for (auto& val_term : val)
      this->insert(term_t(val_term.first, numeric_traits<Scalar_T>::to_scalar_t(val_term.second)));
  }

  /*
   * @brief Construct a multivector, within a given frame, from a given multivector
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param val Value
   * @param frm Frame
   * @param prechecked Bool: true if i has already been checked?
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi(const multivector_t& val,
               const index_set_t frm, const bool prechecked)
  : map_t(_GLUCAT_HASH_N(val.size()))
  {
    if (!prechecked && (val.frame() | frm) != frm)
      throw error_t("multivector_t(val,frm): cannot initialize with value outside of frame");
    for (auto& val_term : val)
      this->insert(val_term);
  }

  /*
   * @brief Construct a multivector from an index set and a scalar coordinate
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param ist Value
   * @param crd Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi(const index_set_t ist, const Scalar_T& crd)
  : map_t(_GLUCAT_HASH_N(1))
  {
    if (crd != Scalar_T(0))
      this->insert(term_t(ist, crd));
  }

  /*
   * @brief Construct a multivector, within a given frame, from an index set and a scalar coordinate
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param ist Value
   * @param crd Value
   * @param frm Frame
   * @param prechecked Bool: true if i has already been checked?
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi(const index_set_t ist, const Scalar_T& crd,
               const index_set_t frm, const bool prechecked)
  : map_t(_GLUCAT_HASH_N(1))
  {
    if (!prechecked && (ist | frm) != frm)
      throw error_t("multivector_t(ist,crd,frm): cannot initialize with value outside of frame");
    if (crd != Scalar_T(0))
      this->insert(term_t(ist, crd));
  }

  /*
   * @brief Construct a multivector from a scalar (within a frame, if given)
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param scr Value
   * @param frm Frame
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi(const Scalar_T& scr, const index_set_t frm)
  : map_t(_GLUCAT_HASH_N(1))
  {
    if (scr != Scalar_T(0))
      this->insert(term_t(index_set_t(), scr));
  }

  /*
   * @brief Construct a multivector from an int (within a frame, if given)
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param scr Value
   * @param frm Frame
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi(const int scr, const index_set_t frm)
  : map_t(_GLUCAT_HASH_N(1))
  {
    if (scr != Scalar_T(0))
      this->insert(term_t(index_set_t(), Scalar_T(scr)));
  }

  /*
   * @brief Construct a multivector, within a given frame, from a given vector
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param vec Value
   * @param frm Frame
   * @param prechecked Bool: true if i has already been checked?
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi(const vector_t& vec,
               const index_set_t frm, const bool prechecked)
  : map_t(_GLUCAT_HASH_N(vec.size()))
  {
    if (!prechecked && index_t(vec.size()) != frm.count())
      throw error_t("multivector_t(vec,frm): cannot initialize with vector not matching frame");
    auto idx = frm.min();
    const auto frm_end = frm.max()+1;
    for (auto& crd : vec)
    {
      *this += term_t(index_set_t(idx), crd);
      for (
        ++idx;
        idx != frm_end && !frm[idx];
        ++idx)
        ;
    }
  }

  /*
   * @brief Construct a multivector from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param str Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
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

  /*
   * @brief Construct a multivector, within a given frame, from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param str Value
   * @param frm Frame
   * @param prechecked Bool: true if i has already been checked?
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi(const std::string& str, const index_set_t frm, const bool prechecked)
  : map_t(_GLUCAT_HASH_N(0))
  {
    if (prechecked)
      *this = multivector_t(str);
    else
      *this = multivector_t(multivector_t(str), frm, false);
  }

  /*
   * @brief Construct a multivector from a matrix_multi_t
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param val Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  template< typename Other_Scalar_T, typename Other_Tune_P >
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  framed_multi(const matrix_multi<Other_Scalar_T,LO,HI,Other_Tune_P>& val)
  : map_t(_GLUCAT_HASH_N(1))
  {
    if (val == Other_Scalar_T(0))
      return;

    const auto dim = val.m_matrix.nbr_rows();
    using traits_t = numeric_traits<Scalar_T>;
    if (dim == 1)
    {
      this->insert(term_t(index_set_t(), traits_t::to_scalar_t(val.m_matrix(0, 0))));
      return;
    }
    using Tuning_Values_P = typename Tune_P::tuning_values_p;
    if (dim >= Tuning_Values_P::inv_fast_dim_threshold)
      try
      {
        *this = (val.template fast_framed_multi<Scalar_T,Tune_P>()).truncated();
        return;
      }
      catch (const std::exception& e)
      { } // Fall back to the slow algorithm


    const auto val_norm = traits_t::to_scalar_t(val.norm());
    if (traits_t::isNaN_or_isInf(val_norm))
    {
      *this = multivector_t(traits_t::NaN());
      return;
    }
    const auto frm = val.frame();
    const auto algebra_dim = set_value_t(1) << frm.count();
    auto result = multivector_t(
      _GLUCAT_HASH_SIZE_T(std::min<size_t>(algebra_dim, val.m_matrix.nnz())));
    for (auto
        stv = set_value_t(0);
        stv != algebra_dim;
        stv++)
    {
      const auto ist = index_set_t(stv, frm, true);
      const auto crd =
        traits_t::to_scalar_t(val.basis_element(ist).template inner<Other_Scalar_T>(val.m_matrix));
      if (crd != Scalar_T(0))
        result.insert(term_t(ist, crd));
    }
    *this = result.truncated();
  }

  /*
   * @brief Test for equality of multivectors
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param rhs Right hand side
   * @return True if equal
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator==  (const multivector_t& rhs) const -> bool
  {
    if (this->size() != rhs.size())
      return false;
    const auto rhs_end = rhs.end();
    for (auto& this_term : *this)
    {
      const const_iterator& rhs_it = rhs.find(this_term.first);
      if (rhs_it == rhs_end || rhs_it->second != this_term.second)
        return false;
    }
    return true;
  }

  /*
   * @brief Test for equality of multivector and scalar
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param scr Value
   * @return True if equal
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator==  (const Scalar_T& scr) const -> bool
  {
    switch (this->size())
    {
    case 0:
      return scr == Scalar_T(0);
    case 1:
      {
        const auto& this_it = this->begin();
        return this_it->first == index_set_t() && this_it->second == scr;
      }
    default:
      return false;
    }
  }

  /*
   * @brief Add scalar
   * @details
   *
   * Usage example:
   * Location: glucat/clifford_algebra_imp.h:237
   *
   * @code
   *
   * return result += scr;
   * @endcode
   *
   * @param scr Scalar
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator+= (const Scalar_T& scr) -> multivector_t&
  {
    *this += term_t(index_set_t(), scr);
    return *this;
  }

  /*
   * @brief Geometric sum
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param rhs Right hand side
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator+= (const multivector_t& rhs) -> multivector_t&
  { // simply add terms
    for (auto& rhs_term : rhs)
      *this += rhs_term;
    return *this;
  }

  _GLUCAT_CLIFFORD_ALGEBRA_ASSIGNMENT_OPERATIONS_IMP(framed_multi)

  /*
   * @brief Subtract scalar
   * @details
   *
   * Usage example:
   * Location: glucat/clifford_algebra_imp.h:291
   *
   * @code
   *
   * return result -= scr;
   * @endcode
   *
   * @param scr Scalar
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator-= (const Scalar_T& scr) -> multivector_t&
  {
    *this += term_t(index_set_t(), -scr);
    return *this;
  }

  /*
   * @brief Geometric difference
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param rhs Right hand side
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator-= (const multivector_t& rhs) -> multivector_t&
  {
    for (auto& rhs_term : rhs)
      *this += term_t(rhs_term.first, -(rhs_term.second));
    return *this;
  }

  /*
   * @brief Unary minus
   * @details
   *
   * Usage example:
   * Location: test11/peg11.h:277
   *
   * @code
   *
   * transcendtest(-m_(1));
   * @endcode
   *
   * @return Unary minus
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator- () const -> multivector_t
  { // multiply coordinates of all terms by -1
    auto result = *this;
    for (auto& result_term : result)
      const_cast<Scalar_T&>(result_term.second) *= Scalar_T(-1);
    return result;
  }

  /*
   * @brief Product of multivector and scalar
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param scr Value
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
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
      for (auto& this_term : *this)
        const_cast<Scalar_T&>(this_term.second) *= scr;
    return *this;
  }

  /*
   * @brief Geometric product
   * @details
   *
   * Usage example:
   * Location: test00/peg00.h:114
   *
   * @code
   *
   * rhs = (a_r * b_s)(index_t(std::abs(r-s)));
   * @endcode
   *
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Product
   *
   * @par Example:
   * @code
   * clifford<>("{1}") * clifford<>("{2}"); // Returns {1,2}
   * @endcode
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  operator* (const framed_multi<Scalar_T,LO,HI,Tune_P>& lhs, const framed_multi<Scalar_T,LO,HI,Tune_P>& rhs) -> framed_multi<Scalar_T,LO,HI,Tune_P>
  {
    using multivector_t = framed_multi<Scalar_T,LO,HI,Tune_P>;
    using traits_t = numeric_traits<Scalar_T>;
    using term_t = typename multivector_t::term_t;
    using index_set_t = typename multivector_t::index_set_t;

    if (lhs.isnan() || rhs.isnan())
      return traits_t::NaN();

    const double lhs_size = lhs.size();
    const double rhs_size = rhs.size();
    const auto our_frame = lhs.frame() | rhs.frame();
    const auto frm_count = our_frame.count();
    const auto algebra_dim = set_value_t(1) << frm_count;
    using Tuning_Values_P = typename Tune_P::tuning_values_p;
    const auto direct_mult = lhs_size * rhs_size <= double(algebra_dim);
    if (direct_mult)
    { // If we have a sparse multiply, store the result directly
      auto result = multivector_t(
        _GLUCAT_HASH_SIZE_T(size_t(std::min(lhs_size * rhs_size, double(algebra_dim)))));
      for (auto& lhs_term : lhs)
        for (auto& rhs_term : rhs)
          result += term_t(lhs_term) * term_t(rhs_term);
      return result;
    }
    else if (frm_count < index_t(Tuning_Values_P::mult_matrix_threshold))
    {
      const set_value_t dim = algebra_dim;
      std::vector<Scalar_T> l_vec(dim, Scalar_T(0));
      for (auto& term : lhs) l_vec[term.first.value_of_fold(our_frame)] = term.second;
      std::vector<Scalar_T> r_vec(dim, Scalar_T(0));
      for (auto& term : rhs) r_vec[term.first.value_of_fold(our_frame)] = term.second;

      std::vector<index_set_t> ists(dim);
      std::vector<set_value_t> l_indices, r_indices;
      for (set_value_t k = 0; k < dim; ++k)
      {
        ists[k] = index_set_t(k, our_frame, true);
        if (l_vec[k] != Scalar_T(0)) l_indices.push_back(k);
        if (r_vec[k] != Scalar_T(0)) r_indices.push_back(k);
      }

      std::vector<Scalar_T> res_vec(dim, Scalar_T(0));
      for (auto i : l_indices)
      {
        const auto& ist_i = ists[i];
        for (auto j : r_indices)
        {
          const auto& ist_j = ists[j];
          Scalar_T sign = traits_t::to_scalar_t(ist_i.sign_of_mult(ist_j));
          res_vec[i ^ j] += sign * l_vec[i] * r_vec[j];
        }
      }
      auto result = multivector_t(_GLUCAT_HASH_SIZE_T(dim));
      for (set_value_t k = 0; k < dim; ++k)
        if (res_vec[k] != Scalar_T(0))
          result.insert(term_t(index_set_t(k, our_frame, true), res_vec[k]));
      return result.truncated();
    }
    else
    { // Past a certain threshold, the matrix algorithm is fastest
      using matrix_multi_t = typename multivector_t::matrix_multi_t;
      return multivector_t(matrix_multi_t(lhs, our_frame, true) *
                           matrix_multi_t(rhs, our_frame, true));
    }
  }

  /*
   * @brief Geometric product
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param rhs Right hand side
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator*= (const multivector_t& rhs) -> multivector_t&
  { return *this = *this * rhs; }

  /*
   * @brief Outer product
   * @details
   *
   * Usage example:
   * Location: test00/peg00.h:176
   *
   * @code
   *
   * lhs = (a ^ b) ^ c;
   * @endcode
   *
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Outer product
   *
   * @par Example:
   * @code
   * clifford<>("{1}") ^ clifford<>("{2}"); // Returns {1,2}
   * @endcode
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  operator^ (const framed_multi<Scalar_T,LO,HI,Tune_P>& lhs, const framed_multi<Scalar_T,LO,HI,Tune_P>& rhs) -> framed_multi<Scalar_T,LO,HI,Tune_P>
  { // Arvind Raja's original reference:
    // "old clical, outerproduct(p,q:pterm):pterm in file compmod.pas"

    if (lhs.empty() || rhs.empty())
      return Scalar_T(0);

    using multivector_t = framed_multi<Scalar_T,LO,HI,Tune_P>;
    using index_set_t = typename multivector_t::index_set_t;
    using term_t = typename multivector_t::term_t;

    const auto our_frame = lhs.frame() | rhs.frame();
    const auto algebra_dim = set_value_t(1) << our_frame.count();
    const double lhs_size = lhs.size();
    const double rhs_size = rhs.size();
    using Tuning_Values_P = typename Tune_P::tuning_values_p;
    const auto frm_count = our_frame.count();
    if (lhs_size * rhs_size <= double(algebra_dim))
    {
      const auto empty_set = index_set_t();
      auto result = multivector_t(
        _GLUCAT_HASH_SIZE_T(size_t(std::min(lhs_size * rhs_size, double(algebra_dim)))));
      for (auto& lhs_term : lhs)
        for (auto& rhs_term : rhs)
          if ((lhs_term.first & rhs_term.first) == empty_set)
            result += term_t(lhs_term) * term_t(rhs_term);
      return result;
    }
    else if (frm_count < index_t(Tuning_Values_P::products_matrix_threshold))
    {
      const set_value_t dim = algebra_dim;
      std::vector<Scalar_T> l_vec(dim, Scalar_T(0));
      for (auto& term : lhs) l_vec[term.first.value_of_fold(our_frame)] = term.second;
      std::vector<Scalar_T> r_vec(dim, Scalar_T(0));
      for (auto& term : rhs) r_vec[term.first.value_of_fold(our_frame)] = term.second;

      std::vector<set_value_t> unfolded_bits(dim);
      std::vector<unsigned long> h_vals(dim);
      for (set_value_t k = 0; k < dim; ++k)
      {
        unfolded_bits[k] = index_set_t(k, our_frame, true).to_set_value();
        h_vals[k] = inverse_reversed_gray(unfolded_bits[k]);
      }

      std::vector<Scalar_T> res_vec(dim, Scalar_T(0));
      for (set_value_t k = 0; k < dim; ++k)
      {
        for (set_value_t i = k; ; i = (i - 1) & k)
        {
          set_value_t j = k ^ i;
          if (l_vec[i] != Scalar_T(0) && r_vec[j] != Scalar_T(0))
          {
            unsigned long uthis = unfolded_bits[i];
            unsigned long h = h_vals[j];
            unsigned long neg = inverse_gray(uthis & h);
            Scalar_T sign = (neg & 1) ? Scalar_T(-1) : Scalar_T(1);
            res_vec[k] += sign * l_vec[i] * r_vec[j];
          }
          if (i == 0) break;
        }
      }
      auto result = multivector_t(_GLUCAT_HASH_SIZE_T(dim));
      for (set_value_t k = 0; k < dim; ++k)
        if (res_vec[k] != Scalar_T(0))
          result.insert(term_t(index_set_t(k, our_frame, true), res_vec[k]));
      return result;
    }
    else
    {
      using matrix_multi_t = typename multivector_t::matrix_multi_t;
      return multivector_t(matrix_multi_t(lhs, our_frame, true) ^
                           matrix_multi_t(rhs, our_frame, true));
    }
  }

  /*
   * @brief Add multivector
   * @details
   *
   * Usage example:
   * Location: glucat/framed_multi_imp.h:635
   *
   * @code
   *
   * result += lhs_term * rhs_term;
   * @endcode
   *
   * @param rhs Right hand side
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator^= (const multivector_t& rhs) -> multivector_t&
  { return *this = *this ^ rhs; }

  /*
   * @brief Inner product
   * @details
   *
   * Usage example:
   * Location: test00/peg00.h:113
   *
   * @code
   *
   * lhs = a_r & b_s;
   * @endcode
   *
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Inner product
   *
   * @par Example:
   * @code
   * clifford<>("{1}") & clifford<>("{1}"); // Returns 1
   * @endcode
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  operator& (const framed_multi<Scalar_T,LO,HI,Tune_P>& lhs, const framed_multi<Scalar_T,LO,HI,Tune_P>& rhs) -> framed_multi<Scalar_T,LO,HI,Tune_P>
  { // Arvind Raja's original reference:
    // "old clical, innerproduct(p,q:pterm):pterm in file compmod.pas"

    if (lhs.empty() || rhs.empty())
      return Scalar_T(0);

    using multivector_t = framed_multi<Scalar_T,LO,HI,Tune_P>;
    using index_set_t = typename multivector_t::index_set_t;
    using term_t = typename multivector_t::term_t;

    const auto our_frame = lhs.frame() | rhs.frame();
    const auto algebra_dim = set_value_t(1) << our_frame.count();
    const double lhs_size = lhs.size();
    const double rhs_size = rhs.size();
    using Tuning_Values_P = typename Tune_P::tuning_values_p;
    const auto frm_count = our_frame.count();
    if (lhs_size * rhs_size <= double(algebra_dim))
    {
      auto result = multivector_t(
        _GLUCAT_HASH_SIZE_T(size_t(std::min(lhs_size * rhs_size, double(algebra_dim)))));
      const auto empty_set = index_set_t();
      for (auto& lhs_term : lhs)
      {
        const auto lhs_ist = lhs_term.first;
        if (lhs_ist != empty_set)
          for (auto& rhs_term : rhs)
          {
            const auto rhs_ist = rhs_term.first;
            if (rhs_ist != empty_set)
            {
              const auto our_ist = lhs_ist | rhs_ist;
              if ((lhs_ist == our_ist) || (rhs_ist == our_ist))
                result += term_t(lhs_term) * term_t(rhs_term);
            }
          }
      }
      return result;
    }
    else if (frm_count < index_t(Tuning_Values_P::products_matrix_threshold))
    {
      const set_value_t dim = algebra_dim;
      std::vector<Scalar_T> l_vec(dim, Scalar_T(0));
      for (auto& term : lhs) l_vec[term.first.value_of_fold(our_frame)] = term.second;
      std::vector<Scalar_T> r_vec(dim, Scalar_T(0));
      for (auto& term : rhs) r_vec[term.first.value_of_fold(our_frame)] = term.second;

      std::vector<set_value_t> unfolded_bits(dim);
      std::vector<unsigned long> h_vals(dim);
      for (set_value_t k = 0; k < dim; ++k)
      {
        unfolded_bits[k] = index_set_t(k, our_frame, true).to_set_value();
        h_vals[k] = inverse_reversed_gray(unfolded_bits[k]);
      }

      std::vector<Scalar_T> res_vec(dim, Scalar_T(0));
      for (set_value_t i = 1; i < dim; ++i)
      {
        if (l_vec[i] == Scalar_T(0)) continue;
        unsigned long uthis = unfolded_bits[i];
        // j is submask of i
        for (set_value_t j = i; ; j = (j - 1) & i)
        {
          if (j == 0) break;
          if (r_vec[j] != Scalar_T(0))
          {
            unsigned long urhs = unfolded_bits[j];
            unsigned long h = h_vals[j];
            unsigned long k_val = inverse_gray(uthis & h);
            unsigned long q_val = inverse_gray((uthis & urhs) >> -index_set_t::v_lo);
            unsigned long neg = k_val ^ q_val;
            Scalar_T sign = (neg & 1) ? Scalar_T(-1) : Scalar_T(1);
            res_vec[i ^ j] += sign * l_vec[i] * r_vec[j];
          }
          if (i == 0) break;
        }
        // j is strict superset of i
        for (set_value_t j = (i + 1) | i; j < dim; j = (j + 1) | i)
        {
          if (r_vec[j] != Scalar_T(0))
          {
            unsigned long urhs = unfolded_bits[j];
            unsigned long h = h_vals[j];
            unsigned long k_val = inverse_gray(uthis & h);
            unsigned long q_val = inverse_gray((uthis & urhs) >> -index_set_t::v_lo);
            unsigned long neg = k_val ^ q_val;
            Scalar_T sign = (neg & 1) ? Scalar_T(-1) : Scalar_T(1);
            res_vec[i ^ j] += sign * l_vec[i] * r_vec[j];
          }
        }
      }
      auto result = multivector_t(_GLUCAT_HASH_SIZE_T(dim));
      for (set_value_t k = 0; k < dim; ++k)
        if (res_vec[k] != Scalar_T(0))
          result.insert(term_t(index_set_t(k, our_frame, true), res_vec[k]));
      return result;
    }
    else
    {
      using matrix_multi_t = typename multivector_t::matrix_multi_t;
      return multivector_t(matrix_multi_t(lhs, our_frame, true) &
                           matrix_multi_t(rhs, our_frame, true));
    }
  }

  /*
   * @brief Inner product
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param rhs Right hand side
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator&= (const multivector_t& rhs) -> multivector_t&
  { return *this = *this & rhs; }

  /*
   * @brief Left contraction
   * @details
   *
   * Usage example:
   * Location: test00/peg00.h:205
   *
   * @code
   *
   * rhs = (a_1 % b) + (a_1 ^ b);
   * @endcode
   *
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Result
   *
   * @par Example:
   * @code
   * clifford<>("{1,2}") % clifford<>("{1}"); // Returns -{2}
   * @endcode
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  operator% (const framed_multi<Scalar_T,LO,HI,Tune_P>& lhs, const framed_multi<Scalar_T,LO,HI,Tune_P>& rhs) -> framed_multi<Scalar_T,LO,HI,Tune_P>
  {
    // Reference: Leo Dorst, "Honing geometric algebra for its use in the computer sciences",
    // in Geometric Computing with Clifford Algebras, ed. G. Sommer,
    // Springer 2001, Chapter 6, pp. 127-152.
    // http://staff.science.uva.nl/~leo/clifford/index.html

    if (lhs.empty() || rhs.empty())
      return Scalar_T(0);

    using multivector_t = framed_multi<Scalar_T,LO,HI,Tune_P>;
    using index_set_t = typename multivector_t::index_set_t;
    using term_t = typename multivector_t::term_t;

    const auto our_frame = lhs.frame() | rhs.frame();
    const auto algebra_dim = set_value_t(1) << our_frame.count();
    const double lhs_size = lhs.size();
    const double rhs_size = rhs.size();
    using Tuning_Values_P = typename Tune_P::tuning_values_p;
    const auto frm_count = our_frame.count();
    if (lhs_size * rhs_size <= double(algebra_dim))
    {
      auto result = multivector_t(
        _GLUCAT_HASH_SIZE_T(size_t(std::min(lhs_size * rhs_size, double(algebra_dim)))));
      for (auto& rhs_term : rhs)
      {
        const auto rhs_ist = rhs_term.first;
        for (auto& lhs_term : lhs)
        {
          const index_set_t lhs_ist = lhs_term.first;
          if ((lhs_ist | rhs_ist) == rhs_ist)
            result += term_t(lhs_term) * term_t(rhs_term);
        }
      }
      return result;
    }
    else if (frm_count < index_t(Tuning_Values_P::products_matrix_threshold))
    {
      const set_value_t dim = algebra_dim;
      std::vector<Scalar_T> l_vec(dim, Scalar_T(0));
      for (auto& term : lhs) l_vec[term.first.value_of_fold(our_frame)] = term.second;
      std::vector<Scalar_T> r_vec(dim, Scalar_T(0));
      for (auto& term : rhs) r_vec[term.first.value_of_fold(our_frame)] = term.second;

      std::vector<set_value_t> unfolded_bits(dim);
      std::vector<unsigned long> h_vals(dim);
      for (set_value_t k = 0; k < dim; ++k)
      {
        unfolded_bits[k] = index_set_t(k, our_frame, true).to_set_value();
        h_vals[k] = inverse_reversed_gray(unfolded_bits[k]);
      }

      std::vector<Scalar_T> res_vec(dim, Scalar_T(0));
      for (set_value_t i = 0; i < dim; ++i)
      {
        if (l_vec[i] == Scalar_T(0)) continue;
        unsigned long uthis = unfolded_bits[i];
        // j is superset of i
        for (set_value_t j = i; j < dim; j = (j + 1) | i)
        {
          if (r_vec[j] != Scalar_T(0))
          {
            unsigned long urhs = unfolded_bits[j];
            unsigned long h = h_vals[j];
            unsigned long k_val = inverse_gray(uthis & h);
            unsigned long q_val = inverse_gray((uthis & urhs) >> -index_set_t::v_lo);
            unsigned long neg = k_val ^ q_val;
            Scalar_T sign = (neg & 1) ? Scalar_T(-1) : Scalar_T(1);
            res_vec[i ^ j] += sign * l_vec[i] * r_vec[j];
          }
        }
      }
      auto result = multivector_t(_GLUCAT_HASH_SIZE_T(dim));
      for (set_value_t k = 0; k < dim; ++k)
        if (res_vec[k] != Scalar_T(0))
          result.insert(term_t(index_set_t(k, our_frame, true), res_vec[k]));
      return result;
    }
    else
    {
      using matrix_multi_t = typename multivector_t::matrix_multi_t;
      return multivector_t(matrix_multi_t(lhs, our_frame, true) %
                           matrix_multi_t(rhs, our_frame, true));
    }
  }

  /*
   * @brief Left contraction
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param rhs Right hand side
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator%= (const multivector_t& rhs) -> multivector_t&
  { return *this = *this % rhs; }

  /*
   * @brief Scalar product: [HS] (1.44) star(a, b) = scalar(a * b) = <ab>_0
   * @details
   *
   * Usage example:
   * Location: test00/peg00.h:260
   *
   * @code
   *
   * const scalar_t scalar_lhs = star(a, b);
   * @endcode
   *
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  star(const framed_multi<Scalar_T,LO,HI,Tune_P>& lhs, const framed_multi<Scalar_T,LO,HI,Tune_P>& rhs) -> Scalar_T
  {
    auto result = Scalar_T(0);
    const auto small_star_large = lhs.size() < rhs.size();
    const auto* smallp =
      small_star_large
      ? &lhs
      : &rhs;
    const auto* largep =
      small_star_large
      ? &rhs
      : &lhs;

    for (auto& small_term : *smallp)
    {
      const auto small_ist = small_term.first;
      const auto large_crd = (*largep)[small_ist];
      if (large_crd != Scalar_T(0))
        result += small_ist.sign_of_square() * small_term.second * large_crd;
    }
    return result;
  }
 
  /*
   * @brief Hestenes inner product: [H] (1.10) hstar(a, b) = scalar(reverse(a) * b) = <a†b>_0
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  hstar(const framed_multi<Scalar_T,LO,HI,Tune_P>& lhs, const framed_multi<Scalar_T,LO,HI,Tune_P>& rhs) -> Scalar_T
  {
    return scalar(lhs.reverse() * rhs);
  }

  /*
   * @brief Quotient of multivector and scalar
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param scr Value
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
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
      for (auto& this_term : *this)
        const_cast<Scalar_T&>(this_term.second) /= scr;
    return *this;
  }

  /*
   * @brief Geometric quotient
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Quotient
   *
   * @par Example:
   * @code
   * clifford<>("2{1}") / clifford<>("{1}"); // Returns 2
   * @endcode
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  operator/ (const framed_multi<Scalar_T,LO,HI,Tune_P>& lhs, const framed_multi<Scalar_T,LO,HI,Tune_P>& rhs) -> framed_multi<Scalar_T,LO,HI,Tune_P>
  {
    using multivector_t = framed_multi<Scalar_T,LO,HI,Tune_P>;
    using traits_t = numeric_traits<Scalar_T>;
    using matrix_multi_t = typename multivector_t::matrix_multi_t;

    if (rhs == Scalar_T(0))
      return traits_t::NaN();

    const auto our_frame = lhs.frame() | rhs.frame();
    return multivector_t(matrix_multi_t(lhs, our_frame, true) / matrix_multi_t(rhs, our_frame, true));
  }

  /*
   * @brief Geometric quotient
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param rhs Right hand side
   * @return Geometric quotient *this / rhs
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator/= (const multivector_t& rhs) -> multivector_t&
  { return *this = *this / rhs; }

  /*
   * @brief Transformation via twisted adjoint action
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return rhs transformed by lhs
   *
   * @par Example:
   * @code
   * clifford<>("{2}") | clifford<>("{1,2}"); // Returns -{2}
   * @endcode
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  operator| (const framed_multi<Scalar_T,LO,HI,Tune_P>& lhs, const framed_multi<Scalar_T,LO,HI,Tune_P>& rhs) -> framed_multi<Scalar_T,LO,HI,Tune_P>
  {
    using multivector_t = framed_multi<Scalar_T,LO,HI,Tune_P>;
    using matrix_multi_t = typename multivector_t::matrix_multi_t;

    return multivector_t(matrix_multi_t(rhs) * matrix_multi_t(lhs) / matrix_multi_t(rhs.involute()));
  }

  /*
   * @brief Transformation via twisted adjoint action
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param rhs Right hand side
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator|= (const multivector_t& rhs) -> multivector_t&
  { return *this = *this | rhs; }

  /*
   * @brief Clifford multiplicative inverse
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  inv() const -> multivector_t
  {
    auto result = matrix_multi_t(Scalar_T(1), this->frame());
    return multivector_t(result /= matrix_multi_t(*this));
  }

  /*
   * @brief Move assignment
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param other Other matrix
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator= (framed_multi&& other) noexcept(std::is_nothrow_move_assignable_v<Scalar_T>) -> multivector_t&
  {
    map_t::operator=(std::move(other));
    return *this;
  }

  /*
   * @brief Integer power of multivector: *this to the m
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param m Matrix
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  pow(int m) const -> multivector_t
  { return glucat::pow(*this, m); }

  /*
   * @brief Outer product power of multivector
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param m Matrix
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  outer_pow(int m) const -> multivector_t
  {
    if (m < 0)
      throw error_t("outer_pow(int): negative exponent");
    auto result = multivector_t(Scalar_T(1));
    auto a = *this;
    for (;
        m != 0;
        m >>= 1, a = a ^ a)
      if (m & 1)
        result ^= a;
    return result;
  }

  /*
   * @brief Frame of multivector: union of index sets of terms
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  frame() const -> index_set_t
  {
    auto result = index_set_t();
    for (auto& this_term : *this)
      result |= this_term.first;
    return result;
  }

  /*
   * @brief Grade of multivector: maximum of the grades of each term
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Grade
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  grade() const -> index_t
  {
    auto result = index_t(0);
    for (auto& this_term : *this)
      result = std::max(result, this_term.first.count());
    return result;
  }

  /*
   * @brief Subscripting: map from index set to scalar coordinate
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param ist Value
   * @return Element reference
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator[] (const index_set_t ist) const -> Scalar_T
  {
    const auto& this_it = this->find(ist);
    if (this_it == this->end())
      return Scalar_T(0);
    else
      return this_it->second;
  }

  /*
   * @brief Grading: part where each term is a grade-vector
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Element
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator() (index_t grade) const -> multivector_t
  {
    if ((grade < 0) || (grade > HI-LO))
      return Scalar_T(0);
    else
    {
      auto result = multivector_t();
      for (auto& this_term : *this)
        if (this_term.first.count() == grade)
          result += this_term;
      return result;
    }
  }

  /*
   * @brief Scalar part
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Scalar part
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  scalar() const -> Scalar_T
  { return (*this)[index_set_t()]; }

  /*
   * @brief Pure part
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Pure part
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  pure() const -> multivector_t
  { return *this - this->scalar(); }

  /*
   * @brief Even part, sum of the even grade terms
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Even part
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  even() const -> multivector_t
  { // even part of x, sum of the pure(count) with even count
    auto result = multivector_t();
    for (auto& this_term : *this)
      if ((this_term.first.count() % 2) == 0)
        result.insert(this_term);
    return result;
  }

  /*
   * @brief Odd part, sum of the odd grade terms
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Odd part
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  odd() const -> multivector_t
  { // even part of x, sum of the pure(count) with even count
    auto result = multivector_t();
    for (auto& this_term : *this)
      if ((this_term.first.count() % 2) == 1)
        result.insert(this_term);
    return result;
  }

  /*
   * @brief Vector part of multivector, as a vector_t
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Vector part
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  vector_part() const -> vector_t
  { return this->vector_part(this->frame(), true); }

  /*
   * @brief Vector part of multivector, as a vector_t
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param frm Frame
   * @param prechecked Bool: true if i has already been checked?
   * @return True if successful or condition met
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  vector_part(const index_set_t frm, const bool prechecked) const -> vector_t
  {
    if (!prechecked && (this->frame() | frm) != frm)
      throw error_t("vector_part(frm): value is outside of requested frame");
    auto result = vector_t();
    result.reserve(frm.count());
    const auto frm_end = frm.max()+1;
    for (auto
        idx  = frm.min();
        idx != frm_end;
        ++idx)
      // Frame may contain indices which do not correspond to a grade 1 term but
      // frame cannot omit any index corresponding to a grade 1 term
      if (frm[idx])
        result.push_back((*this)[index_set_t(idx)]);
    return result;
  }

  /*
   * @brief Main involution, each {i} is replaced by -{i} in each term
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Main involution
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  involute() const -> multivector_t
  {
    auto result = *this;
    for (auto& result_term : result)
    { // for a k-vector u, involute(u) == (-1)^k * u
      if ((result_term.first.count() % 2) == 1)
        const_cast<Scalar_T&>(result_term.second) *= Scalar_T(-1);
    }
    return result;
  }

  /*
   * @brief Reversion, order of {i} is reversed in each term
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Reverse
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  reverse() const -> multivector_t
  {
    auto result = *this;
    for (auto& result_term : result)
      // For a k-vector u, reverse(u) = { -u, k == 2,3 (mod 4)
      //                                {  u, k == 0,1 (mod 4)
      switch (result_term.first.count() % 4)
      {
      case 2:
      case 3:
        const_cast<Scalar_T&>(result_term.second) *= Scalar_T(-1);
        break;
      default:
        break;
      }
    return result;
  }

  /*
   * @brief Conjugation, conj == reverse o involute == involute o reverse
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Clifford conjugate
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  conj() const -> multivector_t
  {
    auto result = *this;
    // Use explicit iterator to avoid deduction issues with robin_map
    using map_t = typename framed_multi<Scalar_T,LO,HI,Tune_P>::map_t;
    for (typename map_t::iterator it = result.begin(); it != result.end(); ++it)
    {
      auto& result_term = *it;
      // For a k-vector u, conj(u) = { -u, k == 1,2 (mod 4)
      //                             {  u, k == 0,3 (mod 4)
      switch (result_term.first.count() % 4)
      {
      case 1:
      case 2:
        const_cast<Scalar_T&>(result_term.second) *= Scalar_T(-1);
        break;
      default:
        break;
      }
    }
    return result;
  }

  /*
   * @brief Quadratic form := scalar part of rev(x)*x
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Quadratic form
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  quad() const -> Scalar_T
  {
    // scalar(conj(x)*x) = 2*quad(even(x)) - quad(x)
    // ref: old clical: quadfunction(p:pter):pterm in file compmod.pas
    auto result = Scalar_T(0);
    for (auto& this_term : *this)
    {
      const auto sign =
        (this_term.first.count_neg() % 2)
        ? -Scalar_T(1)
        :  Scalar_T(1);
      result += sign * (this_term.second) * (this_term.second);
    }
    return result;
  }

  /*
   * @brief Norm squared := sum of norm squared of coordinates
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Norm
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  norm() const -> Scalar_T
  {
    using traits_t = numeric_traits<Scalar_T>;

    auto result = Scalar_T(0);
    for (auto& this_term : *this)
    {
      const auto abs_crd = traits_t::abs(this_term.second);
      result +=  abs_crd * abs_crd;
    }
    return result;
  }

  /*
   * @brief Maximum of absolute values of components of multivector: multivector infinity norm
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Maximum absolute value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  max_abs() const -> Scalar_T
  {
    using traits_t = numeric_traits<Scalar_T>;

    auto result = Scalar_T(0);
    for (auto& this_term : *this)
    {
      const auto abs_crd = traits_t::abs(this_term.second);
      if (abs_crd > result)
        result = abs_crd;
    }
    return result;
  }

  /*
   * @brief Random multivector within a frame
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param frm Frame
   * @param fill Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  random(const index_set_t frm, Scalar_T fill) -> multivector_t
  {
    using multivector_t = framed_multi<Scalar_T,LO,HI,Tune_P>;
    using index_set_t = typename multivector_t::index_set_t;
    using term_t = typename multivector_t::term_t;

    using random_generator_t = random_generator<Scalar_T>;
    auto& generator = random_generator_t::generator();

    fill =
      (fill < Scalar_T(0))
      ? Scalar_T(0)
      : (fill > Scalar_T(1))
        ? Scalar_T(1)
        : fill;
    const auto algebra_dim = set_value_t(1) << frm.count();
    using traits_t = numeric_traits<Scalar_T>;
    const auto mean_abs = traits_t::sqrt(Scalar_T(double(algebra_dim)));
    auto result = multivector_t();
    for (auto
        stv = set_value_t(0);
        stv != algebra_dim;
        ++stv)
      if (generator.uniform() < fill)
      {
        const auto& result_crd = generator.normal() / mean_abs;
        result.insert(term_t(index_set_t(stv, frm, true), result_crd));
      }
    return result;
  }

  /*
   * @brief Write multivector to output
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param msg Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  void
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  write(const std::string& msg) const
  { std::cout << msg << std::endl << "  " << (*this) << std::endl; }

  /*
   * @brief Write multivector to file
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param ofile Value
   * @param msg Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  void
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  write(std::ofstream& ofile, const std::string& msg) const
  {
    if (!ofile)
      throw error_t("write(ofile,msg): cannot write to output file");
    ofile << msg << std::endl << "  " << (*this) << std::endl;
  }

  /*
   * @brief Sorted range for use with output
   * @details
   * @tparam Map_T
   * @tparam Sorted_Map_T
   */
  template< typename Map_T,typename Sorted_Map_T >
  class sorted_range
  {
  public:
    using map_t = Map_T;
    using sorted_map_t = Sorted_Map_T;
    using sorted_iterator = typename Sorted_Map_T::const_iterator;

    sorted_range (Sorted_Map_T &sorted_val, const Map_T& val)
    {
      for (auto& val_term : val)
        sorted_val.insert(val_term);
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

  /*
   * @brief Write multivector to output
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param os Output stream
   * @param val Value
   * @return Output stream
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  operator<< (std::ostream& os, const framed_multi<Scalar_T,LO,HI,Tune_P>& val) -> std::ostream&
  {
    using limits_t = std::numeric_limits<Scalar_T>;
    if (val.empty())
      os << 0;
    else if (val.isnan())
      os << limits_t::quiet_NaN();
    else if (val.isinf())
    {
      const Scalar_T& inf = limits_t::infinity();
      os << (scalar(val) < 0.0 ? -inf : inf);
    }
    else
    {
      using traits_t = numeric_traits<Scalar_T>;
      using multivector_t = framed_multi<Scalar_T,LO,HI,Tune_P>;
      Scalar_T truncation;
      switch (os.flags() & std::ios::floatfield)
      {
        case std::ios_base::scientific:
          truncation = Scalar_T(1) / traits_t::pow(Scalar_T(10), int(os.precision()) + 1);
          break;
        case std::ios_base::fixed:
          truncation = Scalar_T(1) / (traits_t::pow(Scalar_T(10), int(os.precision())) * val.max_abs());
          break;
        case std::ios_base::fixed | std::ios_base::scientific:
          truncation = multivector_t::default_truncation;
          break;
        default:
          truncation = Scalar_T(1) / traits_t::pow(Scalar_T(10), int(os.precision()));
          break;
      }
      auto truncated_val = val.truncated(truncation);
      if (truncated_val.empty())
        os << 0;
      else
      {
        using map_t = typename multivector_t::map_t;
        using sorted_map_t = typename multivector_t::sorted_map_t;
        auto sorted_val = sorted_map_t();
        const auto sorted_val_range = sorted_range< map_t, sorted_map_t >(sorted_val, truncated_val);
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
    }
    return os;
  }

  /*
   * @brief Write term to output
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @param os Output stream
   * @param term Value
   * @return Output stream
   */
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator<< (std::ostream& os, const std::pair< const index_set<LO,HI>, Scalar_T >& term) -> std::ostream&
  {
    const auto second_as_double = numeric_traits<Scalar_T>::to_double(term.second);
    const auto use_double =
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
        auto tol = std::pow(10.0,-os.precision());
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

  /*
   * @brief Read multivector from input
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param s Value
   * @param val Value
   * @return Input stream
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  operator>> (std::istream& s, framed_multi<Scalar_T,LO,HI,Tune_P> & val) -> std::istream&
  { // Input looks like 1.0-2.0{1,2}+3.2{3,4}.
    using multivector_t = framed_multi<Scalar_T,LO,HI,Tune_P>;
    // Parsing variables.
    auto local_val = multivector_t();
    auto c = 0;
    // Parsing control variables.
    auto negative = false;
    auto expect_term = true;
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
      auto coordinate = Scalar_T(1);
      // Default index set is empty.
      auto ist = index_set<LO,HI>();
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

  /*
   * @brief Number of terms
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return Number of terms
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  nbr_terms () const -> size_type
  { return this->size(); }

  /*
   * @brief Insert a term into a multivector.
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param term The term to insert.
   * @return Reference to this.
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  operator+= (const term_t& term) -> multivector_t&
  { // Do not insert terms with 0 coordinate
    if (term.second != Scalar_T(0))
    {
      const auto& this_it = this->find(term.first);
      if (this_it == this->end())
        this->insert(term);
      else if (this_it->second + term.second == Scalar_T(0))
        // Erase term if resulting coordinate is 0
        this->erase(this_it);
      else
        const_cast<Scalar_T&>(this_it->second) += term.second;
    }
    return *this;
  }

  /*
   * @brief Check if a multivector contains any IEEE Inf values
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return True if condition met
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  isinf() const -> bool
  {
    using traits_t = numeric_traits<Scalar_T>;

    if (std::numeric_limits<Scalar_T>::has_infinity)
      for (auto& this_term : *this)
        if (traits_t::isInf(this_term.second))
          return true;
    return false;
  }

  /*
   * @brief Check if a multivector contains any IEEE NaN values
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @return True if successful or condition met
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  isnan() const -> bool
  {
    using traits_t = numeric_traits<Scalar_T>;

    if (std::numeric_limits<Scalar_T>::has_quiet_NaN)
      for (auto& this_term : *this)
        if (traits_t::isNaN(this_term.second))
          return true;
    return false;
  }

  /*
   * @brief Remove all terms with relative size smaller than limit
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param limit Truncation limit
   * @return Truncated mulivector
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  truncated(const Scalar_T& limit) const -> multivector_t
  {
    using traits_t = numeric_traits<Scalar_T>;

    if (this->isnan() || this->isinf())
      return *this;
    const auto truncation = traits_t::abs(limit);
    const auto top = max_abs();
    auto result = multivector_t();
    if (top != Scalar_T(0))
      for (auto& this_term : *this)
        if (traits_t::abs(this_term.second) > top * truncation)
          result.insert(this_term);
    return result;
  }

  /*
   * @brief Subalgebra isomorphism: fold each term within the given frame
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param frm Frame
   * @return Multivector folded from frame
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  fold(const index_set_t frm) const -> multivector_t
  {
    if (frm.is_contiguous())
      return *this;
    else
    {
      auto result = multivector_t();
      for (auto& this_term : *this)
        result.insert(term_t(this_term.first.fold(frm), this_term.second));
      return result;
    }
  }

  /*
   * @brief Subalgebra isomorphism: unfold each term within the given frame
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param frm Frame
   * @return Multivector unfolded into frame
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  unfold(const index_set_t frm) const -> multivector_t
  {
    if (frm.is_contiguous())
      return *this;
    else
    {
      auto result = multivector_t();
      for (auto& this_term : *this)
        result.insert(term_t(this_term.first.unfold(frm), this_term.second));
      return result;
    }
  }

  /*
   * @brief Subalgebra isomorphism: R_{p,q} to R_{p-4,q+4}
   * @details
   * @param p Positive index in R_{p,q}
   * @param q Negative index in R_{p,q}
   * @return Resulting multivector
   */
  // Reference: [L] 16.4 Periodicity of 8, p216
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  centre_pm4_qp4(index_t& p, index_t& q) -> multivector_t&
  {
    // We add 4 to q by subtracting 4 from p
    if (q+4 > -LO)
      throw error_t("centre_pm4_qp4(p,q): LO is too high to represent this value");
    if (this->frame().max() > p-4)
    {
      using index_pair_t = typename index_set_t::index_pair_t;
      const auto pm3210 = index_set_t(index_pair_t(p-3,p), true);
      const auto qm4321 = index_set_t(index_pair_t(-q-4,-q-1), true);
      const auto& tqm4321 = term_t(qm4321, Scalar_T(1));
      auto result = multivector_t();
      for (auto& this_term : *this)
      {
        const auto ist = this_term.first;
        if (ist.max() > p-4)
        {
          auto var_term = var_term_t();
          for (auto
              n = index_t(0);
              n != index_t(4);
              ++n)
            if (ist[n+p-3])
              var_term *= term_t(index_set_t(n-q-4), Scalar_T(1)) * tqm4321;
          // Mask out {p-3}..{p}
          result.insert(term_t(ist & ~pm3210, this_term.second) *
                        term_t(var_term.first, var_term.second));
        }
        else
          result.insert(this_term);
      }
      *this = result;
    }
    p -=4; q += 4;
    return *this;
  }

  /*
   * @brief Subalgebra isomorphism: R_{p,q} to R_{p+4,q-4}
   * @details
   * @param p Positive index in R_{p,q}
   * @param q Negative index in R_{p,q}
   * @return Resulting multivector
   */
  // Reference: [L] 16.4 Periodicity of 8, p216
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  centre_pp4_qm4(index_t& p, index_t& q) -> multivector_t&
  {
    // We add 4 to p by subtracting 4 from q
    if (p+4 > HI)
      throw error_t("centre_pp4_qm4(p,q): HI is too low to represent this value");
    if (this->frame().min() < -q+4)
    {
      using index_pair_t = typename index_set_t::index_pair_t;
      const auto qp0123 = index_set_t(index_pair_t(-q,-q+3), true);
      const auto pp1234 = index_set_t(index_pair_t(p+1,p+4), true);
      const auto& tpp1234 = term_t(pp1234, Scalar_T(1));
      auto result = multivector_t();
      for (auto& this_term : *this)
      {
        index_set_t ist = this_term.first;
        if (ist.min() < -q+4)
        {
          auto var_term = var_term_t();
          for (auto
              n = index_t(0);
              n != index_t(4);
              ++n)
            if (ist[n-q])
              var_term *= term_t(index_set_t(n+p+1), Scalar_T(1)) * tpp1234;
          // Mask out {-q}..{-q+3}
          result.insert(term_t(var_term.first, var_term.second) *
                        term_t(ist & ~qp0123, this_term.second));
        }
        else
          result.insert(this_term);
      }
      *this = result;
    }
    p +=4; q -= 4;
    return *this;
  }

  /*
   * @brief Subalgebra isomorphism: R_{p,q} to R_{q+1,p-1}
   * @details
   * @param p Positive index in R_{p,q}
   * @param q Negative index in R_{p,q}
   * @return Resulting multivector
   */
  // Reference: [P] Proposition 15.20, p 131
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  centre_qp1_pm1(index_t& p, index_t& q) -> multivector_t&
  {
    if (q+1 > HI)
      throw error_t("centre_qp1_pm1(p,q): HI is too low to represent this value");
    if (p-1 > -LO)
      throw error_t("centre_qp1_pm1(p,q): LO is too high to represent this value");
    const auto qp1 = index_set_t(q+1);
    const auto& tqp1 = term_t(qp1, Scalar_T(1));
    auto result = multivector_t();
    for (auto& this_term : *this)
    {
      const auto ist = this_term.first;
      auto var_term = var_term_t(index_set_t(), this_term.second);
      for (auto
          n = -q;
          n != p;
          ++n)
        if (n != 0 && ist[n])
          var_term *= term_t(index_set_t(-n) | qp1, Scalar_T(1));
      if (p != 0 && ist[p])
        var_term *= tqp1;
      result.insert(term_t(var_term.first, var_term.second));
    }
    index_t orig_p = p;
    p = q+1;
    q = orig_p-1;
    return *this = result;
  }

  /*
   * @brief Divide multivector *this into quotient with terms divisible by index set, and remainder
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param ist Index set
   * @return Pair containing quotient and remainder
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  divide(const index_set_t ist) const -> framed_pair_t
  {
    auto quo = multivector_t();
    auto rem = multivector_t();
    for (auto& this_term : *this)
      if ((this_term.first | ist) == this_term.first)
        quo.insert(term_t(this_term.first ^ ist, this_term.second));
      else
        rem.insert(this_term);
    return framed_pair_t(quo, rem);
  }

  /*
   * @brief Generalized FFT from multivector_t to matrix_t
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param level Recursion level
   * @param odd Bool: If true, take the odd part, otherwise take the even part
   * @return The matrix representing the multivector *this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  fast(const index_t level, const bool odd) const -> matrix_t
  {
    // Assume val is already folded and centred

    if (this->empty())
    {
      const auto dim = matrix_index_t(1) << level;
      auto result = matrix_t(dim, dim);
      result.clear();
      return result;
    }
    if (level == 0)
      return matrix::unit<matrix_t>(1) * this->scalar();

    auto I = matrix::unit<matrix_t>(2);
    auto J = matrix_t(2,2);
    J.zeros(); // Ensure zeroed
    J(0,1)  = Scalar_T(-1);
    J(1,0)  = Scalar_T( 1);
    auto K = J;
    K(0,1)  = Scalar_T( 1);
    auto JK = I;
    JK(0,0) = Scalar_T(-1);

    const auto ist_mn = index_set_t(-level);
    const auto ist_pn = index_set_t(level);
    if (level == 1)
    {
      const auto& val_mn = (*this)[ist_mn];
      const auto& val_pn = (*this)[ist_pn];
      const auto& val_scalar = this->scalar();
      const auto& val_mnpn = (*this)[ist_mn ^ ist_pn];

      if (odd)
        return matrix_t(J * val_mn + K * val_pn);
      else
        return matrix_t(I * val_scalar + JK * val_mnpn);
    }
    else
    {
      const auto& pair_mn = this->divide(ist_mn);
      const auto& quo_mn = pair_mn.first;
      const auto& rem_mn = pair_mn.second;
      const auto& pair_quo_mnpn = quo_mn.divide(ist_pn);
      const auto& val_mnpn = pair_quo_mnpn.first;
      const auto& val_mn   = pair_quo_mnpn.second;
      const auto& pair_rem_mnpn = rem_mn.divide(ist_pn);
      const auto& val_pn   = pair_rem_mnpn.first;
      const auto& val_1    = pair_rem_mnpn.second;
      if (odd)
        return - JK.kron(val_1.fast   (level-1, 1))
               + I.kron(val_mnpn.fast(level-1, 1))
               + J.kron(val_mn.fast  (level-1, 0))
               + K.kron(val_pn.fast  (level-1, 0));
      else
        return   I.kron(val_1.fast   (level-1, 0))
               + JK.kron(val_mnpn.fast(level-1, 0))
               + K.kron(val_mn.fast  (level-1, 1))
               - J.kron(val_pn.fast  (level-1, 1));
    }
  }

  /*
   * @brief Use generalized FFT to construct a matrix_multi_t
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param frm Frame
   * @return The matrix_multi_t value representing the multivector *this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  template< typename Other_Scalar_T, typename Other_Tune_P >
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  fast_matrix_multi(const index_set_t frm) const -> matrix_multi<Other_Scalar_T,LO,HI,Other_Tune_P>
  {
    // Fold val
    auto val = this->fold(frm);
    auto p = frm.count_pos();
    auto q = frm.count_neg();
    const auto bott_offset = gen::offset_to_super[pos_mod(p - q, 8)];
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
    const auto& ev_val = val.even();
    const auto& od_val = val.odd();

    auto ev_res = ev_val.fast(level, 0);
    auto od_res = od_val.fast(level, 1);

    return matrix_multi<Other_Scalar_T,LO,HI,Other_Tune_P>(ev_res + od_res, frm);
  }

  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  auto
  framed_multi<Scalar_T,LO,HI,Tune_P>::
  fast_framed_multi() const -> multivector_t
  { return *this; }

  /*
   * @brief Coordinate of product of terms
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
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @param lhs Left hand side term
   * @param rhs Right hand side term
   * @return The coordinate of lhs * rhs as multivectors
   */
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  static
  auto
  crd_of_mult(const std::pair<const index_set<LO,HI>, Scalar_T>& lhs,
              const std::pair<const index_set<LO,HI>, Scalar_T>& rhs) -> Scalar_T
  { return lhs.first.sign_of_mult(rhs.first) * lhs.second * rhs.second; }

  /*
   * @brief Product of two terms
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
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @param lhs Left hand side term
   * @param rhs Right hand side term
   * @return Clifford product of lhs and rhs as multivectors
   */
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  auto
  operator* (const std::pair<const index_set<LO,HI>, Scalar_T>& lhs,
             const std::pair<const index_set<LO,HI>, Scalar_T>& rhs) -> std::pair<const index_set<LO,HI>, Scalar_T>
  {
    using term_t = std::pair<const index_set<LO,HI>, Scalar_T>;
    return term_t(lhs.first ^ rhs.first, crd_of_mult(lhs, rhs));
  }

  /*
   * @brief Square root of multivector with specified complexifier
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param val Multivector
   * @param i Complexifier -- commuting sqrt of -1
   * @param prechecked Bool: true if i has already been checked?
   * @return Square root of val
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  sqrt(const framed_multi<Scalar_T,LO,HI,Tune_P>& val, const framed_multi<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> framed_multi<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    check_complex(val, i, prechecked);

    const auto realval = val.scalar();
    if (val == realval)
    {
      if (realval < Scalar_T(0))
        return i * traits_t::sqrt(-realval);
      else
        return traits_t::sqrt(realval);
    }
    using multivector_t = framed_multi<Scalar_T,LO,HI,Tune_P>;
    using matrix_multi_t = typename multivector_t::matrix_multi_t;
    return multivector_t(sqrt(matrix_multi_t(val), matrix_multi_t(i), prechecked));
  }

  /*
   * @brief Exponential of multivector
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param val Multivector
   * @return Exponential of val
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  exp(const framed_multi<Scalar_T,LO,HI,Tune_P>& val) -> framed_multi<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    const auto s = scalar(val);
    if (val == s)
      return traits_t::exp(s);

    using Tuning_Values_P = typename Tune_P::tuning_values_p;
    const double size = val.size();
    const auto frm_count = val.frame().count();
    const auto algebra_dim = set_value_t(1) << frm_count;

    using multivector_t = framed_multi<Scalar_T,LO,HI,Tune_P>;

    if( (size * size <= double(algebra_dim)) || (frm_count < index_t(Tuning_Values_P::mult_matrix_threshold)))
    {
      using tune_same_p = typename Tune_P::tuning_same_p;
      const precision_t function_precision = Tune_P::function_precision;

      if constexpr (function_precision == precision_demoted)
      {
        using demoted_scalar_t = typename traits_t::demoted::type;
        using demoted_multivector_t = framed_multi<demoted_scalar_t,LO,HI,tune_same_p>;

        const auto& demoted_val = demoted_multivector_t(val);
        return multivector_t(clifford_exp(demoted_val));
      }
      else if constexpr (function_precision == precision_promoted)
      {
        using promoted_scalar_t = typename traits_t::promoted::type;
        using promoted_multivector_t = framed_multi<promoted_scalar_t,LO,HI,tune_same_p>;

        const auto& promoted_val = promoted_multivector_t(val);
        return multivector_t(clifford_exp(promoted_val));
      }
      else
        return clifford_exp(val);
    }
    else
    {
      using matrix_multi_t = matrix_multi<Scalar_T,LO,HI,Tune_P>;
      return multivector_t(exp(matrix_multi_t(val)));
    }
  }

  /*
   * @brief Natural logarithm of multivector with specified complexifier
   * @details
   * @tparam Scalar_T Scalar type
   * @tparam LO Low index limit
   * @tparam HI High index limit
   * @tparam Tune_P Tuning policy
   * @param val Multivector
   * @param i Complexifier -- commuting sqrt of -1
   * @param prechecked Bool: true if i has already been checked?
   * @return Logarithm of val
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  log(const framed_multi<Scalar_T,LO,HI,Tune_P>& val, const framed_multi<Scalar_T,LO,HI,Tune_P>& i, bool prechecked) -> framed_multi<Scalar_T,LO,HI,Tune_P>
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();

    check_complex(val, i, prechecked);

    const auto realval = val.scalar();
    if (val == realval)
    {
      if (realval < Scalar_T(0))
        return i * traits_t::pi() + traits_t::log(-realval);
      else
        return traits_t::log(realval);
    }
    using multivector_t = framed_multi<Scalar_T,LO,HI,Tune_P>;
    using matrix_multi_t = typename multivector_t::matrix_multi_t;
    return multivector_t(log(matrix_multi_t(val), matrix_multi_t(i), prechecked));
  }
}
#ifdef GLUCAT_DOCTEST
#include <iostream>
#include <iomanip>
#include <numbers>
#include <filesystem>
#include <system_error>
#include <chrono>

TEST_CASE("framed_multi<Scalar_T, LO, HI, Tune_P>") {
  using namespace glucat;
  using fm_t = glucat::framed_multi<double, -8, 8>;
  using T = double;
  using is_t = fm_t::index_set_t;

  SUBCASE("Metadata") {
    CHECK(fm_t::classname() == "framed_multi");
  }

  SUBCASE("Constructor and string representation") {
    fm_t f1(2.0);
    std::ostringstream oss1;
    oss1 << f1;
    CHECK(oss1.str() == "2");

    fm_t f2("2{1,2,3}");
    std::ostringstream oss2;
    oss2 << f2;
    CHECK(oss2.str() == "2{1,2,3}");

    fm_t f3("-{1}");
    std::ostringstream oss3;
    oss3 << f3;
    CHECK(oss3.str() == "-{1}");
  }

  SUBCASE("Geometric operations") {
    fm_t e1("{1}");
    fm_t e2("{2}");
    fm_t e3("{3}");
    CHECK((e1 * e2) == fm_t("{1,2}"));
    CHECK((e2 * e1) == fm_t("-{1,2}"));
    CHECK((e1 * e1) == fm_t(1.0));

    // HS (1.21a): (a_r * b_s)(|r-s|) == a_r & b_s
    CHECK(scalar(e1 * e2) == 0.0);
    CHECK(scalar(e1 * e1) == 1.0);

    // HS (1.25a): (a ^ b) ^ c == a ^ (b ^ c)
    CHECK(((e1 ^ e2) ^ e3) == (e1 ^ (e2 ^ e3)));

    // HS (1.31): a_1 * b == (a_1 & b) + (a_1 ^ b)
    fm_t b = fm_t("1+{1}+{2}+{1,2}");
    CHECK((e1 * b) == ((e1 & b) + (e1 ^ b)));

    // HS (1.44): star(a, b) == scalar(a * b)
    CHECK(star(e1, e2) == scalar(e1 * e2));
    CHECK(star(e1, e1) == scalar(e1 * e1));
  }

  SUBCASE("Arithmetic and approximate equality") {
    fm_t m1(1.0);
    fm_t m2("{1}");
    CHECK((m1 + m2) == fm_t("1+{1}"));
    CHECK((m1 - m2) == fm_t("1-{1}"));
    CHECK((m1 * 2.0) == fm_t(2.0));
    CHECK((2.0 * m1) == fm_t(2.0));
    CHECK(approx_equal(m1, fm_t(1.0 + 1e-15)));
    CHECK(m1 != m2);
    CHECK(m1 != 0.0);
    CHECK(0.0 != m1);
  }

  SUBCASE("Transcendental functions") {
    fm_t x("{1,2}");
    const double pi = std::numbers::pi;

    // exp and log
    fm_t e_x = exp(x * (pi/4.0));
    // exp({1,2}*pi/4) = cos(pi/4) + {1,2}*sin(pi/4) = (1 + {1,2})/sqrt(2)
    CHECK(approx_equal(e_x, (fm_t(1.0) + fm_t("{1,2}")) / std::sqrt(2.0)));
    CHECK(approx_equal(exp(log(e_x)), e_x));

    // sin and cos
    fm_t s_x = sin(x);
    fm_t c_x = cos(x);
    // sin^2 + cos^2 = 1
    CHECK(approx_equal(s_x*s_x + c_x*c_x, fm_t(1.0)));
    CHECK(approx_equal(cos(acos(fm_t(0.5))), fm_t(0.5)));

    // sinh and cosh
    fm_t sh_x = sinh(x);
    fm_t ch_x = cosh(x);
    // cosh^2 - sinh^2 = 1
    CHECK(approx_equal(ch_x*ch_x - sh_x*sh_x, fm_t(1.0)));
    CHECK(approx_equal(cosh(acosh(fm_t(2.0))), fm_t(2.0)));

    // tan and tanh
    CHECK(approx_equal(tan(atan(fm_t(0.5))), fm_t(0.5)));
    CHECK(approx_equal(tanh(atanh(fm_t(0.5))), fm_t(0.5)));

    // peg11.h identities
    fm_t A = fm_t("0.5{1}+0.5{1,2}");
    CHECK(approx_equal(exp(A) * exp(-A), fm_t(1.0)));
    CHECK(approx_equal(cosh(A) + sinh(A), exp(A)));
    CHECK(approx_equal(cos(A) + complexifier(A)*sin(A), exp(complexifier(A)*A)));
    CHECK(approx_equal(cos(A)*tan(A), sin(A)));
    CHECK(approx_equal(cosh(A)*tanh(A), sinh(A)));
    CHECK(approx_equal(sqrt(fm_t(4.0)), fm_t(2.0)));
  }

  SUBCASE("Adversarial and edge cases") {
    // NaN handling
    fm_t n(glucat::numeric_traits<double>::NaN());
    CHECK(n.isnan());
    CHECK(exp(n).isnan());
    CHECK(log(n).isnan());
    CHECK(log(n, fm_t(glucat::numeric_traits<double>::NaN()), false).isnan());

    // Purely scalar multivector
    fm_t s(T(2.0), is_t());
    CHECK(approx_equal(exp(s), fm_t(std::exp(2.0))));
    CHECK(approx_equal(log(s), fm_t(std::log(2.0))));
    CHECK(approx_equal(sqrt(s), fm_t(std::sqrt(2.0))));
    // For log of scalar, we need a valid complexifier
    fm_t i("{-1}");
    CHECK(approx_equal(log(s, i, false), fm_t(std::log(2.0))));

    // Zero
    fm_t zero(T(0.0), is_t());
    // Log of zero returns NaN in GluCat
    CHECK(log(zero).isnan());

    // Large multivector to trigger matrix-based exp
    fm_t large = fm_t("1+{1}+{2}+{3}+{4}+{5}+{1,2}+{1,3}+{1,4}+{1,5}+{2,3}+{2,4}+{2,5}+{3,4}+{3,5}+{4,5}");
    CHECK(approx_equal(exp(large), exp(glucat::matrix_multi<double, -8, 8>(large))));

    // Negative log for framed_multi
    fm_t neg(T(-1.0), is_t());
    CHECK(approx_equal(log(neg, fm_t("{-1}"), false), fm_t("{-1}") * T(std::numbers::pi)));
  }
  SUBCASE("Transcendental identities (random)") {
    using namespace glucat;
    using index_set_t = fm_t::index_set_t; using is_t = fm_t::index_set_t; index_set_t frm = index_set_t();
    const double fill = 0.5;
    for (index_t i = 1; i <= 7; ++i) {
      frm |= index_set_t(i);
      frm |= index_set_t(-i);

      fm_t a = fm_t::random(frm, fill);

      // exp(a) * exp(-a) == 1
      CHECK(approx_equal(exp(a) * exp(-a), fm_t(T(1.0), is_t())));

      // cosh(a) + sinh(a) == exp(a)
      CHECK(approx_equal(cosh(a) + sinh(a), exp(a)));

      // sqrt(a) * sqrt(a) == a
      // Note: sqrt might not always return a value that squares back to a due to branch cuts,
      // but for many random a it should work or be close.
      fm_t s = sqrt(a);
      if (!s.isnan() && !s.isinf())
        CHECK(approx_equal(s * s, a));
    }
  }

  SUBCASE("Geometric algebra identities (random)") {
    using index_set_t = fm_t::index_set_t; index_set_t frm = index_set_t();
    const double fill = 0.5;
    for (index_t i = 1; i <= 7; ++i) {
      frm |= index_set_t(i);
      frm |= index_set_t(-i);

      fm_t a = fm_t::random(frm, fill);
      fm_t b = fm_t::random(frm, fill);
      fm_t c = fm_t::random(frm, fill);

      // [HS] (1.25a): (a ^ b) ^ c == a ^ (b ^ c)
      CHECK(approx_equal((a ^ b) ^ c, a ^ (b ^ c)));

      // [HS] (1.31): a_1 * b == (a_1 & b) + (a_1 ^ b)
      fm_t a_1 = a(1);
      CHECK(approx_equal(a_1 * b, (a_1 & b) + (a_1 ^ b)));

      // [HS] (1.44): star(a, b) == scalar(a * b)
      CHECK(star(a, b) == doctest::Approx(numeric_traits<T>::to_double(scalar(a * b))));

      // [HS] (1.21a): (a_r * b_s)(|r-s|) == a_r & b_s (sampled)
      index_t r = frm.count() / 2;
      index_t s = frm.count() / 2;
      fm_t a_r = a(r);
      fm_t b_s = b(s);
      CHECK(approx_equal(a_r & b_s, (a_r * b_s)(index_t(std::abs(r-s)))));
    }
  }

  SUBCASE("Mixed-precision and assignment interchangeability") {
    using fm_f_t = glucat::framed_multi<float, -8, 8>;
    using fm_d_t = glucat::framed_multi<double, -8, 8>;

    fm_f_t f_f1("{1}");
    fm_f_t f_f2("{2}");
    fm_d_t f_d1("{1}");
    fm_d_t f_d2("{2}");

    // Compound assignment (matching scalar)
    f_f1 = fm_f_t("{1}");
    f_f1 += f_f2;
    CHECK(f_f1 == fm_f_t("{1}+{2}"));

    f_f1 = fm_f_t("{1}");
    f_f1 -= f_f2;
    CHECK(f_f1 == fm_f_t("{1}-{2}"));

    f_f1 = fm_f_t("{1}");
    f_f1 *= f_f2;
    CHECK(f_f1 == (fm_f_t("{1}") * f_f2));

    f_f1 = fm_f_t("{1}");
    f_f1 /= f_f2;
    CHECK(f_f1 == (fm_f_t("{1}") / f_f2));

    f_f1 = fm_f_t("{1}");
    f_f1 ^= f_f2;
    CHECK(f_f1 == (fm_f_t("{1}") ^ f_f2));

    f_f1 = fm_f_t("{1}");
    f_f1 &= f_f2;
    CHECK(f_f1 == (fm_f_t("{1}") & f_f2));

    f_f1 = fm_f_t("{1}");
    f_f1 %= f_f2;
    CHECK(f_f1 == (fm_f_t("{1}") % f_f2));

    f_f1 = fm_f_t("{1}");
    f_f1 |= f_f2;
    CHECK(f_f1 == (fm_f_t("{1}") | f_f2));

    // Interchangeability: A += B is interchangeable with A = A + B (matching scalar)
    fm_f_t A("{1}");
    fm_f_t B("{2}");
    fm_f_t C = A;
    C += B;
    A = A + B;
    CHECK(A == C);

    // Verify that mixed-precision assignment still works if we use the constructor explicitly
    f_f1 = fm_f_t(f_d1);
    CHECK(f_f1 == fm_f_t("{1}"));
  }

  SUBCASE("Mixed-representation assignment") {
    using mm_d_t = glucat::matrix_multi<double, -8, 8>;
    fm_t f("{1,2}");
    mm_d_t m("{1,2}");

    // Matrix to Framed (matching scalar)
    fm_t f2;
    f2 = m;
    CHECK(f2 == f);
  }

  SUBCASE("Advanced Constructor Gaps") {
    using index_set_t = fm_t::index_set_t;
    using scalar_t = fm_t::scalar_t;
    index_set_t ist(1);
    index_set_t frm(1);
    scalar_t crd(1.0);

    // Test framed_multi constructor with prechecked
    fm_t m1(ist, crd, frm, true);
    CHECK(m1.frame() == frm);

    fm_t m2(ist, crd, frm, false);
    CHECK(m2.frame() == frm);

    // Trigger exception path for constructor
    index_set_t invalid_ist(2);
    CHECK_THROWS_AS(fm_t(invalid_ist, crd, frm, false), glucat_error);
  }

  SUBCASE("Fast Path (32 terms)") {
    using index_set_t = fm_t::index_set_t;
    using scalar_t = fm_t::scalar_t;
    index_set_t large_frm;
    for(int i=1; i<=6; ++i) large_frm |= index_set_t(i); // 6 generators

    fm_t f;
    for(int i=0; i<32; ++i) {
       index_set_t ist;
       for(int j=0; j<6; ++j) if (i & (1 << j)) ist |= index_set_t(j+1);
       f += typename fm_t::term_t(ist, scalar_t(i+1));
    }

    CHECK(f.nbr_terms() >= 16);

    // Explicit conversion triggers constructor paths
    fm_t a(f);
    fm_t b(f, large_frm);

    fm_t res = a * b;
    res = a ^ b;
    res = a & b;
    res = a % b;
  }

  SUBCASE("Exceptions") {
    fm_t f1(1.0);
    // Parsing errors
    CHECK_THROWS(fm_t("{invalid}"));

    // Negative exponent
    CHECK_THROWS(f1.outer_pow(-1));

    // Construction with value outside of frame
    using is_t = fm_t::index_set_t;
    CHECK_THROWS(fm_t(is_t(1), T(1.0), is_t(), false));
  }

  SUBCASE("I/O and Formatting") {
    fm_t f(1.23456789);
    std::ostringstream oss;
    oss << std::setprecision(4) << std::fixed << f;
    CHECK(oss.str().find("1.2346") != std::string::npos);

    oss.str("");
    oss << std::scientific << f;
    CHECK(oss.str().find("1.234") != std::string::npos);
    CHECK(oss.str().find("e+00") != std::string::npos);
  }

  SUBCASE("Complexifier and Log") {
    fm_t f(T(1.0), is_t());
    // Log of a scalar doesn't need a complexifier if scalar > 0
    CHECK(log(f) == fm_t(T(0.0), is_t()));

    // Log of -1 needs a complexifier
    fm_t f_neg(T(-1.0));
    fm_t complexifier = glucat::complexifier(f_neg);
    // We expect this to work if i*j squares to -1
    CHECK_NOTHROW(log(f_neg, complexifier));
  }

  SUBCASE("Clifford Algebra Operations") {
    fm_t f1("{1}");
    fm_t f2("{2}");
    fm_t f12 = f1 * f2;
    fm_t f_mix("1+{1}+{1,2}");

    // Involute
    CHECK(f1.involute() == -f1);
    CHECK(f12.involute() == f12); // (-e1)*(-e2) = e12
    CHECK(f_mix.involute() == fm_t("1-{1}+{1,2}"));

    // Reverse
    CHECK(f1.reverse() == f1);
    CHECK(f12.reverse() == -f12); // e2*e1 = -e12
    CHECK(f_mix.reverse() == fm_t("1+{1}-{1,2}"));

    // Conjugate
    CHECK(f_mix.conj() == fm_t("1-{1}-{1,2}"));

    // Norm and Quad
    CHECK(f_mix.quad() == doctest::Approx(numeric_traits<T>::to_double(T(3.0))));
    CHECK(f_mix.norm() == doctest::Approx(numeric_traits<T>::to_double(T(3.0))));
    CHECK(f_mix.max_abs() == doctest::Approx(numeric_traits<T>::to_double(T(1.0))));
  }

  SUBCASE("Projections and Parts") {
    fm_t f_mix("1+{1}+{1,2}");

    CHECK(f_mix.pure() == fm_t("{1}+{1,2}"));
    CHECK(f_mix.even() == fm_t("1+{1,2}"));
    CHECK(f_mix.odd() == fm_t("{1}"));

    fm_t f_vec("1+2{1}+3{2}+4{1,2}");
    auto v = f_vec.vector_part();
    CHECK(v.size() == 2);
    CHECK(v[0] == doctest::Approx(numeric_traits<T>::to_double(T(2.0))));
    CHECK(v[1] == doctest::Approx(numeric_traits<T>::to_double(T(3.0))));

    using is_t = fm_t::index_set_t;
    auto v2 = f_vec.vector_part(is_t("{-1,1,2}"));
    CHECK(v2.size() == 3);
    CHECK(v2[0] == doctest::Approx(numeric_traits<T>::to_double(T(0.0))));
    CHECK(v2[1] == doctest::Approx(numeric_traits<T>::to_double(T(2.0))));
    CHECK(v2[2] == doctest::Approx(numeric_traits<T>::to_double(T(3.0))));
  }

  SUBCASE("Numerical Stability and Truncation") {
    using is_t = fm_t::index_set_t;
    fm_t f_small("1e8+{1}+1e-8{1,2}");
    CHECK(f_small.truncated(T(1.0e-6)) == fm_t(T(1.0e8), is_t()));

    fm_t f_nan(std::numeric_limits<T>::quiet_NaN(), is_t());
    fm_t f_inf(std::numeric_limits<T>::infinity(), is_t());

    CHECK(f_nan.isnan());
    CHECK(!f_nan.isinf());
    CHECK(f_inf.isinf());
    CHECK(!f_inf.isnan());
  }

  SUBCASE("Extended I/O") {
    fm_t f("1+{1}");
    f.write("Test prefix"); // Mostly for coverage

    auto temp_path = std::filesystem::temp_directory_path() / ("test_io_" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count()) + ".txt");
    struct Cleanup { std::filesystem::path p; ~Cleanup() { std::error_code ec; if(std::filesystem::exists(p, ec)) std::filesystem::remove(p, ec); } } cleanup{temp_path};

    std::ofstream ofs(temp_path);
    f.write(ofs, "File prefix");
    ofs.close();
    CHECK(std::filesystem::exists(temp_path));
  }

  SUBCASE("More Clifford Operations") {
    fm_t f1("{1}");
    fm_t f_mix("1+{1}+{1,2}");
    using is_t = fm_t::index_set_t;

    // Grade and Frame
    CHECK(f_mix.grade() == 2);
    CHECK(f_mix.frame() == is_t("{1,2}"));

    // Subscripting
    CHECK(f_mix[is_t("{1,2}")] == doctest::Approx(numeric_traits<T>::to_double(T(1.0))));

    // Pow and Inv
    fm_t f_inv = f1.inv();
    CHECK(f_inv == f1);
    CHECK(f1.pow(2) == fm_t(T(1.0), is_t()));

    // Outer power
    CHECK(f1.outer_pow(2) == fm_t(T(0.0), is_t()));
  }

  SUBCASE("Defensive and Edge Case Coverage") {
    using is_t = fm_t::index_set_t;
    using mm_t = matrix_multi<T, -8, 8>;

    // Outside-frame initialization (Line 181)
    CHECK_THROWS_AS(fm_t(fm_t("{1,2}"), is_t("{1}"), false), glucat_error);

    // Malformed string parsing (Line 320)
    CHECK_THROWS_AS(fm_t("1.0{1} garbage"), glucat_error);

    // NaN propagation from matrix (Lines 381-385)
    // Use a 2x2 matrix (dim=2) to avoid the "dim == 1" fast return at line 367
    // but stay below the default fast-path threshold of 4.
    mm_t m_nan_multi = mm_t(std::numeric_limits<T>::quiet_NaN(), is_t()) * mm_t(is_t("{1}"), 1.0);
    fm_t f_nan(m_nan_multi);
    CHECK(f_nan.isnan());

    // Zero matrix conversion (Line 361)
    mm_t m_zero = mm_t(T(0.0), is_t());
    fm_t f_zero(m_zero);
    CHECK(f_zero.nbr_terms() == 0);

    // 1x1 matrix conversion (Lines 366-369)
    mm_t m_1x1(T(3.14), is_t()); // scalar is effectively 1x1 matrix
    fm_t f_1x1(m_1x1);
    CHECK(f_1x1[is_t()] == doctest::Approx(numeric_traits<T>::to_double(T(3.14))));

    // Fast Path Conversion (Line 372)
    // Default tuning policy has threshold = 4.
    // 4 generators -> 8x8 matrix (dim=8), which triggers the fast path.
    mm_t m_8x8 = mm_t::random(is_t("{1,2,3,4}"), 1.0);
    fm_t f_fast(m_8x8);
    CHECK(f_fast.frame() == is_t("{1,2,3,4}"));

    // Equality and Comparison Edge Cases (Lines 419, 425, 449)
    fm_t f_a("{1}");
    fm_t f_b("{2}");
    fm_t f_empty;
    CHECK_FALSE(f_a == f_b);      // Different terms, same size
    CHECK_FALSE(f_a == f_empty);  // Different sizes
    CHECK(f_empty == T(0.0));     // Scalar zero equality (Line 449)
    CHECK_FALSE(f_a == T(0.0));   // Scalar zero inequality (Line 456)

    // Arithmetic Edge Cases (Lines 598, 601, 644)
    T nan = std::numeric_limits<T>::quiet_NaN();
    fm_t f_val("{1}");

    // operator*= with NaN (Line 598)
    fm_t f_nan_res = f_val;
    f_nan_res *= nan;
    CHECK(f_nan_res.isnan());

    // operator*= 0.0 with NaN (Line 601)
    fm_t f_nan_val(nan, is_t());
    f_nan_val *= T(0.0);
    CHECK(f_nan_val.isnan());

    // operator* with NaN (Line 644)
    CHECK((f_nan_val * f_val).isnan());
    CHECK((f_val * f_nan_val).isnan());
  }
}
#endif

#endif  // _GLUCAT_FRAMED_MULTI_IMP_H
