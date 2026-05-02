#ifndef _GLUCAT_MATRIX_MULTI_IMP_H
#define _GLUCAT_MATRIX_MULTI_IMP_H
/**************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_multi_imp.h : Implement the matrix representation of a multivector
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

#include "glucat/matrix_multi.h"

#include "glucat/scalar.h"
#include "glucat/generation.h"
#include "glucat/matrix.h"

#include <iomanip>
#include <array>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <limits>

namespace glucat
{
  // References for algorithms:
  // [CHKL]:
  // [L]: Pertti Lounesto, "Clifford algebras and spinors", Cambridge UP, 1997.
  // [MB]: Beatrice Meini, "The Matrix Square Root From a New Functional Perspective:
  // Theoretical Results and Computational Issues", SIAM Journal on
  // Matrix Analysis and Applications 26(2):362-376, 2004.
  // [P]: Ian R. Porteous, "Clifford algebras and the classical groups", Cambridge UP, 1995.

  /*
   * @brief Class name used in messages
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline std::string_view
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  classname()
  { return "matrix_multi"; }

  /*
   * @brief Determine the log2 dim corresponding to signature p, q
   * @details
   * @param p First index
   * @param q Second index
   * @return Result
   */
  // Reference: [P] Table 15.27, p 133
  inline index_t
  offset_level(const index_t p, const index_t q)
  {
    // Offsets between the log2 of the matrix dimension for the current signature
    // and that of the real superalgebra
    static const std::array<int, 8> offset_log2_dim = {0, 1, 0, 1, 1, 2, 1, 1};
    const index_t bott = pos_mod(p-q, 8);
    return (p+q)/2 + offset_log2_dim[bott];
  }

  /*
   * @brief Determine the matrix dimension of the fold of a subalegbra
   * @details
   * @param sub Value
   */
  // Reference: [P] Table 15.27, p 133
  template< typename Matrix_Index_T, int LO, int HI >
  inline Matrix_Index_T
  folded_dim( const index_set<LO,HI>& sub )
  { return 1 << offset_level(sub.count_pos(), sub.count_neg()); }

  /*
   * @brief Default constructor
   * @details
   *
   * Usage example:
   * Location: glucat/framed_multi_imp.h:1136
   *
   * @code
   *
   * auto result = matrix_multi_t(Scalar_T(1), this->frame());
   * @endcode
   *
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi()
  : m_frame( index_set_t() ),
    m_matrix( matrix_t( 1, 1 ) )
  { this->m_matrix.zeros(); }

  /*
   * @brief Move constructor
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param other Other matrix
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(matrix_multi&& other) noexcept(std::is_nothrow_move_constructible_v<Scalar_T>)
  : m_frame(std::move(other.m_frame)),
    m_matrix(std::move(other.m_matrix))
  { }

  /*
   * @brief Construct a multivector from a multivector with a different scalar type
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  template< typename Other_Scalar_T, typename Other_Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const matrix_multi<Other_Scalar_T,LO,HI,Other_Tune_P>& val)
  : m_frame( val.m_frame ), m_matrix( val.m_matrix.nbr_rows(), val.m_matrix.nbr_cols() )
  {
    this->m_matrix.zeros();
    for (matrix_index_t i = 0; i < val.m_matrix.nbr_rows(); ++i)
      for (matrix_index_t j = 0; j < val.m_matrix.nbr_cols(); ++j)
        this->m_matrix(i, j) = numeric_traits<Scalar_T>::to_scalar_t(val.m_matrix(i, j));
  }

  /*
   * @brief Construct a multivector, within a given frame, from a given multivector
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   * @param frm Value
   * @param prechecked Already checked?
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  template< typename Other_Scalar_T, typename Other_Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const matrix_multi<Other_Scalar_T,LO,HI,Other_Tune_P>& val, const index_set_t frm, const bool prechecked)
  : m_frame( frm )
  {
    if (frm != val.m_frame)
      *this = multivector_t(framed_multi_t(val), frm);
    else
    {
      const matrix_index_t dim = folded_dim<matrix_index_t>(frm);
      this->m_matrix.zeros(dim, dim);
      for (matrix_index_t i = 0; i < val.m_matrix.nbr_rows(); ++i)
        for (matrix_index_t j = 0; j < val.m_matrix.nbr_cols(); ++j)
          this->m_matrix(i, j) = numeric_traits<Scalar_T>::to_scalar_t(val.m_matrix(i, j));
    }
  }

  /*
   * @brief Construct a multivector, within a given frame, from a given multivector
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   * @param frm Value
   * @param prechecked Already checked?
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const multivector_t& val, const index_set_t frm, const bool prechecked)
  : m_frame( frm )
  {
    if (frm != val.m_frame)
      *this = multivector_t(framed_multi_t(val), frm);
    else
      this->m_matrix = val.m_matrix;
  }

  /*
   * @brief Construct a multivector from an index set and a scalar coordinate
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param ist Value
   * @param crd Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const index_set_t ist, const Scalar_T& crd)
  : m_frame( ist )
  {
    const auto dim = folded_dim<matrix_index_t>(this->m_frame);
    this->m_matrix.zeros(dim, dim);
    *this += term_t(ist, crd);
  }

  /*
   * @brief Construct a multivector, within a given frame, from an index set and a scalar coordinate
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param ist Value
   * @param crd Value
   * @param frm Value
   * @param prechecked Already checked?
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const index_set_t ist, const Scalar_T& crd, const index_set_t frm, const bool prechecked)
  : m_frame( frm )
  {
    if (!prechecked && (ist | frm) != frm)
      throw error_t("multivector_t(ist,crd,frm): cannot initialize with value outside of frame");
    const matrix_index_t dim = folded_dim<matrix_index_t>(frm);
    this->m_matrix.zeros(dim, dim);
    *this += term_t(ist, crd);
  }

  /*
   * @brief Construct a multivector from a scalar (within a frame, if given)
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param scr Value
   * @param frm Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const Scalar_T& scr, const index_set_t frm)
  : m_frame( frm )
  {
    const auto dim = folded_dim<matrix_index_t>(frm);
    this->m_matrix.zeros(dim, dim);
    *this += term_t(index_set_t(), scr);
  }

  /*
   * @brief Construct a multivector from an int (within a frame, if given)
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param scr Value
   * @param frm Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const int scr, const index_set_t frm)
  { *this = multivector_t(Scalar_T(scr), frm); }

  /*
   * @brief Construct a multivector, within a given frame, from a given vector
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param vec Value
   * @param frm Value
   * @param prechecked Already checked?
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const vector_t& vec,
               const index_set_t frm, const bool prechecked)
  : m_frame( frm )
  {
    if (!prechecked && index_t(vec.size()) != frm.count())
      throw error_t("multivector_t(vec,frm): cannot initialize with vector not matching frame");
    const auto dim = folded_dim<matrix_index_t>(frm);
    this->m_matrix.zeros(dim, dim);
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
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param str Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const std::string& str)
  { *this = multivector_t(framed_multi_t(str)); }

  /*
   * @brief Construct a multivector, within a given frame, from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param str Value
   * @param frm Value
   * @param prechecked Already checked?
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const std::string& str, const index_set_t frm, const bool prechecked)
  { *this = multivector_t(framed_multi_t(str), frm, prechecked); }

  /*
   * @brief Construct a multivector from a framed_multi_t
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  template< typename Other_Scalar_T, typename Other_Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const framed_multi<Other_Scalar_T,LO,HI,Other_Tune_P>& val)
  : m_frame( val.frame() )
  {
    using Tuning_Values_P = typename Tune_P::tuning_values_p;
    if (val.size() >= Tuning_Values_P::fast_size_threshold)
      try
      {
        auto tmp = val.template fast_matrix_multi<Scalar_T,Tune_P>(this->m_frame);

        *this = tmp;
        return;
      }
      catch (const glucat_error& e)
      { }
    const auto dim = folded_dim<matrix_index_t>(this->m_frame);
    if (dim == 0) {

    }
    this->m_matrix.zeros(dim, dim);

    for (auto& val_term : val)
      *this += val_term;
  }

  /*
   * @brief Construct a multivector, within a given frame, from a framed_multi_t
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param framed_val Value
   * @param frm Value
   * @param prechecked Already checked?
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  template< typename Other_Scalar_T, typename Other_Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const framed_multi<Other_Scalar_T,LO,HI,Other_Tune_P>& framed_val, const index_set_t frm, const bool prechecked)
  {
    using Tuning_Values_P = typename Tune_P::tuning_values_p;
    const auto val = framed_val.truncated();
    const auto our_frame = val.frame() | frm;
    if (val.size() >= Tuning_Values_P::fast_size_threshold)
      try
      {
        *this = val.template fast_matrix_multi<Scalar_T,Tune_P>(our_frame);
        return;
      }
      catch (const glucat_error& e)
      { }
    this->m_frame = our_frame;
    const auto dim = folded_dim<matrix_index_t>(our_frame);
    this->m_matrix.zeros(dim, dim);

    for (auto& val_term : val)
      *this += val_term;
  }

  /*
   * @brief Construct a multivector within a given frame from a given matrix
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param mtx Value
   * @param frm Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  template< typename Matrix_T >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const Matrix_T& mtx, const index_set_t frm)
  : m_frame( frm )
  {
    if constexpr (requires { this->m_matrix = mtx; })
      this->m_matrix = mtx;
    else
    {
      const matrix_index_t mtx_nbr_rows = matrix::nbr_rows(mtx);
      const matrix_index_t mtx_nbr_cols = matrix::nbr_cols(mtx);

      this->m_matrix.set_size(mtx_nbr_rows, mtx_nbr_cols);
      for (matrix_index_t i = 0; i < mtx_nbr_rows; ++i)
        for (matrix_index_t j = 0; j < mtx_nbr_cols; ++j)
          this->m_matrix(i, j) = numeric_traits<Scalar_T>::to_scalar_t(mtx(i, j));
    }
  }

  /*
   * @brief Construct a multivector within a given frame from a given matrix
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param mtx Value
   * @param frm Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  matrix_multi(const matrix_t& mtx, const index_set_t frm)
  : m_frame( frm ), m_matrix( mtx )
  { }

  /*
   * @brief Find a common frame for operands of a binary operator
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @param lhs_reframed Value
   * @param rhs_reframed Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline index_set<LO,HI>
  reframe (const matrix_multi<Scalar_T,LO,HI,Tune_P>& lhs,    const matrix_multi<Scalar_T,LO,HI,Tune_P>& rhs,
                 matrix_multi<Scalar_T,LO,HI,Tune_P>& lhs_reframed, matrix_multi<Scalar_T,LO,HI,Tune_P>& rhs_reframed)
  {
    using index_set_t = index_set<LO, HI>;
    using multivector_t = matrix_multi<Scalar_T,LO,HI,Tune_P>;
    using framed_multi_t = typename multivector_t::framed_multi_t;
    // Determine the initial common frame
    index_set_t our_frame = lhs.m_frame | rhs.m_frame;
    framed_multi_t framed_lhs;
    framed_multi_t framed_rhs;
    if ((lhs.m_frame != our_frame) || (rhs.m_frame != our_frame))
    {
      // The common frame may expand as a result of the transform to framed_multi_t
      framed_lhs = framed_multi_t(lhs);
      framed_rhs = framed_multi_t(rhs);
      our_frame |= framed_lhs.frame() | framed_rhs.frame();
    }
    // Do the reframing only where necessary
    if (lhs.m_frame != our_frame)
      lhs_reframed = multivector_t(framed_lhs, our_frame, true);
    if (rhs.m_frame != our_frame)
      rhs_reframed = multivector_t(framed_rhs, our_frame, true);
    return our_frame;
  }

  /*
   * @brief Test for equality of multivectors
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param rhs Right hand side
   * @return True if equal
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline bool
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator== (const multivector_t& rhs) const
  {
    // Ensure that there is no aliasing
    if (this == &rhs)
      return true;

    // Operate only within a common frame
    multivector_t lhs_reframed;
    multivector_t rhs_reframed;
    const index_set_t our_frame = reframe(*this, rhs, lhs_reframed, rhs_reframed);
    const multivector_t& lhs_ref = (this->m_frame == our_frame)
      ? *this
      : lhs_reframed;
    const multivector_t& rhs_ref = (rhs.m_frame == our_frame)
      ? rhs
      : rhs_reframed;

    return (lhs_ref.m_matrix - rhs_ref.m_matrix).norm_inf() == 0;
  }

  /*
   * @brief Test for equality of multivector and scalar
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param scr Value
   * @return True if equal
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline bool
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator== (const Scalar_T& scr) const
  {
    if (scr != Scalar_T(0))
      return *this == multivector_t(scr, this->m_frame);
    else if (this->m_matrix.norm_inf() != 0)
      return false;
    else
    {
      const matrix_index_t dim = this->m_matrix.nbr_rows();
      return !(dim == 1 && this->isnan());
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
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator+= (const Scalar_T& scr)
  { return *this += term_t(index_set_t(), scr); }

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
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator+= (const multivector_t& rhs)
  {
    // Ensure that there is no aliasing
    if (this == &rhs)
      return *this *= Scalar_T(2);

    // Operate only within a common frame
    multivector_t rhs_reframed;
    const index_set_t our_frame = reframe(*this, rhs, *this, rhs_reframed);
    const multivector_t& rhs_ref = (rhs.m_frame == our_frame)
      ? rhs
      : rhs_reframed;

    this->m_matrix += rhs_ref.m_matrix;
    return *this;
  }

  _GLUCAT_CLIFFORD_ALGEBRA_ASSIGNMENT_OPERATIONS_IMP(matrix_multi)

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
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator-= (const Scalar_T& scr)
  { return *this += term_t(index_set_t(), -scr); }

  /*
   * @brief Geometric difference
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param rhs Right hand side
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator-= (const multivector_t& rhs)
  {
    // Ensure that there is no aliasing
    if (this == &rhs)
      return *this = Scalar_T(0);

    // Operate only within a common frame
    multivector_t rhs_reframed;
    const index_set_t our_frame = reframe(*this, rhs, *this, rhs_reframed);
    const multivector_t& rhs_ref = (rhs.m_frame == our_frame)
      ? rhs
      : rhs_reframed;

    this->m_matrix -= rhs_ref.m_matrix;
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
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator- () const
  { return multivector_t(-(this->m_matrix), this->m_frame); }

  /*
   * @brief Product of multivector and scalar
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param scr Value
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator*= (const Scalar_T& scr)
  { // multiply coordinates of all terms by scalar

    using traits_t = numeric_traits<Scalar_T>;
    if (traits_t::isNaN_or_isInf(scr) || this->isnan())
      return *this = traits_t::NaN();
    if (scr == Scalar_T(0))
      *this = Scalar_T(0);
    else
      this->m_matrix *= scr;
    return *this;
  }

  /*
   * @brief Geometric product
   * @details
   *
   * Usage example:
   * Location: glucat/framed_multi_imp.h:1101
   *
   * @code
   *
   * return matrix_multi_t(rhs) * matrix_multi_t(lhs) / matrix_multi_t(rhs.involute());
   * @endcode
   *
   * @par Example:
   * @code
   * clifford<>("{1}") * clifford<>("{2}"); // Returns {1,2}
   * @endcode
   *
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Product
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline matrix_multi<Scalar_T,LO,HI,Tune_P>
  operator* (const matrix_multi<Scalar_T,LO,HI,Tune_P>& lhs, const matrix_multi<Scalar_T,LO,HI,Tune_P>& rhs)
  {
    using multivector_t = matrix_multi<Scalar_T,LO,HI,Tune_P>;
    using index_set_t = typename multivector_t::index_set_t;

    if (lhs.isnan() || rhs.isnan())
      return multivector_t(numeric_traits<Scalar_T>::NaN());

    // Operate only within a common frame
    multivector_t lhs_reframed;
    multivector_t rhs_reframed;
    const index_set_t our_frame = reframe(lhs, rhs, lhs_reframed, rhs_reframed);
    const multivector_t& lhs_ref = (lhs.m_frame == our_frame)
      ? lhs
      : lhs_reframed;
    const multivector_t& rhs_ref = (rhs.m_frame == our_frame)
      ? rhs
      : rhs_reframed;

    using matrix_t = typename multivector_t::matrix_t;
    using matrix_index_t = typename matrix_t::size_type;

    const matrix_index_t dim = lhs_ref.m_matrix.nbr_rows();
    multivector_t result = multivector_t(matrix_t(dim, dim), our_frame);
    result.m_matrix.zeros();
    result.m_matrix = lhs_ref.m_matrix * rhs_ref.m_matrix;
    if (result.m_matrix.nbr_rows() == 0) {

    }
    return result;
  }

  /*
   * @brief Geometric product
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param rhs Right hand side
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator*= (const multivector_t& rhs)
  { return *this = *this * rhs; }

  /*
   * @brief Outer product
   * @details
   * @par Example:
   * @code
   * clifford<>("{1}") ^ clifford<>("{2}"); // Returns {1,2}
   * @endcode
   *
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Outer product
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline matrix_multi<Scalar_T,LO,HI,Tune_P>
  operator^ (const matrix_multi<Scalar_T,LO,HI,Tune_P>& lhs, const matrix_multi<Scalar_T,LO,HI,Tune_P>& rhs)
  {
    using multivector_t = matrix_multi<Scalar_T,LO,HI,Tune_P>;
    using framed_multi_t = typename multivector_t::framed_multi_t;
    return multivector_t(framed_multi_t(lhs) ^ framed_multi_t(rhs));
  }

  /*
   * @brief Outer product
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param rhs Right hand side
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator^= (const multivector_t& rhs)
  { return *this = *this ^ rhs; }

  /*
   * @brief Inner product
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
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
  inline matrix_multi<Scalar_T,LO,HI,Tune_P>
  operator& (const matrix_multi<Scalar_T,LO,HI,Tune_P>& lhs, const matrix_multi<Scalar_T,LO,HI,Tune_P>& rhs)
  {
    using multivector_t = matrix_multi<Scalar_T,LO,HI,Tune_P>;
    using framed_multi_t = typename multivector_t::framed_multi_t;
    return multivector_t(framed_multi_t(lhs) & framed_multi_t(rhs));
  }

  /*
   * @brief Inner product
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param rhs Right hand side
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator&= (const multivector_t& rhs)
  { return *this = *this & rhs; }

  /*
   * @brief Left contraction
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
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
  inline matrix_multi<Scalar_T,LO,HI,Tune_P>
  operator% (const matrix_multi<Scalar_T,LO,HI,Tune_P>& lhs, const matrix_multi<Scalar_T,LO,HI,Tune_P>& rhs)
  {
    using multivector_t = matrix_multi<Scalar_T,LO,HI,Tune_P>;
    using framed_multi_t = typename multivector_t::framed_multi_t;
    return multivector_t(framed_multi_t(lhs) % framed_multi_t(rhs));
  }

  /*
   * @brief Left contraction
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param rhs Right hand side
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator%= (const multivector_t& rhs)
  { return *this = *this % rhs; }

  /*
   * @brief Hestenes scalar product
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
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline Scalar_T
  star(const matrix_multi<Scalar_T,LO,HI,Tune_P>& lhs, const matrix_multi<Scalar_T,LO,HI,Tune_P>& rhs)
  { return (lhs * rhs).scalar(); }

  /*
   * @brief Quotient of multivector and scalar
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param scr Value
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator/= (const Scalar_T& scr)
  { return *this *= Scalar_T(1)/scr; }

  /*
   * @brief Geometric quotient
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
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
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
  matrix_multi<Scalar_T,LO,HI,Tune_P>
  operator/ (const matrix_multi<Scalar_T,LO,HI,Tune_P>& lhs, const matrix_multi<Scalar_T,LO,HI,Tune_P>& rhs)
  {
    using traits_t = numeric_traits<Scalar_T>;

    using multivector_t = matrix_multi<Scalar_T,LO,HI,Tune_P>;
    if (lhs.isnan() || rhs.isnan())
      return multivector_t(traits_t::NaN());

    if (rhs == Scalar_T(0))
      return multivector_t(traits_t::NaN());

    // Operate only within a common frame
    multivector_t lhs_reframed;
    multivector_t rhs_reframed;
    const auto our_frame = reframe(lhs, rhs, lhs_reframed, rhs_reframed);
    const auto& lhs_ref = (lhs.m_frame == our_frame)
      ? lhs
      : lhs_reframed;
    const auto& rhs_ref = (rhs.m_frame == our_frame)
      ? rhs
      : rhs_reframed;

    // Solve result == lhs_ref/rhs_ref <=> result*rhs_ref == lhs_ref
    // We now solve X == B/A
    // (where X == result, B == lhs_ref.m_matrix and A == rhs_ref.m_matrix)
    // X == B/A <=> X*A == B <=> AT*XT == BT
    // So, we solve AT*XT == BT

    using matrix_t = typename multivector_t::matrix_t;

    const auto& AT = matrix_t(rhs_ref.m_matrix.t());
    const auto& BT = matrix_t(lhs_ref.m_matrix.t());
    matrix_t XT(AT.nbr_rows(), AT.nbr_cols());

    // Solve AT * XT = BT
    if (matrix::solve(XT, AT, BT))
      return multivector_t(XT.t(), our_frame);
    else
      return multivector_t(traits_t::NaN());
  }

  /*
   * @brief Geometric quotient
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param rhs Right hand side
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator/= (const multivector_t& rhs)
  { return *this = *this / rhs; }

  /*
   * @brief Transformation via twisted adjoint action
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @return Bitwise OR
   *
   * @par Example:
   * @code
   * clifford<>("{2}") | clifford<>("{1,2}"); // Returns -{2}
   * @endcode
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline matrix_multi<Scalar_T,LO,HI,Tune_P>
  operator| (const matrix_multi<Scalar_T,LO,HI,Tune_P>& lhs, const matrix_multi<Scalar_T,LO,HI,Tune_P>& rhs)
  { return rhs * lhs / rhs.involute(); }

  /*
   * @brief Transformation via twisted adjoint action
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param rhs Right hand side
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator|= (const multivector_t& rhs)
  { return *this = rhs * *this / rhs.involute(); }

  /*
   * @brief Clifford multiplicative inverse
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  inv() const
  { return multivector_t(Scalar_T(1), this->m_frame) / *this; }

  /*
   * @brief Integer power of multivector: *this to the m
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param m Matrix
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  pow(int m) const
  { return glucat::pow(*this, m); }

  /*
   * @brief Move assignment
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param other Other matrix
   * @return Reference to this
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator= (matrix_multi&& other) noexcept(std::is_nothrow_move_assignable_v<Scalar_T>)
  {
    this->m_frame = std::move(other.m_frame);
    this->m_matrix = std::move(other.m_matrix);
    return *this;
  }

  /*
   * @brief Outer product power of multivector
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param m Matrix
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  outer_pow(int m) const
  {
    if (m < 0)
      throw error_t("outer_pow(m): negative exponent");
    framed_multi_t a(*this);
    return multivector_t(a.outer_pow(m));
  }

  /*
   * @brief Grade of multivector: maximum of the grades of each term
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Grade
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline index_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  grade() const
  { return framed_multi_t(*this).grade(); }

  /*
   * @brief Frame of multivector: union of index sets of terms
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::index_set_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  frame() const
  { return this->m_frame; }

  /*
   * @brief Subscripting: map from index set to scalar coordinate
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param ist Value
   * @return Element reference
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline Scalar_T
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator[] (const index_set_t ist) const
  {
    // Use matrix inner product only if ist is in frame
    if ( (ist | this->m_frame) == this->m_frame)
      return (this->basis_element(ist)).template inner<Scalar_T>(this->m_matrix);
    else
      return Scalar_T(0);
  }

  /*
   * @brief Grading: part where each term is a grade-vector
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Element
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator() (index_t grade) const
  {
    if ((grade < 0) || (grade > HI-LO))
      return 0;
    else
      return multivector_t((framed_multi_t(*this))(grade));
  }

  /*
   * @brief Scalar part
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Scalar part
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline Scalar_T
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  scalar() const
  {
    const matrix_index_t dim = this->m_matrix.nbr_rows();
    return this->m_matrix.trace() / Scalar_T( double(dim) );
  }

  /*
   * @brief Pure part
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Pure part
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  pure() const
  { return *this - this->scalar(); }

  /*
   * @brief Even part, sum of the even grade terms
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Even part
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  even() const
  { return multivector_t(framed_multi_t(*this).even()); }

  /*
   * @brief Odd part, sum of the odd grade terms
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Odd part
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  odd() const
  { return multivector_t(framed_multi_t(*this).odd()); }

  /*
   * @brief Vector part of multivector, as a vector_t
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Vector part
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::vector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  vector_part() const
  { return this->vector_part(this->frame(), true); }

  /*
   * @brief Vector part of multivector, as a vector_t
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param frm Value
   * @param prechecked Already checked?
   * @return True if successful or condition met
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::vector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  vector_part(const index_set_t frm, const bool prechecked) const
  {
    if (!prechecked && (this->frame() | frm) != frm)
      throw error_t("vector_part(frm): value is outside of requested frame");
    vector_t result;
    // If we need to enlarge the frame we may as well use a framed_multi_t
    if (this->frame() != frm)
      return framed_multi_t(*this).vector_part(frm, true);

    const auto begin_index = frm.min();
    const auto end_index = frm.max()+1;
    for (auto
        idx = begin_index;
        idx != end_index;
        ++idx)
      if (frm[idx])
        // Frame may contain indices which do not correspond to a grade 1 term but
        // frame cannot omit any index corresponding to a grade 1 term
        result.push_back(
          (this->basis_element(index_set_t(idx))).template inner<Scalar_T>(this->m_matrix));
    return result;
  }

  /*
   * @brief Main involution, each {i} is replaced by -{i} in each term
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Main involution
   *
   * @par Example:
   * @code
   * involute(clifford<>("{1}")); // Returns -{1}
   * @endcode
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  involute() const
  { return multivector_t(framed_multi_t(*this).involute()); }

  /*
   * @brief Reversion, order of {i} is reversed in each term
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Reverse
   *
   * @par Example:
   * @code
   * reverse(clifford<>("{1}*{2}")); // Returns -{1}*{2}
   * @endcode
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  reverse() const
  { return multivector_t(framed_multi_t(*this).reverse()); }

  /*
   * @brief Conjugation, conj == reverse o involute == involute o reverse
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Clifford conjugate
   *
   * @par Example:
   * @code
   * conj(clifford<>("{1}")); // Returns -{1}
   * @endcode
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  conj() const
  { return multivector_t(framed_multi_t(*this).conj()); }

  /*
   * @brief Quadratic form := scalar part of rev(x)*x
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Quadratic form
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline Scalar_T
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  quad() const
  { // scalar(conj(x)*x) = 2*quad(even(x)) - quad(x)
    // Arvind Raja ref: "old clical: quadfunction(p:pter):pterm in file compmod.pas"
    return framed_multi_t(*this).quad();
  }

  /*
   * @brief Scalar_T norm squared= sum of norm squared of coordinates
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Norm
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline Scalar_T
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  norm() const
  {
    const matrix_index_t dim = this->m_matrix.nbr_rows();
    return this->m_matrix.norm_frob2() / Scalar_T( double(dim) );
  }

  /*
   * @brief Maximum of absolute values of components of multivector: multivector infinity norm
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Maximum absolute value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline Scalar_T
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  max_abs() const
  { return framed_multi_t(*this).max_abs(); }

  /*
   * @brief Random multivector within a frame
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param frm Value
   * @param fill Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  random(const index_set<LO,HI> frm, Scalar_T fill)
  {
    return multivector_t(framed_multi<Scalar_T,LO,HI,Tune_P>::random(frm, fill));
  }

  /*
   * @brief Write multivector to output
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param msg Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  void
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  write(const std::string& msg) const
  { framed_multi_t(*this).write(msg); }

  /*
   * @brief Write out multivector to file
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param ofile Value
   * @param msg Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  void
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  write(std::ofstream& ofile, const std::string& msg) const
  {
    if (!ofile)
      throw error_t("write(ofile,msg): cannot write to output file");
    framed_multi_t(*this).write(ofile, msg);
  }

  /*
   * @brief Write multivector to output
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param os Output stream
   * @param val Value
   * @return Output stream
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline std::ostream&
  operator<< (std::ostream& os, const matrix_multi<Scalar_T,LO,HI,Tune_P>& val)
  {
    os << typename matrix_multi<Scalar_T,LO,HI,Tune_P>::framed_multi_t(val);
    return os;
  }

  /*
   * @brief Read multivector from input
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param s Value
   * @param val Value
   * @return Input stream
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline std::istream&
  operator>> (std::istream& s, matrix_multi<Scalar_T,LO,HI,Tune_P>& val)
  { // Input looks like 1.0-2.0{1,2}+3.2{3,4}
    framed_multi<Scalar_T,LO,HI,Tune_P> local;
    s >> local;
    // If s.bad() then we have a corrupt input
    // otherwise we are fine and can copy the resulting matrix_multi
    if (!s.bad())
      val = matrix_multi<Scalar_T,LO,HI,Tune_P>(local);
    return s;
  }

  /*
   * @brief Check if a multivector contains any infinite values
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return True if successful or condition met
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline bool
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  isinf() const
  {
    if (std::numeric_limits<Scalar_T>::has_infinity)
      return this->m_matrix.isinf();
    else
      return false;
  }

  /*
   * @brief Check if a multivector contains any IEEE NaN values
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return True if successful or condition met
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline bool
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  isnan() const
  {
    if (std::numeric_limits<Scalar_T>::has_quiet_NaN)
      return this->m_matrix.isnan();
    else
      return false;
  }

  /*
   * @brief Number of rows
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Number of rows
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::matrix_index_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  nbr_rows() const
  { return matrix::nbr_rows(this->m_matrix); }

  /*
   * @brief Number of columns
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Number of columns
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::matrix_index_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  nbr_cols() const
  { return matrix::nbr_cols(this->m_matrix); }

  /*
   * @brief Remove all terms with relative size smaller than limit
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param limit Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  truncated(const Scalar_T& limit) const
  { return multivector_t(framed_multi_t(*this).truncated(limit)); }

  /*
   * @brief Add a term.
   * @param term The term to add.
   * @return Reference to this.
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t&
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  operator+= (const term_t& term)
  {
    if (term.second != Scalar_T(0))
      this->m_matrix += matrix_t(this->basis_element(term.first)) * term.second;
    return *this;
  }

  /*
   * @brief Inverse generalized Fast Fourier Transform
   * @details
   * @tparam Multivector_T
   * @tparam Matrix_T
   * @param X Value
   * @param level Recursion level
   * @return Result
   */
  template<typename Multivector_T, typename Matrix_T>
  inline Multivector_T
  fast(const Matrix_T& X, index_t level)
  {
    using framed_multi_t = Multivector_T;

    using index_set_t = typename framed_multi_t::index_set_t;
    using Scalar_T = typename framed_multi_t::scalar_t;
    using matrix_t = Matrix_T;
    using matrix_scalar_t = typename matrix_t::value_type;
    using traits_t = numeric_traits<Scalar_T>;

    if (level == 0)
      return framed_multi_t(traits_t::to_scalar_t(X(0,0)));

    if (X.norm_inf() == 0)
      return Scalar_T(0);


    const matrix_t&  I = matrix::unit<matrix_t>(2);
    matrix_t J(2,2);
    J.zeros();
    J(0,1)  = matrix_scalar_t(-1);
    J(1,0)  = matrix_scalar_t( 1);
    matrix_t K = J;
    K(0,1)  = matrix_scalar_t( 1);
    matrix_t JK = I;
    JK(0,0) = matrix_scalar_t(-1);

    const index_set_t ist_mn   = index_set_t(-level);
    const index_set_t ist_pn   = index_set_t(level);
    const index_set_t ist_mnpn = ist_mn | ist_pn;
    if (level == 1)
    {
      using term_t = typename framed_multi_t::term_t;
      const Scalar_T i_x  = traits_t::to_scalar_t(I.nork(X, false)(0, 0));
      const Scalar_T j_x  = traits_t::to_scalar_t(J.nork(X, false)(0, 0));
      const Scalar_T k_x  = traits_t::to_scalar_t(K.nork(X, false)(0, 0));
      const Scalar_T jk_x = traits_t::to_scalar_t(JK.nork(X, false)(0, 0));

      framed_multi_t result = i_x;
      result += term_t(ist_mn,   j_x);  // j_x *  mn;
      result += term_t(ist_pn,   k_x);  // k_x *  pn;
      return result += term_t(ist_mnpn, jk_x); // jk_x * mnpn;
    }
    else
    {
      const framed_multi_t& mn   = framed_multi_t(ist_mn);
      const framed_multi_t& pn   = framed_multi_t(ist_pn);
      const framed_multi_t& mnpn = framed_multi_t(ist_mnpn);
      const framed_multi_t& i_x  = fast<framed_multi_t, matrix_t>
                                       (I.nork(X, false), level-1);
      const framed_multi_t& j_x  = fast<framed_multi_t, matrix_t>
                                       (J.nork(X, false), level-1);
      const framed_multi_t& k_x  = fast<framed_multi_t, matrix_t>
                                       (K.nork(X, false), level-1);
      const framed_multi_t& jk_x = fast<framed_multi_t, matrix_t>
                                       (JK.nork(X, false), level-1);
      framed_multi_t
             result  =  i_x.even() - jk_x.odd();
             result += (j_x.even() - k_x.odd()) * mn;
             result += (k_x.even() - j_x.odd()) * pn;
      return result += (jk_x.even() - i_x.odd()) * mnpn;
    }
  }

  /*
   * @brief Use generalized FFT to construct a matrix_multi_t
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param frm Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::multivector_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  fast_matrix_multi(const index_set_t frm) const
  {
    if (this->m_frame == frm)
      return *this;
    else
      return (this->template fast_framed_multi<Scalar_T,Tune_P>()).template fast_matrix_multi<Scalar_T,Tune_P>(frm);
  }

  /*
   * @brief Use inverse generalized FFT to construct a framed_multi_t
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  template< typename Other_Scalar_T, typename Other_Tune_P >
  auto
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  fast_framed_multi() const -> framed_multi<Other_Scalar_T,LO,HI,Other_Tune_P>
  {
    // Determine the amount of off-centering needed
    index_t p = this->m_frame.count_pos();
    index_t q = this->m_frame.count_neg();

    const index_t bott = pos_mod(p-q, 8);
    p += std::max(gen::offset_to_super[bott],index_t(0));
    q -= std::min(gen::offset_to_super[bott],index_t(0));

    const index_t orig_p = p;
    const index_t orig_q = q;
    while (p-q > 4)
      { p -= 4; q += 4; }
    while (p-q < -3)
      { p += 4; q -= 4; }
    if (p-q > 1)
    {
      index_t old_p = p;
      p = q+1;
      q = old_p-1;
    }
    const index_t level = (p+q)/2;
#if defined(_GLUCAT_MATRIX_MULTI_IMP_DEBUG)
    std::cout << "DEBUG: In matrix_multi_t::fast_framed_multi, dim(m_matrix) is " << (this->m_matrix).nbr_cols() << std::endl;
    std::cout << "DEBUG: m_matrix is" << std::endl;
    std::cout << (this->m_matrix) << std::endl;
#endif
    // Do the inverse fast transform
    using framed_multi_t = framed_multi<Other_Scalar_T,LO,HI,Other_Tune_P>;
    framed_multi_t val = fast<framed_multi_t, matrix_t>(this->m_matrix, level);

    // Off-centre val
    switch (pos_mod(orig_p-orig_q, 8))
    {
    case 2:
    case 3:
    case 4:
      val.centre_qp1_pm1(p, q);
      break;
    default:
      break;
    }
    if (orig_p-orig_q > 4)
      while (p != orig_p)
        val.centre_pp4_qm4(p, q);
    if (orig_p-orig_q < -3)
      while (p != orig_p)
        val.centre_pm4_qp4(p, q);

    // Return unfolded val
    return val.unfold(this->m_frame);
  }

  /*
   * @brief Table of basis elements used as a cache by basis_element()
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Matrix_T
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Matrix_T >
  class basis_table :
  public std::map< std::pair< const index_set<LO,HI>, const index_set<LO,HI> >,
                   Matrix_T* >
  {
  public:
    /*
     * @brief Single instance of basis table
     * @details
     * @return Result
     */
    static basis_table& basis() { static basis_table b; return b;}
  private:
    /*
     * @brief Friend declaration to avoid compiler warning:
     * @details
     */
    /*
     * @brief "... only defines a private destructor and has no friends"
     * @details
     */
    /*
     * @brief Ref: Carlos O'Ryan, ACE http://doc.ece.uci.edu
     * @details
     * @return Result
     */
    friend class friend_for_private_destructor;
    // Enforce singleton
    // Reference: A. Alexandrescu, "Modern C++ Design", Chapter 6
    basis_table() = default;
    ~basis_table() = default;
  public:
    basis_table(const basis_table&) = delete;
    basis_table& operator= (const basis_table&) = delete;
  };

  /*
   * @brief Create a basis element matrix within the current frame
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param ist Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::basis_matrix_t
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  basis_element(const index_set_t& ist) const
  {
    using index_set_pair_t = std::pair<const index_set_t, const index_set_t>;
    const auto& unfolded_pair = index_set_pair_t(ist, this->m_frame);

    using basis_table_t = basis_table<Scalar_T, LO, HI, basis_matrix_t>;
    auto& basis_cache = basis_table_t::basis();

    using Tuning_Values_P = typename Tune_P::tuning_values_p;
    const auto frame_count = this->m_frame.count();
    const auto use_cache = frame_count <= index_t(Tuning_Values_P::basis_max_count);

    if (use_cache)
    {
      const auto basis_it = basis_cache.find(unfolded_pair);
      if (basis_it != basis_cache.end())
        return *(basis_it->second);
    }
    const auto folded_set = ist.fold(this->m_frame);
    const auto folded_frame = this->m_frame.fold();
    const auto& folded_pair = index_set_pair_t(folded_set, folded_frame);
    using basis_pair_t = std::pair<const index_set_pair_t, basis_matrix_t *>;
    if (use_cache)
    {
      const auto basis_it = basis_cache.find(folded_pair);
      if (basis_it != basis_cache.end())
      {
        auto* result_ptr = basis_it->second;
        basis_cache.insert(basis_pair_t(unfolded_pair, result_ptr));
        return *result_ptr;
      }
    }
    const auto folded_max = folded_frame.max();
    const auto folded_min = folded_frame.min();
    const auto p = std::max(folded_max,           index_t(0));
    const auto q = std::max(index_t(-folded_min), index_t(0));
    const auto* e = (gen::generator_table<basis_matrix_t>::generator())(p, q);
    const auto dim = matrix_index_t(1) << offset_level(p, q);
    auto result = matrix::unit<basis_matrix_t>(dim);
    for (auto
        k = folded_min;
        k <= folded_max;
        ++k)
      if (folded_set[k])
        result = result * e[k];
    if (use_cache)
    {
      auto* result_ptr = new basis_matrix_t(result);
      basis_cache.insert(basis_pair_t(folded_pair, result_ptr));
      basis_cache.insert(basis_pair_t(unfolded_pair, result_ptr));
    }
    return result;
  }

  /*
   * @brief Pade' approximation
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @tparam Size
   * @param numer Value
   * @param denom Value
   * @param X Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P, const size_t Size >
  inline
  static
  auto
  pade_approx(
    const std::array<Scalar_T, Size>& numer,
    const std::array<Scalar_T, Size>& denom,
    const matrix_multi<Scalar_T,LO,HI,Tune_P>& X) -> matrix_multi<Scalar_T,LO,HI,Tune_P>
  {
    // Pade' approximation
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [GL], Section 11.3, p572-576.

    using multivector_t = matrix_multi<Scalar_T,LO,HI,Tune_P>;
    using traits_t = numeric_traits<Scalar_T>;

    if (X.isnan())
      return traits_t::NaN();

    // Array size is assumed to be even
    const auto nbr_even_powers = Size/2 - 1;

    // Create an array of even powers
    auto XX = std::vector<multivector_t>(nbr_even_powers);
    XX[0] = X * X;
    XX[1] = XX[0] * XX[0];
    for (auto
      k = size_t(2);
      k != nbr_even_powers;
      ++k)
      XX[k] = XX[k-2] * XX[1];

    // Calculate numerator N and denominator D
    auto N = multivector_t(numer[1]);
    for (auto
        k = size_t(0);
        k != nbr_even_powers;
        ++k)
      N += XX[k] * numer[2*k + 3];
    N *= X;
    N += numer[0];
    for (auto
        k = size_t(0);
        k != nbr_even_powers;
        ++k)
      N += XX[k] * numer[2*k + 2];
    auto D = multivector_t(denom[1]);
    for (auto
        k = size_t(0);
        k != nbr_even_powers;
        ++k)
      D += XX[k] * denom[2*k + 3];
    D *= X;
    D += denom[0];
    for (auto
        k = size_t(0);
        k != nbr_even_powers;
        ++k)
      D += XX[k] * denom[2*k + 2];
    return N / D;
  }

  /*
   * @brief Single step of product form of Denman-Beavers square root iteration
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param M Value
   * @param Y Value
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline
  static
  void
  db_step(matrix_multi<Scalar_T,LO,HI,Tune_P>& M, matrix_multi<Scalar_T,LO,HI,Tune_P>& Y)
  {
    // Reference: [CHKL]
    const auto& invM = M.inv();
    M = ((M + invM)/Scalar_T(2) + Scalar_T(1)) / Scalar_T(2);
    Y *= (invM + Scalar_T(1)) / Scalar_T(2);
  }

  /*
   * @brief Product form of Denman-Beavers square root iteration
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   * @param norm_tol Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline matrix_multi<Scalar_T,LO,HI,Tune_P>
  newton_schulz (
          const matrix_multi<Scalar_T,LO,HI,Tune_P>& X,
          Scalar_T norm_tol=std::pow(std::numeric_limits<Scalar_T>::epsilon(), 4))
  {
    // Reference: [CHKL]
    using multivector_t = matrix_multi<Scalar_T, LO, HI, Tune_P>;
    if (X == Scalar_T(0))
      return X;

    using Tuning_Values_P = typename Tune_P::tuning_values_p;
    static const auto sqrt_max_steps = Tuning_Values_P::db_sqrt_max_steps;
    auto M = X;
    auto Y = X;

    for (auto
        step = 0;
        step != sqrt_max_steps && (M - Scalar_T(1)).norm() > norm_tol;
        ++step)
    {
      if (Y.isnan())
        return multivector_t(numeric_traits<Scalar_T>::NaN());
      db_step(M, Y);
    }

    return Y;
  }

  /*
   * @brief Cyclic reduction square root iteration
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   * @param norm_Y_tol Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline matrix_multi<Scalar_T,LO,HI,Tune_P>
  refined_newton_schulz (
          const matrix_multi<Scalar_T,LO,HI,Tune_P>& X,
          Scalar_T norm_Y_tol=std::pow(std::numeric_limits<Scalar_T>::epsilon(), 1))
  {
    // Reference: [MB]
    using multivector_t = matrix_multi<Scalar_T, LO, HI, Tune_P>;
    if (X == Scalar_T(0))
      return X;

    using Tuning_Values_P = typename Tune_P::tuning_values_p;
    static const auto sqrt_max_steps = Tuning_Values_P::cr_sqrt_max_steps;
    auto Z = Scalar_T(2) * (Scalar_T(1) + X);
    auto Y = Scalar_T(1) - X;
    auto norm_Y = Y.norm();
    for (auto
        step = 0;
        step != sqrt_max_steps && norm_Y > norm_Y_tol;
        ++step)
    {
      const auto old_norm_Y = norm_Y;
      Y = (-Y / Z) * Y;
      norm_Y = Y.norm();
      if (Y.isnan() || (norm_Y > old_norm_Y * Scalar_T(2)))
        return multivector_t(numeric_traits<Scalar_T>::NaN());

      Z += Y * Scalar_T(2);
    }
    const auto result = Z / Scalar_T(4);

    return result;
  }

  /*
   * @brief Inverse power method for inversion
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param X Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline matrix_multi<Scalar_T,LO,HI,Tune_P>
  inverse_power_method (
    const matrix_multi<Scalar_T,LO,HI,Tune_P>& X)
  {
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [CHKL]
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [Z], Pade1
    return X.inv();
  }
}

namespace pade {
  /*
   * @brief Coefficients of numerator polynomials of Pade approximations produced by Pade1(sqrt(1+x),x,n,n)
   * @details
   */
  // Reference: [Z], Pade1
  template< typename Scalar_T >
  struct pade_sqrt_numer
  {
    using array = std::array<Scalar_T, 14>;
    static const array numer;
  };
  template< typename Scalar_T >
  const typename pade_sqrt_numer<Scalar_T>::array pade_sqrt_numer<Scalar_T>::numer =
  {
        1.0,               27.0/4.0,         81.0/4.0,       2277.0/64.0,
    10395.0/256.0,      32319.0/1024.0,    8721.0/512.0,    26163.0/4096.0,
    53703.0/32768.0,    36465.0/131072.0,  3861.0/131072.0,  7371.0/4194304.0,
      819.0/16777216.0,    27.0/67108864.0
  };

  /*
   * @brief Coefficients of denominator polynomials of Pade approximations produced by Pade1(sqrt(1+x),x,n,n)
   * @details
   */
  // Reference: [Z], Pade1
  template< typename Scalar_T >
  struct pade_sqrt_denom
  {
    using array = std::array<Scalar_T, 14>;
    static const array denom;
  };
  template< typename Scalar_T >
  const typename pade_sqrt_denom<Scalar_T>::array pade_sqrt_denom<Scalar_T>::denom =
  {
        1.0,               25.0/4.0,         69.0/4.0,       1771.0/64.0,
     7315.0/256.0,      20349.0/1024.0,    4845.0/512.0,    12597.0/4096.0,
    21879.0/32768.0,    12155.0/131072.0,  1001.0/131072.0,  1365.0/4194304.0,
       91.0/16777216.0,     1.0/67108864.0
  };

  template< >
  struct pade_sqrt_numer<float>
  {
    using array = std::array<float, 10>;
    static const array numer;
  };
  const typename pade_sqrt_numer<float>::array pade_sqrt_numer<float>::numer =
  {
       1.0,            19.0/4.0,        19.0/2.0,      665.0/64.0,
    1729.0/256.0,    2717.0/1024.0,    627.0/1024.0,   627.0/8192.0,
     285.0/65536.0,    19.0/262144.0
  };
  template< >
  struct pade_sqrt_denom<float>
  {
    using array = std::array<float, 10>;
    static const array denom;
  };
  const typename pade_sqrt_denom<float>::array pade_sqrt_denom<float>::denom =
  {
       1.0,            17.0/4.0,        15.0/2.0,      455.0/64.0,
    1001.0/256.0,    1287.0/1024.0,    231.0/1024.0,   165.0/8192.0,
      45.0/65536,       1.0/262144.0
  };

  template< >
  struct pade_sqrt_numer<long double>
  {
    using array = std::array<long double, 18>;
    static const array numer;
  };
  const typename pade_sqrt_numer<long double>::array pade_sqrt_numer<long double>::numer =
  {
        1.0L,                   35.0L/4.0L,             35.0L,               5425.0L/64.0L,
    35525.0L/256.0L,        166257.0L/1024.0L,      143325.0L/1024.0L,     740025.0L/8192.0L,
  2877875.0L/65536.0L,     4206125.0L/262144.0L,    572033.0L/131072.0L,  1820105.0L/2097152.0L,
  1028755.0L/8388608.0L,    395675.0L/33554432.0L,   24225.0L/33554432.0L,   6783.0L/268435456.0L,
     1785.0L/4294967296.0L,     35.0L/17179869184.0L
  };
  template< >
  struct pade_sqrt_denom<long double>
  {
    using array = std::array<long double, 18>;
    static const array denom;
  };
  const typename pade_sqrt_denom<long double>::array pade_sqrt_denom<long double>::denom =
  {
        1.0L,                   33.0L/4.0L,             31.0L,               4495.0L/64.0L,
    27405.0L/256.0L,        118755.0L/1024.0L,       94185.0L/1024.0L,     444015.0L/8192.0L,
  1562275.0L/65536.0L,     2042975.0L/262144.0L,    245157.0L/131072.0L,   676039.0L/2097152.0L,
   323323.0L/8388608.0L,    101745.0L/33554432.0L,    4845.0L/33554432.0L,    969.0L/268435456.0L,
      153.0L/4294967296.0L,      1.0L/17179869184.0L
  };

#if defined(_GLUCAT_USE_QD)
  template< >
  struct pade_sqrt_numer<dd_real>
  {
    using array = std::array<dd_real, 22>;
    static const array numer;
  };
  const typename pade_sqrt_numer<dd_real>::array pade_sqrt_numer<dd_real>::numer =
  {
          dd_real("1"),                               dd_real("43")/dd_real("4"),
        dd_real("215")/dd_real("4"),               dd_real("10621")/dd_real("64"),
      dd_real("90687")/dd_real("256"),            dd_real("567987")/dd_real("1024"),
     dd_real("168861")/dd_real("256"),           dd_real("1246355")/dd_real("2048"),
    dd_real("7228859")/dd_real("16384"),        dd_real("16583853")/dd_real("65536"),
    dd_real("7538115")/dd_real("65536"),       dd_real("173376645")/dd_real("4194304"),
  dd_real("195747825")/dd_real("16777216"),    dd_real("171655785")/dd_real("67108864"),
   dd_real("14375115")/dd_real("33554432"),     dd_real("14375115")/dd_real("268435456"),
   dd_real("20764055")/dd_real("4294967296"),    dd_real("5167525")/dd_real("17179869184"),
     dd_real("206701")/dd_real("17179869184"),     dd_real("76153")/dd_real("274877906944"),
       dd_real("3311")/dd_real("1099511627776") ,     dd_real("43")/dd_real("4398046511104")
  };
  template< >
  struct pade_sqrt_denom<dd_real>
  {
    using array = std::array<dd_real, 22>;
    static const array denom;
  };
  const typename pade_sqrt_denom<dd_real>::array pade_sqrt_denom<dd_real>::denom =
  {
          dd_real("1"),                               dd_real("41")/dd_real("4"),
        dd_real("195")/dd_real("4"),                dd_real("9139")/dd_real("64"),
      dd_real("73815")/dd_real("256"),            dd_real("435897")/dd_real("1024"),
     dd_real("121737")/dd_real("256"),            dd_real("840565")/dd_real("2048"),
    dd_real("4539051")/dd_real("16384"),         dd_real("9641775")/dd_real("65536"),
    dd_real("4032015")/dd_real("65536"),        dd_real("84672315")/dd_real("4194304"),
   dd_real("86493225")/dd_real("16777216"),     dd_real("67863915")/dd_real("67108864"),
    dd_real("5014575")/dd_real("33554432"),      dd_real("4345965")/dd_real("268435456"),
    dd_real("5311735")/dd_real("4294967296"),    dd_real("1081575")/dd_real("17179869184"),
      dd_real("33649")/dd_real("17179869184"),      dd_real("8855")/dd_real("274877906944"),
        dd_real("231")/dd_real("1099511627776"),       dd_real("1")/dd_real("4398046511104")
  };

  template< >
  struct pade_sqrt_numer<qd_real>
  {
    using array = std::array<qd_real, 34>;
    static const array numer;
  };
  const typename pade_sqrt_numer<qd_real>::array pade_sqrt_numer<qd_real>::numer =
  {
              qd_real("1"),                                            qd_real("67")/qd_real("4"),
            qd_real("134"),                                         qd_real("43617")/qd_real("64"),
         qd_real("633485")/qd_real("256"),                        qd_real("6992857")/qd_real("1024"),
       qd_real("15246721")/qd_real("1024"),                     qd_real("215632197")/qd_real("8192"),
     qd_real("2518145487")/qd_real("65536"),                  qd_real("12301285425")/qd_real("262144"),
     qd_real("6344873535")/qd_real("131072"),                 qd_real("89075432355")/qd_real("2097152"),
   qd_real("267226297065")/qd_real("8388608"),               qd_real("687479618945")/qd_real("33554432"),
   qd_real("379874182975")/qd_real("33554432"),             qd_real("1443521895305")/qd_real("268435456"),
  qd_real("9425348845815")/qd_real("4294967296"),          qd_real("13195488384141")/qd_real("17179869184"),
   qd_real("987417498133")/qd_real("4294967296"),           qd_real("8055248011085")/qd_real("137438953472"),
  qd_real("6958363175533")/qd_real("549755813888"),         qd_real("5056698705201")/qd_real("2199023255552"),
   qd_real("766166470485")/qd_real("2199023255552"),         qd_real("766166470485")/qd_real("17592186044416"),
   qd_real("623623871325")/qd_real("140737488355328"),       qd_real("203123203803")/qd_real("562949953421312"),
     qd_real("6478601247")/qd_real("281474976710656"),         qd_real("5038912081")/qd_real("4503599627370496"),
      qd_real("719844583")/qd_real("18014398509481984"),         qd_real("71853815")/qd_real("72057594037927936"),
        qd_real("1165197")/qd_real("72057594037927936"),            qd_real("87703")/qd_real("576460752303423488"),
          qd_real("12529")/qd_real("18446744073709551616"),            qd_real("67")/qd_real("73786976294838206464")
  };
  template< >
  struct pade_sqrt_denom<qd_real>
  {
    using array = std::array<qd_real, 34>;
    static const array denom;
  };
  const typename pade_sqrt_denom<qd_real>::array pade_sqrt_denom<qd_real>::denom =
  {
            qd_real("1"),                                              qd_real("65")/qd_real("4"),
          qd_real("126"),                                           qd_real("39711")/qd_real("64"),
       qd_real("557845")/qd_real("256"),                          qd_real("5949147")/qd_real("1024"),
     qd_real("12515965")/qd_real("1024"),                       qd_real("170574723")/qd_real("8192"),
   qd_real("1916797311")/qd_real("65536"),                     qd_real("8996462475")/qd_real("262144"),
   qd_real("4450881435")/qd_real("131072"),                   qd_real("59826782925")/qd_real("2097152"),
 qd_real("171503444385")/qd_real("8388608"),                 qd_real("420696483235")/qd_real("33554432"),
 qd_real("221120793075")/qd_real("33554432"),                qd_real("797168807855")/qd_real("268435456"),
qd_real("4923689695575")/qd_real("4294967296"),             qd_real("6499270398159")/qd_real("17179869184"),
 qd_real("456864812569")/qd_real("4294967296"),             qd_real("3486599885395")/qd_real("137438953472"),
qd_real("2804116503573")/qd_real("549755813888"),           qd_real("1886827875075")/qd_real("2199023255552"),
 qd_real("263012370465")/qd_real("2199023255552"),           qd_real("240141729555")/qd_real("17592186044416"),
 qd_real("176848560525")/qd_real("140737488355328"),          qd_real("51538723353")/qd_real("562949953421312"),
   qd_real("1450433115")/qd_real("281474976710656"),            qd_real("977699359")/qd_real("4503599627370496"),
    qd_real("118183439")/qd_real("18014398509481984"),            qd_real("9652005")/qd_real("72057594037927936"),
       qd_real("121737")/qd_real("72057594037927936"),               qd_real("6545")/qd_real("576460752303423488"),
          qd_real("561")/qd_real("18446744073709551616"),               qd_real("1")/qd_real("73786976294838206464")
  };
#endif
}

namespace glucat
{
  /*
   * @brief Square root of multivector with specified complexifier
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   * @param i Complexifier -- commuting sqrt of -1
   * @param level Recursion level
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  auto
  matrix_sqrt(const matrix_multi<Scalar_T,LO,HI,Tune_P>& val,
              const matrix_multi<Scalar_T,LO,HI,Tune_P>& i,
              const index_t level) -> matrix_multi<Scalar_T,LO,HI,Tune_P>
  {
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [Z], Pade1

    using traits_t = numeric_traits<Scalar_T>;

    if (val.isnan())
      return traits_t::NaN();

    const auto scr_val = val.scalar();
    if (val == scr_val)
    {
      if (scr_val < Scalar_T(0))
        return i * traits_t::sqrt(-scr_val);
      else
        return traits_t::sqrt(scr_val);
    }

    // Scale val towards abs(A) == 1 or towards A == 1 as appropriate
    const auto scale =
      (scr_val != Scalar_T(0) && (val/scr_val - Scalar_T(1)).norm() < Scalar_T(1))
      ? scr_val
      : (scr_val < Scalar_T(0))
        ? -abs(val)
        :  abs(val);
    const auto sqrt_scale = traits_t::sqrt(traits_t::abs(scale));
    if (traits_t::isNaN_or_isInf(sqrt_scale))
      return traits_t::NaN();

    using multivector_t = matrix_multi<Scalar_T,LO,HI,Tune_P>;
    auto rescale = multivector_t(sqrt_scale);
    if (scale < Scalar_T(0))
      rescale = i * sqrt_scale;

    const auto& unitval = val / scale;
    static const auto max_norm = Scalar_T(1.0/4.0);
    auto use_approx_sqrt = true;
    auto use_cr_sqrt = false;
    auto scaled_result = multivector_t();
    static const auto sqrt_2 = traits_t::sqrt(Scalar_T(2));
    if (level == 0)
    {
      // What kind of eigenvalues does the matrix contain?
      const auto genus = unitval.m_matrix.classify_eigenvalues();
      const index_t next_level =
        (genus.m_is_singular)
        ? level
        : level + 1;
      switch (genus.m_eig_case)
      {
      case matrix::neg_real_eigs:
        scaled_result = matrix_sqrt(-i * unitval, i, next_level) * (i + Scalar_T(1)) / sqrt_2;
        use_approx_sqrt = false;
        break;
      case matrix::both_eigs:
        {
          const auto safe_arg = genus.m_safe_arg;
          scaled_result = matrix_sqrt(exp(i*safe_arg) * unitval, i, next_level) * exp(-i*safe_arg / Scalar_T(2));
        }
        use_approx_sqrt = false;
        break;
      default:
        break;
      }
      use_cr_sqrt = genus.m_is_singular;
    }
    if (use_approx_sqrt)
    {
      bool used_db_sqrt = false;
      if ((unitval - Scalar_T(1)).norm() < max_norm)
      {
         // Pade' approximation of square root
         scaled_result = pade_approx(pade::pade_sqrt_numer<Scalar_T>::numer,
                      pade::pade_sqrt_denom<Scalar_T>::denom,
                      unitval - Scalar_T(1));
      }
      else if (use_cr_sqrt)
      {
          scaled_result = refined_newton_schulz(unitval);
      }
      else
      {
          scaled_result = newton_schulz(unitval);
          used_db_sqrt = true;
      }

      if (used_db_sqrt && (scaled_result.isnan() || !approx_equal(pow(scaled_result, 2), unitval)))
      {
          // Fallback to Cyclic Reduction if Denman-Beavers failed or produced poor result
          scaled_result = refined_newton_schulz(unitval);
      }
    }
    const auto frame = (use_approx_sqrt)
      ? val.frame()
      : i.frame();

    if (scaled_result.isnan()) return traits_t::NaN();
    // Return result even if verification fails, as long as it's not NaN.
    return multivector_t(scaled_result * rescale, frame);
  }

  /*
   * @brief Square root of multivector with specified complexifier
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   * @param i Complexifier -- commuting sqrt of -1
   * @param prechecked Already checked?
   * @return True if successful or condition met
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline matrix_multi<Scalar_T,LO,HI,Tune_P>
  sqrt(const matrix_multi<Scalar_T,LO,HI,Tune_P>& val, const matrix_multi<Scalar_T,LO,HI,Tune_P>& i, bool prechecked)
  {
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [Z], Pade1

    using traits_t = numeric_traits<Scalar_T>;

    if (val.isnan())
      return traits_t::NaN();

    using tune_same_p = typename Tune_P::tuning_same_p;
    const precision_t function_precision = Tune_P::function_precision;
    check_complex(val, i, prechecked);

    if constexpr (function_precision == precision_demoted)
    {
      using demoted_scalar_t = typename traits_t::demoted::type;
      using demoted_multivector_t = matrix_multi<demoted_scalar_t,LO,HI,tune_same_p>;

      const auto& demoted_val = demoted_multivector_t(val);
      const auto& demoted_i = demoted_multivector_t(i);

      return matrix_multi<Scalar_T,LO,HI,Tune_P>(matrix_sqrt(demoted_val, demoted_i, 0));
    }
    else if constexpr (function_precision == precision_promoted)
    {
      using promoted_scalar_t = typename traits_t::promoted::type;
      using promoted_multivector_t = matrix_multi<promoted_scalar_t,LO,HI,tune_same_p>;

      const auto& promoted_val = promoted_multivector_t(val);
      const auto& promoted_i = promoted_multivector_t(i);

      return matrix_multi<Scalar_T,LO,HI,Tune_P>(matrix_sqrt(promoted_val, promoted_i, 0));
    }
    else
      return matrix_sqrt(val, i, 0);
  }
}

namespace pade {
  /*
   * @brief Coefficients of numerator polynomials of Pade approximations produced by Pade1(log(1+x),x,n,n)
   * @details
   */
  // Reference: [Z], Pade1
  template< typename Scalar_T >
  struct pade_log_numer
  {
    using array = std::array<Scalar_T, 14>;
    static const array numer;
  };
  template< typename Scalar_T >
  const typename pade_log_numer<Scalar_T>::array pade_log_numer<Scalar_T>::numer =
  {
         0.0,                     1.0,                    6.0,            4741.0/300.0,
      1441.0/60.0,           107091.0/4600.0,          8638.0/575.0,    263111.0/40250.0,
    153081.0/80500.0,        395243.0/1101240.0,      28549.0/688275.0, 605453.0/228813200.0,
    785633.0/10296594000.0, 1145993.0/1873980108000.0
  };

  /*
   * @brief Coefficients of denominator polynomials of Pade approximations produced by Pade1(log(1+x),x,n,n)
   * @details
   */
  // Reference: [Z], Pade1
  template< typename Scalar_T >
  struct pade_log_denom
  {
    using array = std::array<Scalar_T, 14>;
    static const array denom;
  };
  template< typename Scalar_T >
  const typename pade_log_denom<Scalar_T>::array pade_log_denom<Scalar_T>::denom =
  {
         1.0,                    13.0/2.0,             468.0/25.0,        1573.0/50.0,
      1573.0/46.0,            11583.0/460.0,         10296.0/805.0,       2574.0/575.0,
     11583.0/10925.0,           143.0/874.0,           572.0/37145.0,      117.0/148580.0,
        13.0/742900.0,            1.0/10400600.0
  };

  template< >
  struct pade_log_numer<float>
  {
    using array = std::array<float, 10>;
    static const array numer;
  };
  const typename pade_log_numer<float>::array pade_log_numer<float>::numer =
  {
      0.0,            1.0,             4.0,       1337.0/204.0,
    385.0/68.0,    1879.0/680.0,     193.0/255.0,  197.0/1820.0,
    419.0/61880.0, 7129.0/61261200.0
  };
  template< >
  struct pade_log_denom<float>
  {
    using array = std::array<float, 10>;
    static const array denom;
  };
  const typename pade_log_denom<float>::array pade_log_denom<float>::denom =
  {
      1.0,            9.0/2.0,       144.0/17.0,   147.0/17.0,
    441.0/85.0,      63.0/34.0,       84.0/221.0,    9.0/221.0,
      9.0/4862.0,     1.0/48620.0
  };

  template< >
  struct pade_log_numer<long double>
  {
    using array = std::array<long double, 18>;
    static const array numer;
  };
  const typename pade_log_numer<long double>::array pade_log_numer<long double>::numer =
  {
         0.0L,                       1.0L,                           8.0L,                    3835.0L/132.0L,
      8365.0L/132.0L,         11363807.0L/122760.0L,            162981.0L/1705.0L,         9036157.0L/125860.0L,
  18009875.0L/453096.0L,      44211925.0L/2718576.0L,          4149566.0L/849555.0L,      16973929.0L/16020180.0L,
    172459.0L/1068012.0L,    116317061.0L/7025382936.0L,      19679783.0L/18441630207.0L, 23763863.0L/614721006900.0L,
     50747.0L/79318839600.0L, 42142223.0L/14295951736466400.0L
  };
  template< >
  struct pade_log_denom<long double>
  {
    using array = std::array<long double, 18>;
    static const array denom;
  };
  const typename pade_log_denom<long double>::array pade_log_denom<long double>::denom =
  {
         1.0L,                      17.0L/2.0L,                   1088.0L/33.0L,               850.0L/11.0L,
     41650.0L/341.0L,           140777.0L/1023.0L,             1126216.0L/9889.0L,           63206.0L/899.0L,
    790075.0L/24273.0L,          60775.0L/5394.0L,               38896.0L/13485.0L,          21658.0L/40455.0L,
     21658.0L/310155.0L,          4165.0L/682341.0L,               680.0L/2047023.0L,           34.0L/3411705.0L,
        17.0L/129644790.0L,          1.0L/2333606220
  };
#if defined(_GLUCAT_USE_QD)
  template< >
  struct pade_log_numer<dd_real>
  {
    using array = std::array<dd_real, 22>;
    static const array numer;
  };
  const typename pade_log_numer<dd_real>::array pade_log_numer<dd_real>::numer =
  {
          dd_real("0"),                                  dd_real("1"),
         dd_real("10"),                              dd_real("22781")/dd_real("492"),
      dd_real("21603")/dd_real("164"),             dd_real("5492649")/dd_real("21320"),
     dd_real("978724")/dd_real("2665"),            dd_real("4191605")/dd_real("10619"),
   dd_real("12874933")/dd_real("39442"),          dd_real("11473457")/dd_real("54612"),
    dd_real("2406734")/dd_real("22755"),         dd_real("166770367")/dd_real("4004880"),
   dd_real("30653165")/dd_real("2402928"),       dd_real("647746389")/dd_real("215195552"),
   dd_real("25346331")/dd_real("47074027"),      dd_real("278270613")/dd_real("3900419380"),
  dd_real("105689791")/dd_real("15601677520"),   dd_real("606046475")/dd_real("1379188292768"),
     dd_real("969715")/dd_real("53502994116"),    dd_real("11098301")/dd_real("26204577562592"),
     dd_real("118999")/dd_real("26204577562592"), dd_real("18858053")/dd_real("1392249205900512960")
  };
  template< >
  struct pade_log_denom<dd_real>
  {
    using array = std::array<dd_real, 22>;
    static const array denom;
  };
  const typename pade_log_denom<dd_real>::array pade_log_denom<dd_real>::denom =
  {
          dd_real("1"),                                 dd_real("21")/dd_real("2"),
       dd_real("2100")/dd_real("41"),                dd_real("12635")/dd_real("82"),
     dd_real("341145")/dd_real("1066"),            dd_real("1037799")/dd_real("2132"),
   dd_real("11069856")/dd_real("19721"),           dd_real("9883800")/dd_real("19721"),
    dd_real("6918660")/dd_real("19721"),            dd_real("293930")/dd_real("1517"),
    dd_real("1410864")/dd_real("16687"),             dd_real("88179")/dd_real("3034"),
     dd_real("734825")/dd_real("94054"),            dd_real("305235")/dd_real("188108"),
     dd_real("348840")/dd_real("1363783"),           dd_real("40698")/dd_real("1363783"),
       dd_real("6783")/dd_real("2727566"),            dd_real("9975")/dd_real("70916716"),
        dd_real("266")/dd_real("53187537"),              dd_real("7")/dd_real("70916716"),
          dd_real("7")/dd_real("8155422340"),            dd_real("1")/dd_real("538257874440")
  };

  template< >
  struct pade_log_numer<qd_real>
  {
    using array = std::array<qd_real, 34>;
    static const array numer;
  };
  const typename pade_log_numer<qd_real>::array pade_log_numer<qd_real>::numer =
  {
                qd_real("0"),                                                          qd_real("1"),
                qd_real("16"),                                                     qd_real("95201")/qd_real("780"),
             qd_real("30721")/qd_real("52"),                                     qd_real("7416257")/qd_real("3640"),
           qd_real("1039099")/qd_real("195"),                                 qd_real("6097772319")/qd_real("555100"),
        qd_real("1564058073")/qd_real("85400"),                              qd_real("30404640205")/qd_real("1209264"),
         qd_real("725351278")/qd_real("25193"),                            qd_real("4092322670789")/qd_real("147429436"),
     qd_real("4559713849589")/qd_real("201040140"),                        qd_real("5049361751189")/qd_real("320023080"),
       qd_real("74979677195")/qd_real("8000577"),                         qd_real("16569850691873")/qd_real("3481514244"),
     qd_real("1065906022369")/qd_real("515779888"),                      qd_real("335956770855841")/qd_real("438412904800"),
  qd_real("1462444287585964")/qd_real("6041877844275"),                  qd_real("397242326339851")/qd_real("6122436215532"),
    qd_real("64211291334131")/qd_real("4373168725380"),                  qd_real("142322343550859")/qd_real("51080680851480"),
   qd_real("154355972958659")/qd_real("351179680853925"),                qd_real("167483568676259")/qd_real("2937139148960100"),
     qd_real("4230788929433")/qd_real("704913395750424"),                qd_real("197968763176019")/qd_real("392923948371995600"),
    qd_real("10537522306718")/qd_real("319250708052246425"),             qd_real("236648286272519")/qd_real("144249197475035425500"),
   qd_real("260715545088119")/qd_real("4375558990076074573500"),         qd_real("289596255666839")/qd_real("192874640282553367199880"),
     qd_real("8802625510547")/qd_real("361639950529787563499775"),       qd_real("373831661521439")/qd_real("1659204093030665341336967700"),
   qd_real("446033437968239")/qd_real("464577146048586295574350956000"),  qd_real("53676090078349")/qd_real("47386868896955802148583797512000")
      };
  template< >
  struct pade_log_denom<qd_real>
  {
    using array = std::array<qd_real, 34>;
    static const array denom;
  };
  const typename pade_log_denom<qd_real>::array pade_log_denom<qd_real>::denom =
  {
                 qd_real("1"),                                                        qd_real("33")/qd_real("2"),
              qd_real("8448")/qd_real("65"),                                       qd_real("42284")/qd_real("65"),
            qd_real("211420")/qd_real("91"),                                      qd_real("573562")/qd_real("91"),
          qd_real("32119472")/qd_real("2379"),                                  qd_real("92917044")/qd_real("3965"),
         qd_real("603960786")/qd_real("17995"),                                qd_real("144626625")/qd_real("3599"),
        qd_real("2776831200")/qd_real("68381"),                              qd_real("16692542100")/qd_real("478667"),
       qd_real("12241197540")/qd_real("478667"),                              qd_real("1098569010")/qd_real("68381"),
       qd_real("31387686000")/qd_real("3624193"),                             qd_real("9939433900")/qd_real("2479711"),
       qd_real("67091178825")/qd_real("42155087"),                            qd_real("2683647153")/qd_real("4959422"),
       qd_real("19083713088")/qd_real("121505839"),                           qd_real("4708152900")/qd_real("121505839"),
         qd_real("941630580")/qd_real("116546417"),                             qd_real("88704330")/qd_real("62755763"),
          qd_real("12902448")/qd_real("62755763"),                               qd_real("1542684")/qd_real("62755763"),
           qd_real("6427850")/qd_real("2698497809"),                             qd_real("3471039")/qd_real("18889484663"),
           qd_real("8544096")/qd_real("774468871183"),                             qd_real("39556")/qd_real("79027435835"),
            qd_real("118668")/qd_real("7191496660985"),                            qd_real("10230")/qd_real("27327687311743"),
              qd_real("5456")/qd_real("1011124430534491"),                            qd_real("44")/qd_real("1011124430534491"),
                qd_real("11")/qd_real("70778710137414370"),                            qd_real("1")/qd_real("7219428434016265740")
  };
#endif
}

namespace glucat{
  /*
   * @brief Pade' approximation of log
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>
  pade_log(const matrix_multi<Scalar_T,LO,HI,Tune_P>& val)
  {
    // Reference: [GW], Section 4.3, pp318-322
    // Reference: [CHKL]
    // Reference: [GL], Section 11.3, p572-576
    // Reference: [Z], Pade1

    using traits_t = numeric_traits<Scalar_T>;
    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();
    else
    {
      const auto result = pade_approx(
                         pade::pade_log_numer<Scalar_T>::numer,
                         pade::pade_log_denom<Scalar_T>::denom,
                         val - Scalar_T(1));

      return result;
    }
  }

  /*
   * @brief Incomplete square root cascade and Pade' approximation of log
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline matrix_multi<Scalar_T,LO,HI,Tune_P>
  cascade_log(const matrix_multi<Scalar_T,LO,HI,Tune_P>& val)
  {

    // Reference: [CHKL]
    using multivector_t = matrix_multi<Scalar_T,LO,HI,Tune_P>;
    using traits_t = numeric_traits<Scalar_T>;
    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();

    using limits_t = std::numeric_limits<Scalar_T>;
    using Tuning_Values_P = typename Tune_P::tuning_values_p;
    static const auto epsilon = limits_t::epsilon();
    static const auto max_inner_norm = traits_t::pow(epsilon, 2);
    static const auto max_outer_norm = Scalar_T(6.0/limits_t::digits);
    auto Y = val;
    auto E = multivector_t(Scalar_T(0));
    Scalar_T norm_Y_1;
    auto pow_2_outer_step = Scalar_T(1);
    auto pow_4_outer_step = Scalar_T(1);
    int outer_step;

    for (outer_step = 0, norm_Y_1 = (Y - Scalar_T(1)).norm();
        outer_step != Tuning_Values_P::log_max_outer_steps && norm_Y_1 * pow_2_outer_step > max_outer_norm;
        ++outer_step,    norm_Y_1 = (Y - Scalar_T(1)).norm())
    {
      if (Y == Scalar_T(0) || Y.isnan())
        return traits_t::NaN();

      // Incomplete product form of Denman-Beavers square root iteration
      auto M = Y;
      for (auto
          inner_step = 0;
          inner_step != Tuning_Values_P::log_max_inner_steps &&
            (M - Scalar_T(1)).norm() * pow_4_outer_step > max_inner_norm;
          ++inner_step)
        db_step(M, Y);

      E += (M - Scalar_T(1)) * pow_2_outer_step;
      pow_2_outer_step *= Scalar_T(2);
      pow_4_outer_step *= Scalar_T(4);
    }
    if (outer_step == Tuning_Values_P::log_max_outer_steps && norm_Y_1 * pow_2_outer_step > max_outer_norm)
      return traits_t::NaN();
    else
      return pade_log(Y) * pow_2_outer_step - E;
  }

  /*
   * @brief Natural logarithm of multivector with specified complexifier
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   * @param i Complexifier -- commuting sqrt of -1
   * @param level Recursion level
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  matrix_multi<Scalar_T,LO,HI,Tune_P>
  matrix_log( const matrix_multi<Scalar_T,LO,HI,Tune_P>& val,
              const matrix_multi<Scalar_T,LO,HI,Tune_P>& i,
              const index_t level)
  {
    // Scaled incomplete square root cascade and scaled Pade' approximation of log
    // Reference: [CHKL]

    using traits_t = numeric_traits<Scalar_T>;
    if (val == Scalar_T(0) || val.isnan() || val.inv().isnan())
      return traits_t::NaN();

    static const auto pi = traits_t::pi();
    const auto scr_val = val.scalar();
    if (val == scr_val)
    {
      if (scr_val < Scalar_T(0))
        return i * pi + traits_t::log(-scr_val);
      else
        return traits_t::log(scr_val);
    }

    // Scale val towards abs(A) == 1 or towards A == 1 as appropriate
    const auto max_norm = Scalar_T(1.0/9.0);
    const auto scale =
      (scr_val != Scalar_T(0) && (val/scr_val - Scalar_T(1)).norm() < max_norm)
      ? scr_val
      : (scr_val < Scalar_T(0))
        ? -abs(val)
        :  abs(val);
    if (scale == Scalar_T(0))
      return traits_t::NaN();

    using multivector_t = matrix_multi<Scalar_T,LO,HI,Tune_P>;
    const auto log_scale = traits_t::log(traits_t::abs(scale));
    auto rescale = multivector_t(log_scale);
    if (scale < Scalar_T(0))
      rescale = i * pi + log_scale;
    const auto unitval = val/scale;
    if (unitval.inv().isnan())
      return traits_t::NaN();

    auto use_approx_log = true;
    auto scaled_result = multivector_t();
    if (level == 0)
    {
      // What kind of eigenvalues does the matrix contain?
      auto genus = unitval.m_matrix.classify_eigenvalues();
      if (genus.m_is_singular)
        return traits_t::NaN();

      switch (genus.m_eig_case)
      {
        case matrix::neg_real_eigs:
        scaled_result = matrix_log(-i * unitval, i, level + 1) + i * pi/Scalar_T(2);
        use_approx_log = false;
        break;
        case matrix::both_eigs:
        {
          const Scalar_T safe_arg = genus.m_safe_arg;
          scaled_result = matrix_log(exp(i*safe_arg) * unitval, i, level + 1) - i * safe_arg;
        }
        use_approx_log = false;
        break;
      default:
        scaled_result = cascade_log(unitval);
        break;
      }
    }
    else
      scaled_result = cascade_log(unitval);
    const auto frame = (use_approx_log)
    ? val.frame()
    : i.frame();
    return (scaled_result.isnan())
      ? multivector_t(traits_t::NaN())
      : multivector_t(scaled_result + rescale, frame);
  }

  /*
   * @brief Natural logarithm of multivector with specified complexifier
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   * @param i Complexifier -- commuting sqrt of -1
   * @param prechecked Already checked?
   * @return True if successful or condition met
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline matrix_multi<Scalar_T,LO,HI,Tune_P>
  log(const matrix_multi<Scalar_T,LO,HI,Tune_P>& val, const matrix_multi<Scalar_T,LO,HI,Tune_P>& i, bool prechecked)
  {
    using traits_t = numeric_traits<Scalar_T>;

    if (val == Scalar_T(0) || val.isnan())
      return traits_t::NaN();

    check_complex(val, i, prechecked);

    using tune_same_p = typename Tune_P::tuning_same_p;
    const precision_t function_precision = Tune_P::function_precision;

    if constexpr (function_precision == precision_demoted)
    {
      using demoted_scalar_t = typename traits_t::demoted::type;
      using demoted_multivector_t = matrix_multi<demoted_scalar_t,LO,HI,tune_same_p>;

      const auto& demoted_val = demoted_multivector_t(val);
      const auto& demoted_i = demoted_multivector_t(i);

      return matrix_multi<Scalar_T,LO,HI,Tune_P>(matrix_log(demoted_val, demoted_i, 0));
    }
    else if constexpr (function_precision == precision_promoted)
    {
      using promoted_scalar_t = typename traits_t::promoted::type;
      using promoted_multivector_t = matrix_multi<promoted_scalar_t,LO,HI,tune_same_p>;

      const auto& promoted_val = promoted_multivector_t(val);
      const auto& promoted_i = promoted_multivector_t(i);

      return matrix_multi<Scalar_T,LO,HI,Tune_P>(matrix_log(promoted_val, promoted_i, 0));
    }
    else
      return matrix_log(val, i, 0);
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
  inline typename matrix_multi<Scalar_T,LO,HI,Tune_P>::size_type
  matrix_multi<Scalar_T,LO,HI,Tune_P>::
  nbr_terms() const
  {
    return framed_multi_t(*this).nbr_terms();
  }

  /*
   * @brief Exponential of multivector
   * @details
   * @tparam Scalar_T
   * @tparam LO
   * @tparam HI
   * @tparam Tune_P
   * @param val Value
   * @return Result
   */
  template< typename Scalar_T, const index_t LO, const index_t HI, typename Tune_P >
  inline matrix_multi<Scalar_T,LO,HI,Tune_P>
  exp(const matrix_multi<Scalar_T,LO,HI,Tune_P>& val)
  {
    using traits_t = numeric_traits<Scalar_T>;
    if (val.isnan())
      return traits_t::NaN();

    const auto scr_val = val.scalar();
    if (val == scr_val)
      return traits_t::exp(scr_val);

    using tune_same_p = typename Tune_P::tuning_same_p;
    const precision_t function_precision = Tune_P::function_precision;

    if constexpr (function_precision == precision_demoted)
    {
      using demoted_scalar_t = typename traits_t::demoted::type;
      using demoted_multivector_t = matrix_multi<demoted_scalar_t,LO,HI,tune_same_p>;

      const auto& demoted_val = demoted_multivector_t(val);
      return matrix_multi<Scalar_T,LO,HI,Tune_P>(clifford_exp(demoted_val));
    }
    else if constexpr (function_precision == precision_promoted)
    {
      using promoted_scalar_t = typename traits_t::promoted::type;
      using promoted_multivector_t = matrix_multi<promoted_scalar_t,LO,HI,tune_same_p>;

      const auto& promoted_val = promoted_multivector_t(val);
      return matrix_multi<Scalar_T,LO,HI,Tune_P>(clifford_exp(promoted_val));
    }
    else
      return clifford_exp(val);
  }
}
#ifdef GLUCAT_DOCTEST
#include <iostream>
#include <sstream>
#include <cmath>
#include <numbers>
#include <chrono>

TEST_CASE("matrix_multi<Scalar_T, LO, HI, Tune_P>") {
  using namespace glucat;
  using mm_t = glucat::matrix_multi<double, -8, 8>;
  using T = double;
  using is_t = mm_t::index_set_t;

  SUBCASE("Metadata") {
    CHECK(mm_t::classname() == "matrix_multi");
  }

  SUBCASE("Constructor and string representation") {
    mm_t m1(T(2.0), mm_t::index_set_t());
    std::ostringstream oss1;
    oss1 << m1;
    CHECK(oss1.str() == "2");

    mm_t m2("2{1,2,3}");
    std::ostringstream oss2;
    oss2 << m2;
    CHECK(oss2.str() == "2{1,2,3}");

    mm_t m3("-{1}");
    std::ostringstream oss3;
    oss3 << m3;
    CHECK(oss3.str() == "-{1}");
  }

  SUBCASE("Geometric operations") {
    mm_t e1("{1}");
    mm_t e2("{2}");
    mm_t e3("{3}");
    CHECK((e1 * e2) == mm_t("{1,2}"));
    CHECK((e2 * e1) == mm_t("-{1,2}"));
    CHECK((e1 * e1) == mm_t(T(1.0), mm_t::index_set_t()));
    CHECK(scalar(e1 * e2) == T(0.0));
    CHECK(scalar(e1 * e1) == T(1.0));

    // HS (1.25a): (a ^ b) ^ c == a ^ (b ^ c)
    CHECK(((e1 ^ e2) ^ e3) == (e1 ^ (e2 ^ e3)));

    // HS (1.31): a_1 * b == (a_1 & b) + (a_1 ^ b)
    mm_t b = mm_t("1+{1}+{2}+{1,2}");
    CHECK((e1 * b) == ((e1 & b) + (e1 ^ b)));

    // HS (1.44): star(a, b) == scalar(a * b)
    CHECK(star(e1, e2) == scalar(e1 * e2));
    CHECK(star(e1, e1) == scalar(e1 * e1));
  }

  SUBCASE("Arithmetic and approximate equality") {
    mm_t m1(T(1.0));
    mm_t m2("{1}");
    CHECK((m1 + m2) == mm_t("1+{1}"));
    CHECK((m1 - m2) == mm_t("1-{1}"));
    CHECK((m1 * T(2.0)) == mm_t(T(2.0)));
    CHECK((T(2.0) * m1) == mm_t(T(2.0)));
    CHECK(approx_equal(m1, mm_t(T(1.0) + T(1e-15))));
    CHECK(m1 != m2);
    CHECK(m1 != 0.0);
    CHECK(0.0 != m1);
  }

  SUBCASE("Transcendental functions") {
    mm_t x("{1,2}");
    const T pi = T(std::numbers::pi);

    // exp and log
    mm_t e_x = exp(x * (pi/4.0));
    // exp({1,2}*pi/4) = cos(pi/4) + {1,2}*sin(pi/4) = (1 + {1,2})/sqrt(2)
    CHECK(approx_equal(e_x, (mm_t(1.0) + mm_t("{1,2}")) / std::sqrt(2.0)));
    CHECK(approx_equal(exp(log(e_x)), e_x));

    // sin and cos
    mm_t s_x = sin(x);
    mm_t c_x = cos(x);
    // sin^2 + cos^2 = 1
    CHECK(approx_equal(s_x*s_x + c_x*c_x, mm_t(1.0)));
    CHECK(approx_equal(cos(acos(mm_t(0.5))), mm_t(0.5)));

    // sinh and cosh
    mm_t sh_x = sinh(x);
    mm_t ch_x = cosh(x);
    // cosh^2 - sinh^2 = 1
    CHECK(approx_equal(ch_x*ch_x - sh_x*sh_x, mm_t(1.0)));
    CHECK(approx_equal(cosh(acosh(mm_t(2.0))), mm_t(2.0)));

    // tan and tanh
    CHECK(approx_equal(tan(atan(mm_t(0.5))), mm_t(0.5)));
    CHECK(approx_equal(tanh(atanh(mm_t(0.5))), mm_t(0.5)));

    // peg11.h identities
    mm_t A = mm_t("0.5{1}+0.5{1,2}");
    CHECK(approx_equal(exp(A) * exp(-A), mm_t(T(1.0), mm_t::index_set_t())));
    CHECK(approx_equal(cosh(A) + sinh(A), exp(A)));
    CHECK(approx_equal(cos(A) + complexifier(A)*sin(A), exp(complexifier(A)*A)));
    CHECK(approx_equal(cos(A)*tan(A), sin(A)));
    CHECK(approx_equal(cosh(A)*tanh(A), sinh(A)));
    CHECK(approx_equal(sqrt(mm_t(T(4.0))), mm_t(T(2.0))));
  }

  SUBCASE("Adversarial and edge cases") {
    // NaN handling
    mm_t n(glucat::numeric_traits<T>::NaN());
    CHECK(n.isnan());
    CHECK(exp(n).isnan());
    CHECK(log(n).isnan());
    CHECK(log(n, mm_t(glucat::numeric_traits<T>::NaN()), false).isnan());

    // Purely scalar multivector in transcendental functions
    mm_t s(T(2.0), mm_t::index_set_t());
    CHECK(approx_equal(exp(s), mm_t(T(std::exp(2.0)))));
    CHECK(approx_equal(log(s), mm_t(T(std::log(2.0)))));
    CHECK(approx_equal(sqrt(s), mm_t(T(std::sqrt(2.0)))));

    // Zero
    mm_t zero(T(0.0), mm_t::index_set_t());
    // Log of zero returns NaN in GluCat
    CHECK(log(zero).isnan());

    // Negative log with complexifier
    mm_t neg(T(-1.0), mm_t::index_set_t());
    mm_t i("{-1}");
    CHECK(approx_equal(log(neg, i, false), i * T(std::numbers::pi)));
  }
  SUBCASE("Transcendental identities (random)") {
    using namespace glucat;
    using index_set_t = mm_t::index_set_t; index_set_t frm = index_set_t();
    const T fill = T(0.5);
    for (index_t i = 1; i <= 7; ++i) {
      frm |= index_set_t(i);
      frm |= index_set_t(-i);
      
      mm_t a = mm_t::random(frm, fill);
      
      // exp(a) * exp(-a) == 1
      CHECK(approx_equal(exp(a) * exp(-a), mm_t(T(1.0), mm_t::index_set_t())));
      
      // cosh(a) + sinh(a) == exp(a)
      CHECK(approx_equal(cosh(a) + sinh(a), exp(a)));
      
      // sqrt(a) * sqrt(a) == a
      mm_t s = sqrt(a);
      if (!s.isnan() && !s.isinf())
        CHECK(approx_equal(s * s, a));
    }
  }

  SUBCASE("Geometric algebra identities (random)") {
    using index_set_t = mm_t::index_set_t; index_set_t frm = index_set_t();
    const T fill = T(0.5);
    for (index_t i = 1; i <= 7; ++i) {
      frm |= index_set_t(i);
      frm |= index_set_t(-i);
      
      mm_t a = mm_t::random(frm, fill);
      mm_t b = mm_t::random(frm, fill);
      mm_t c = mm_t::random(frm, fill);
      
      // [HS] (1.25a): (a ^ b) ^ c == a ^ (b ^ c)
      CHECK(approx_equal((a ^ b) ^ c, a ^ (b ^ c)));
      
      // [HS] (1.31): a_1 * b == (a_1 & b) + (a_1 ^ b)
      mm_t a_1 = a(1);
      CHECK(approx_equal(a_1 * b, (a_1 & b) + (a_1 ^ b)));
      
      // [HS] (1.44): star(a, b) == scalar(a * b)
      CHECK(star(a, b) == doctest::Approx(scalar(a * b)));
      
      // [HS] (1.21a): (a_r * b_s)(|r-s|) == a_r & b_s (sampled)
      index_t r = frm.count() / 2;
      index_t s = frm.count() / 2;
      mm_t a_r = a(r);
      mm_t b_s = b(s);
      CHECK(approx_equal(a_r & b_s, (a_r * b_s)(index_t(std::abs(r-s)))));
    }
  }

  SUBCASE("Mixed-precision and assignment interchangeability") {
    using mm_f_t = glucat::matrix_multi<float, -8, 8>;
    using mm_d_t = glucat::matrix_multi<double, -8, 8>;

    mm_f_t m_f1("{1}");
    mm_f_t m_f2("{2}");
    mm_d_t m_d1("{1}");
    mm_d_t m_d2("{2}");

    // Compound assignment (matching scalar)
    m_f1 = mm_f_t("{1}");
    m_f1 += m_f2;
    CHECK(m_f1 == mm_f_t("{1}+{2}"));

    m_f1 = mm_f_t("{1}");
    m_f1 -= m_f2;
    CHECK(m_f1 == mm_f_t("{1}-{2}"));

    m_f1 = mm_f_t("{1}");
    m_f1 *= m_f2;
    CHECK(m_f1 == (mm_f_t("{1}") * m_f2));

    m_f1 = mm_f_t("{1}");
    m_f1 /= m_f2;
    CHECK(m_f1 == (mm_f_t("{1}") / m_f2));

    m_f1 = mm_f_t("{1}");
    m_f1 ^= m_f2;
    CHECK(m_f1 == (mm_f_t("{1}") ^ m_f2));

    m_f1 = mm_f_t("{1}");
    m_f1 &= m_f2;
    CHECK(m_f1 == (mm_f_t("{1}") & m_f2));

    m_f1 = mm_f_t("{1}");
    m_f1 %= m_f2;
    CHECK(m_f1 == (mm_f_t("{1}") % m_f2));

    m_f1 = mm_f_t("{1}");
    m_f1 |= m_f2;
    CHECK(m_f1 == (mm_f_t("{1}") | m_f2));

    // Interchangeability: A += B is interchangeable with A = A + B (matching scalar)
    mm_f_t A("{1}");
    mm_f_t B("{2}");
    mm_f_t C = A;
    C += B;
    A = A + B;
    CHECK(A == C);

    // Verify that mixed-precision assignment still works if we use the constructor explicitly
    m_f1 = mm_f_t(m_d1);
    CHECK(m_f1 == mm_f_t("{1}"));
  }

  SUBCASE("Mixed-representation assignment") {
    using fm_d_t = glucat::framed_multi<double, -8, 8>;
    mm_t m("{1,2}");
    fm_d_t f("{1,2}");
    
    // Framed to Matrix (matching scalar)
    mm_t m2;
    m2 = f;
    CHECK(m2 == m);

    // Matrix to Framed
    fm_d_t f2;
    f2 = m2;
    CHECK(f2 == f);
  }

  SUBCASE("Advanced Constructor Gaps") {
    using index_set_t = mm_t::index_set_t;
    using scalar_t = mm_t::scalar_t;
    index_set_t ist(1);
    index_set_t frm(1);
    scalar_t crd(1.0);
    
    // Test matrix_multi constructor with prechecked
    mm_t m1(ist, crd, frm, true);
    CHECK(m1.frame() == frm);
    
    mm_t m2(ist, crd, frm, false);
    CHECK(m2.frame() == frm);
    
    // Trigger exception path for constructor
    index_set_t invalid_ist(2);
    CHECK_THROWS_AS(mm_t(invalid_ist, crd, frm, false), glucat_error);
  }

  SUBCASE("Basis element and specialized members") {
    using index_set_t = mm_t::index_set_t;
    index_set_t frm(1);
    mm_t a("{1}");
    mm_t b_elem(frm);
    CHECK(b_elem == a);
  }

  SUBCASE("Exceptions") {
    mm_t m1(1.0);
    // Negative exponent in outer_pow
    CHECK_THROWS(m1.outer_pow(-1));

    // Construction with value outside of frame
    // mm_t is <-8, 8>, index_set_t(9) is out of frame
    using is_t = mm_t::index_set_t;
    CHECK_THROWS(mm_t(is_t(9), 1.0, is_t(), false));
  }

  SUBCASE("Move Semantics and Reframing") {
    mm_t a(T(1.0), mm_t::index_set_t());
    mm_t b(std::move(a));
    CHECK(b == mm_t(T(1.0)));
    
    mm_t c;
    c = std::move(b);
    CHECK(c == mm_t(T(1.0)));
    
    c = c; // Self-assignment
    CHECK(c == mm_t(T(1.0)));

    // Reframing overlap
    using is_t = mm_t::index_set_t;
    mm_t m1(T(1.0), is_t());
    mm_t m2(is_t("{1}"), T(1.0));
    mm_t sum = m1 + m2; // Triggers reframe
    CHECK(sum.frame() == is_t("{1}"));
    CHECK(sum.scalar() == doctest::Approx(T(1.0)));
    CHECK(sum[is_t("{1}")] == doctest::Approx(T(1.0)));
  }

  SUBCASE("Advanced Matrix Algorithms") {
    mm_t m1(T(1.0), mm_t::index_set_t());
    mm_t m2 = m1.outer_pow(2);
    
    // Test refined_newton_schulz vs newton_schulz
    mm_t a("1+0.1{1}");
    mm_t s1 = glucat::newton_schulz(a, T(1.0));
    mm_t s2 = glucat::refined_newton_schulz(a, T(1.0));
    
    // They should be approximately equal for small perturbations
    CHECK(s1.scalar() == doctest::Approx(s2.scalar()));
  }

  SUBCASE("Numerical Analysis") {
    mm_t f_small("1e8+{1}+1e-8{1,2}");
    CHECK(f_small.truncated(T(1.0e-6)) == mm_t(T(1.0e8), mm_t::index_set_t()));
    
    mm_t f_nan(std::numeric_limits<T>::quiet_NaN());
    mm_t f_inf(std::numeric_limits<T>::infinity());
    
    CHECK(f_nan.isnan());
    CHECK(f_inf.isinf());
  }

  SUBCASE("File I/O") {
    mm_t f("1+{1}");
    f.write("Test prefix"); 
    
    auto temp_path = std::filesystem::temp_directory_path() / ("test_io_mm_" + std::to_string(std::chrono::system_clock::now().time_since_epoch().count()) + ".txt");
    struct Cleanup { std::filesystem::path p; ~Cleanup() { if(std::filesystem::exists(p)) std::filesystem::remove(p); } } cleanup{temp_path};

    std::ofstream ofs(temp_path);
    f.write(ofs, "File prefix");
    ofs.close();
    CHECK(std::filesystem::exists(temp_path));
  }

  SUBCASE("Clifford Algebra Operations") {
    mm_t m1("{1}");
    mm_t m2("{2}");
    mm_t m12 = m1 * m2;
    mm_t m_mix("1+{1}+{1,2}");

    // Involute, Reverse, Conjugate
    CHECK(m1.involute() == -m1);
    CHECK(m12.reverse() == -m12);
    CHECK(m_mix.conj() == mm_t("1-{1}-{1,2}"));

    // Norm, Quad, MaxAbs
    CHECK(m_mix.quad() == doctest::Approx(3.0));
    CHECK(m_mix.norm() == doctest::Approx(3.0));
    CHECK(m_mix.max_abs() == doctest::Approx(1.0));

    // Projections
    CHECK(m_mix.pure() == mm_t("{1}+{1,2}"));
    CHECK(m_mix.even() == mm_t("1+{1,2}"));
    CHECK(m_mix.odd() == mm_t("{1}"));

    // Vector part
    mm_t m_vec("1+2{1}+3{2}+4{1,2}");
    auto v = m_vec.vector_part();
    CHECK(v.size() == 2);
    CHECK(v[0] == doctest::Approx(2.0));
    CHECK(v[1] == doctest::Approx(3.0));

    using is_t = mm_t::index_set_t;
    auto v2 = m_vec.vector_part(is_t("{-1,1,2}"));
    CHECK(v2.size() == 3);
    
    // Grade and Frame
    CHECK(m_mix.grade() == 2);
    CHECK(m_mix.frame() == is_t("{1,2}"));
    
    // Subscripting
    CHECK(m_mix[is_t("{1,2}")] == doctest::Approx(1.0));
    
    // Pow and Inv
    mm_t m_inv = m1.inv();
    CHECK(m_inv == m1);
    CHECK(m1.pow(2) == mm_t(1.0));
  }
}
#endif

#endif  // _GLUCAT_MATRIX_MULTI_IMP_H
