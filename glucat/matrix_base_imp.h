#ifndef _GLUCAT_MATRIX_BASE_IMP_H
#define _GLUCAT_MATRIX_BASE_IMP_H
/**************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_base.h : Declare common matrix class
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2026 by Paul C. Leopardi
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

#include "glucat/errors.h"
#include "glucat/scalar.h"
#include "glucat/matrix.h"

#include <set>
#include <vector>
#include <complex>
#include <type_traits>

#include <algorithm>
#include <limits>

namespace glucat { namespace matrix
{

  // =========================================================================
  // matrix_base Member Definitions
  // =========================================================================

  /*
   * @brief Return const reference to derived class
   * @details
   * @tparam Derived_T
   * @return Derived class
   */
  template< typename Derived_T >
  inline auto
  matrix_base<Derived_T>::
  derived() const -> const Derived_T&
  { return static_cast<const Derived_T&>(*this); }

  /*
   * @brief Return reference to derived class
   * @details
   * @tparam Derived_T
   * @return Derived class
   */
  template< typename Derived_T >
  inline auto
  matrix_base<Derived_T>::
  derived() -> Derived_T&
  { return static_cast<Derived_T&>(*this); }

  /*
   * @brief Generic classify_eigenvalues relies on eigenvalues() member
   * @details
   * @tparam Derived_T
   * @retval result Eigenvalue classification
   */
  template< typename Derived_T >
  inline auto
  matrix_base<Derived_T>::
  classify_eigenvalues() const
  {
    using Scalar_T = typename Derived_T::value_type;
    eig_genus<Derived_T> result;

    auto lambda = derived().eigenvalues(); // Call member

    std::set<double> arg_set;

    const auto dim = lambda.size();
    static const auto epsilon =
      std::max(std::numeric_limits<double>::epsilon(),
               numeric_traits<Scalar_T>::to_double(std::numeric_limits<Scalar_T>::epsilon()));
    static const auto zero_eig_tol = 4096.0 * epsilon;

    bool neg_real_eig_found = false;
    bool imag_eig_found = false;
    bool zero_eig_found = false;

    for (auto
         k = decltype(dim)(0);
         k != dim;
         ++k)
    {
      const auto lambda_k = lambda[k];
      arg_set.insert(std::arg(lambda_k));

      const auto real_lambda_k = std::real(lambda_k);
      const auto imag_lambda_k = std::imag(lambda_k);
      const auto norm_tol = 4096.0 * epsilon * std::norm(lambda_k);

      if (!neg_real_eig_found &&
          real_lambda_k < -epsilon &&
          (imag_lambda_k == 0.0 ||
           imag_lambda_k * imag_lambda_k < norm_tol))
        neg_real_eig_found = true;
      if (!imag_eig_found &&
          imag_lambda_k > epsilon &&
          (real_lambda_k == 0.0 ||
           real_lambda_k * real_lambda_k < norm_tol))
        imag_eig_found = true;
      if (!zero_eig_found &&
          std::norm(lambda_k) < zero_eig_tol)
        zero_eig_found = true;
    }

    if (zero_eig_found)
      result.m_is_singular = true;

    static const auto pi = numeric_traits<double>::pi();
    if (neg_real_eig_found)
    {
      if (imag_eig_found)
        result.m_eig_case = both_eigs;
      else
      {
        result.m_eig_case = neg_real_eigs;
        result.m_safe_arg = Scalar_T(-pi / 2.0);
      }
    }

    if (result.m_eig_case == both_eigs)
    {
      auto arg_it = arg_set.begin();
      auto first_arg = *arg_it;
      auto best_arg = first_arg;
      auto best_diff = 0.0;
      auto previous_arg = first_arg;
      for (++arg_it;
           arg_it != arg_set.end();
           ++arg_it)
      {
        const auto arg_diff = *arg_it - previous_arg;
        if (arg_diff > best_diff)
        {
          best_diff = arg_diff;
          best_arg = previous_arg;
        }
        previous_arg = *arg_it;
      }
      const auto arg_diff = first_arg + 2.0 * pi - previous_arg;
      if (arg_diff > best_diff)
      {
        best_diff = arg_diff;
        best_arg = previous_arg;
      }
      result.m_safe_arg = Scalar_T(pi - (best_arg + best_diff / 2.0));
    }
    return result;
  }
  // =========================================================================
  // nbr_rows and nbr_cols Definitions
  // =========================================================================

  /*
   * @brief Number of rows
   * @details
   * @tparam Matrix_T
   * @param mat Matrix
   * @return Number of rows
   */
  template< typename Matrix_T >
  inline auto nbr_rows(const Matrix_T& mat) -> matrix_index_t
  {
    if constexpr(requires { mat.nbr_rows(); })
      return mat.nbr_rows();
    else if constexpr(requires { mat.n_rows; })
      return mat.n_rows;
    else if constexpr(requires { mat.rows(); })
      return mat.rows();
    else
      static_assert(dependent_false<Matrix_T>::value, "Unsupported matrix type");
  }

  /*
   * @brief Number of columns
   * @details
   * @tparam Matrix_T
   * @param mat Matrix
   * @return Number of columns
   */
  template< typename Matrix_T >
  inline auto nbr_cols(const Matrix_T& mat) -> matrix_index_t
  {
    if constexpr(requires { mat.nbr_cols(); })
      return mat.nbr_cols();
    else if constexpr(requires { mat.n_cols; })
      return mat.n_cols;
    else if constexpr(requires { mat.cols(); })
      return mat.cols();
    else
      static_assert(dependent_false<Matrix_T>::value, "Unsupported matrix type");
  }

  // =========================================================================
  // unit_helper and unit Definitions
  // =========================================================================

   /*
    * @brief Helper struct for unit matrix creation
    * @details
    * @tparam Matrix_T
    */
  template< typename Matrix_T >
  struct unit_helper
  {
    static inline auto apply(matrix_index_t dim) -> Matrix_T
    {
      Matrix_T result(dim, dim);
      // Fallback logic
      if constexpr (requires { result.unit(dim, dim); })
        result.unit(dim, dim);
      else if constexpr (requires { result.eye(dim, dim); })
        result.eye(dim, dim);
      else if constexpr (requires { result.setIdentity(); })
        result.setIdentity();
      else
      {
         // Manual identity
         if constexpr (requires { result.zeros(); })
           result.zeros();
         else if constexpr (requires { result.clear(); })
           result.clear();

         for (matrix_index_t i = 0; i < dim; ++i)
           result(i, i) = static_cast<typename Matrix_T::value_type>(1);
      }
      return result;
    }
  };

  /*
   * @brief Identity matrix
   * @details
   * @tparam Matrix_T
   * @param dim Value
   * @return Identity matrix
   */
  template< typename Matrix_T >
  inline auto unit(const matrix_index_t dim) -> Matrix_T
  {
    return unit_helper<Matrix_T>::apply(dim);
  }

} }

#endif // _GLUCAT_MATRIX_BASE_IMP_H
