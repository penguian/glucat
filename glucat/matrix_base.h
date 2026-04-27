#ifndef _GLUCAT_MATRIX_BASE_H
#define _GLUCAT_MATRIX_BASE_H
/***************************************************************************
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

#include <type_traits>
#include <complex>
#include <cstddef>

namespace glucat { namespace matrix
{
  // Default index to use with matrices
  using matrix_index_t = std::size_t;

  /// Helper for static_assert in templates
  template< typename... > struct dependent_false : std::false_type {};

  /// =========================================================================
  /// Traits
  /// =========================================================================
  /// Helper trait for complex check
  template< typename Scalar_T > struct is_complex_t : std::false_type {};
  template< typename Scalar_T > struct is_complex_t<std::complex<Scalar_T>> : std::true_type {};

  /// Classification of eigenvalues of a matrix
  enum eig_case_t
  {
    safe_eigs,
    neg_real_eigs,
    both_eigs
  };

  ///  Structure containing classification of eigenvalues
  template< typename Matrix_T >
  struct eig_genus
  {
    using Scalar_T = typename Matrix_T::value_type;
    // Is the matrix singular?
    bool m_is_singular = false;
    // What kind of eigenvalues does the matrix contain?
    eig_case_t m_eig_case = safe_eigs;
    // Argument such that exp(pi-m_safe_arg) lies between arguments of eigenvalues
    Scalar_T   m_safe_arg = Scalar_T(0);
  };

  /// Base class providing member functions that delegate to the derived implementation (CRTP Pattern)
  template< typename Derived_T >
  class matrix_base
  {
  public:
    /// Element access
    inline auto operator() (matrix_index_t i, matrix_index_t j) -> decltype(auto)
    { return derived()(i, j); }

    /// Const element access
    inline auto operator() (matrix_index_t i, matrix_index_t j) const -> decltype(auto)
    { return derived()(i, j); }

    // Return const reference to derived class
    auto derived() const -> const Derived_T&;
    // Return reference to derived class
    auto derived() -> Derived_T&;

    // Member functions delegating to namespace matrix implementation
    // defined in matrix_base_imp.h

    // Generic classify_eigenvalues relies on eigenvalues() member
    auto classify_eigenvalues() const -> eig_genus<Derived_T>;
  };

  // Core Operations as Free Functions

  // Number of rows
  template< typename Matrix_T >
  auto nbr_rows(const Matrix_T& mat) -> matrix_index_t;

  // Number of columns
  template< typename Matrix_T >
  auto nbr_cols(const Matrix_T& mat) -> matrix_index_t;

  /// Helper struct for unit matrix creation
  template< typename Matrix_T > struct unit_helper;

  // Identity matrix
  template< typename Matrix_T >
  auto unit(const matrix_index_t dim) -> Matrix_T;

} }

#endif  // _GLUCAT_MATRIX_BASE_H
