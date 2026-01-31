#ifndef _GLUCAT_MATRIX_H
#define _GLUCAT_MATRIX_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix.h : Declare common matrix classes and functions
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
 ***************************************************************************/

#include "glucat/matrix_base.h"

#include "glucat/matrix_eigen.h"
#if defined(_GLUCAT_USE_ARMADILLO)
#include "glucat/matrix_arma.h"
#endif

#include <type_traits>

namespace glucat { namespace matrix
{
  // Traits Specializations

  /// Trait to determine if Scalar_T is natively supported by Armadillo
  template< typename Scalar_T >
  struct is_arma_supported : std::false_type {};

  // Declarations to satisfy the compiler
  template< typename Scalar_T > class eigen_matrix_wrapper;
  template< typename Scalar_T > class eigen_sparse_wrapper;

  // Dense Selector
  /// Dense matrix type selector
  template< typename Scalar_T, bool UseArma = is_arma_supported<Scalar_T>::value >
  struct matrix_type_selector
  {
    using type = eigen_matrix_wrapper<Scalar_T>;
  };

  /// Dense matrix type selector
  template< typename Scalar_T >
  using matrix_t = typename matrix_type_selector<Scalar_T>::type;

  /// Sparse matrix type selector
  template< typename Scalar_T, bool UseArma = is_arma_supported<Scalar_T>::value >
  struct sparse_matrix_type_selector
  {
    using type = eigen_sparse_wrapper<Scalar_T>;
  };

  /// Sparse matrix type selector
  template< typename Scalar_T >
  using sparse_matrix_t = typename sparse_matrix_type_selector<Scalar_T>::type;

  // ===========================================================
  // Matrix Template Classes (Facade)
  // Named dense_matrix to avoid collision with namespace matrix
  // ===========================================================

  /// Dense matrix class
  template< typename Scalar_T >
  class dense_matrix :
  public matrix_type_selector<Scalar_T>::type
  {
  public:
    using Base = typename matrix_type_selector<Scalar_T>::type;
    using Base::Base; // Inherit constructors
    using Base::operator=;

    dense_matrix() = default;
    dense_matrix(const dense_matrix&) = default;
    dense_matrix(dense_matrix&&) = default;
    auto operator= (const dense_matrix&) -> dense_matrix& = default;
    auto operator= (dense_matrix&&) -> dense_matrix& = default;

    template< typename Other_Matrix_T >
    dense_matrix(const Other_Matrix_T& other) : Base(other) {}
  };

  /// Sparse matrix class
  template< typename Scalar_T >
  class sparse_matrix :
  public sparse_matrix_type_selector<Scalar_T>::type
  {
  public:
    using Base = typename sparse_matrix_type_selector<Scalar_T>::type;
    using Base::Base; // Inherit constructors
    using Base::operator=;

    sparse_matrix() = default;
    sparse_matrix(const sparse_matrix&) = default;
    sparse_matrix(sparse_matrix&&) = default;
    auto operator= (const sparse_matrix&) -> sparse_matrix& = default;
    auto operator= (sparse_matrix&&) -> sparse_matrix& = default;

    template< typename Other_Matrix_T >
    sparse_matrix(const Other_Matrix_T& other) : Base(other) {}
  };

} }

#endif  // _GLUCAT_MATRIX_H
