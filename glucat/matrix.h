#ifndef _GLUCAT_MATRIX_H
#define _GLUCAT_MATRIX_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix.h : Declare common matrix functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2012 by Paul C. Leopardi
                         : uBLAS interface contributed by Joerg Walter
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

#include "glucat/glucat_config.h"

// Check if Armadillo should be used
#if defined(_GLUCAT_USE_ARMADILLO)

  #include <armadillo>
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#pragma GCC diagnostic pop
#include <type_traits>
#include <complex>
#include <vector>

#if defined(_GLUCAT_USE_QD)
  #include <qd/dd_real.h>
  #include <qd/qd_real.h>
#endif

namespace glucat {

  // Forward declaration of the wrappers
  template<typename Scalar_T> class eigen_matrix_wrapper;
  template<typename Scalar_T> class arma_matrix_wrapper;
  template<typename Scalar_T> class eigen_sparse_wrapper;

  // Trait to determine if T is natively supported by Armadillo
  template<typename T>
  struct is_arma_supported : std::false_type {};

  #if defined(_GLUCAT_USE_ARMADILLO)
    template<> struct is_arma_supported<float> : std::true_type {};
    template<> struct is_arma_supported<double> : std::true_type {};
    template<> struct is_arma_supported<std::complex<float>> : std::true_type {};
    template<> struct is_arma_supported<std::complex<double>> : std::true_type {};
    
    // Explicitly NOTE: long double, dd_real, qd_real will use the default false_type
    // and thus be dispatched to the Eigen wrapper.
  #endif

  // Dense Selector
  template<typename Scalar_T, bool UseArma = is_arma_supported<Scalar_T>::value>
  struct matrix_type_selector {
    using type = eigen_matrix_wrapper<Scalar_T>;
  };

  #if defined(_GLUCAT_USE_ARMADILLO)
    template<typename Scalar_T>
    struct matrix_type_selector<Scalar_T, true> {
      using type = arma_matrix_wrapper<Scalar_T>;
    };
  #endif

  template<typename Scalar_T>
  using matrix_t = typename matrix_type_selector<Scalar_T>::type;

  // Sparse Selector
  template<typename Scalar_T, bool UseArma = is_arma_supported<Scalar_T>::value>
  struct sparse_matrix_type_selector {
    using type = eigen_sparse_wrapper<Scalar_T>;
  };

  #if defined(_GLUCAT_USE_ARMADILLO)
    template<typename Scalar_T>
    struct sparse_matrix_type_selector<Scalar_T, true> {
      using type = arma::SpMat<Scalar_T>;
    };
  #endif

  template<typename Scalar_T>
  using sparse_matrix_t = typename sparse_matrix_type_selector<Scalar_T>::type;

  namespace matrix
  {
    // ... existing function declarations can remain or be updated if signatures change ...
    // For now, we are focusing on the backend type definition.
    // The previous implementation used uBLAS matrix expressions. 
    // The new design moves towards an Armadillo-like interface.
    // We should keep the namespaces and general structure.
    
    // NOTE: The previous matrix.h contained declarations for functions like kron, etc.
    // We will keep them but they might need template parameter adjustments if they assume uBLAS.
    // However, the task specifically asked for "matrix.h" interface updates.
    // The bulk of the logic is likely moving to free functions compatible with the new matrix_t.
    
    // We will re-declare the essential functions consistent with the new types.
    
    /// Kronecker tensor product of matrices - as per Matlab kron
    template< typename LHS_T, typename RHS_T >
    auto
    kron(const LHS_T& lhs, const RHS_T& rhs) -> const
    RHS_T;

    /// Sparse Kronecker tensor product of monomial matrices
    template< typename LHS_T, typename RHS_T >
    auto
    mono_kron(const LHS_T& lhs, const RHS_T& rhs) -> const
    RHS_T;

    /// Left inverse of Kronecker product
    template< typename LHS_T, typename RHS_T >
    auto
    nork(const LHS_T& lhs, const RHS_T& rhs, const bool mono = true) -> const
    RHS_T;

    /// Left inverse of Kronecker product where lhs is a signed permutation matrix
    template< typename LHS_T, typename RHS_T >
    auto
    signed_perm_nork(const LHS_T& lhs, const RHS_T& rhs) -> const
    RHS_T;

    /// Number of non-zeros
    template< typename Matrix_T >
    auto
    nnz(const Matrix_T& m) -> typename Matrix_T::size_type; // size_type check?

    /// Infinite
    template< typename Matrix_T >
    auto
    isinf(const Matrix_T& m) -> bool;

    /// Not a Number
    template< typename Matrix_T >
    auto
    isnan(const Matrix_T& m) -> bool;

    /// Unit matrix - as per Matlab eye
    template< typename Matrix_T >
    auto
    unit(const size_t dim) -> const
    Matrix_T; // Simplified signature from uBLAS

    // Monomial and Prod functions - simplified/generic
    template< typename LHS_T, typename RHS_T >
    auto
    mono_prod(const LHS_T& lhs,
              const RHS_T& rhs) -> const
    decltype(lhs * rhs); // Assume operator* works

    template< typename LHS_T, typename RHS_T >
    auto
    sparse_prod(const LHS_T& lhs,
                const RHS_T& rhs) -> const
    decltype(lhs * rhs);

    template< typename LHS_T, typename RHS_T >
    auto
    prod(const LHS_T& lhs,
         const RHS_T& rhs) -> const
    decltype(lhs * rhs);

    /// Inner product: sum(x(i,j)*y(i,j))/x.nrows()
    template< typename Scalar_T, typename LHS_T, typename RHS_T >
    auto
    inner(const LHS_T& lhs, const RHS_T& rhs) -> Scalar_T;

    /// Square of Frobenius norm
    template< typename Matrix_T >
    auto
    norm_frob2(const Matrix_T& val) -> typename Matrix_T::elem_type; // Use elem_type for Arma/Wrapper compatibility

    /// Matrix trace
    template< typename Matrix_T >
    auto
    trace(const Matrix_T& val) -> typename Matrix_T::elem_type;

    /// Eigenvalues of a matrix
    template< typename Matrix_T >
    auto
    eigenvalues(const Matrix_T& val) -> std::vector< std::complex<double> >;

    /// Classification of eigenvalues of a matrix
    using eig_case_t = enum {
      safe_eigs,
      neg_real_eigs,
      both_eigs};

    ///  Structure containing classification of eigenvalues
    template< typename Matrix_T >
    struct eig_genus
    {
      using Scalar_T = typename Matrix_T::elem_type; // elem_type
      /// Is the matrix singular?
      bool m_is_singular = false;
      /// What kind of eigenvalues does the matrix contain?
      eig_case_t m_eig_case = safe_eigs;
      /// Argument such that exp(pi-m_safe_arg) lies between arguments of eigenvalues
      Scalar_T   m_safe_arg = Scalar_T(0);
    };

    /// Classify the eigenvalues of a matrix
    template< typename Matrix_T >
    auto
    classify_eigenvalues(const Matrix_T& val) -> eig_genus<Matrix_T>;
  }
}

#include "glucat/matrix_imp.h" // For wrapper implementation and function definitions

#endif  // _GLUCAT_MATRIX_H
