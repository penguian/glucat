#ifndef _GLUCAT_MATRIX_EIGEN_H
#define _GLUCAT_MATRIX_EIGEN_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix.h : Declare common matrix functions
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>
#pragma GCC diagnostic pop

#include <type_traits>
#include <complex>
#include <vector>
#include <iostream>

namespace glucat { namespace matrix
{
  // =========================================================================
  // eigen_matrix_wrapper
  // =========================================================================

  /// Wrapper for Eigen matrix (forward)
  template< typename Scalar_T > class eigen_matrix_wrapper;

  /// Wrapper for Eigen sparse matrix (forward)
  template< typename Scalar_T > class eigen_sparse_wrapper;

  // Output to stream (forward)
  template< typename Scalar_T >
  auto operator<< (std::ostream& os, const eigen_matrix_wrapper<Scalar_T>& m) -> std::ostream&;

  /// Wrapper for Eigen matrix
  template< typename Scalar_T >
  class eigen_matrix_wrapper :
  public matrix_base<eigen_matrix_wrapper<Scalar_T>>
  {
  public:
    using MatrixType = Eigen::Matrix<Scalar_T, Eigen::Dynamic, Eigen::Dynamic>;

    using value_type = Scalar_T;
    using size_type = typename MatrixType::Index;

    MatrixType m_mat;

    /// Constructors
    /// Default constructor
    eigen_matrix_wrapper() = default;

    // Armadillo constructor (rows, cols)
    eigen_matrix_wrapper(matrix_index_t rows, matrix_index_t cols);

    // Constructor from Eigen expressions (e.g. m * s)
    template< typename Derived_T >
    eigen_matrix_wrapper(const Eigen::MatrixBase<Derived_T>& other);

    /// Generic Interop Constructor (e.g. from Armadillo matrix)
    template< typename Other_Matrix_T >
    explicit eigen_matrix_wrapper(const Other_Matrix_T& other);

    // Constructor from eigen_sparse_wrapper
    template< typename Other_Scalar_T >
    explicit eigen_matrix_wrapper(const eigen_sparse_wrapper<Other_Scalar_T>& other);

    // Copy constructor
    // Copy constructor
    eigen_matrix_wrapper(const eigen_matrix_wrapper& other);

    // Move constructor
    // Move constructor
    eigen_matrix_wrapper(eigen_matrix_wrapper&& other) noexcept;

    // Assignment
    // Copy assignment
    auto operator= (const eigen_matrix_wrapper& other) -> eigen_matrix_wrapper&;
    // Move assignment
    auto operator= (eigen_matrix_wrapper&& other) noexcept -> eigen_matrix_wrapper&;

    // Generic Interop Assignment
    template< typename Other_Matrix_T >
    auto operator= (const Other_Matrix_T& other) -> eigen_matrix_wrapper&;

    // Constructor from Eigen
    eigen_matrix_wrapper(const MatrixType& m);
    // Constructor from Eigen (move)
    eigen_matrix_wrapper(MatrixType&& m);

    // Set size
    void set_size(matrix_index_t rows, matrix_index_t cols);

    // Resize
    void resize(matrix_index_t rows, matrix_index_t cols, bool preserve = false);

    // Helpers
    // Number of rows
    auto nbr_rows() const -> matrix_index_t;
    // Number of columns
    auto nbr_cols() const -> matrix_index_t;

    // Clear
    void clear();

    // Set to zero
    void zeros();
    // Set size then set to zero
    void zeros(matrix_index_t rows, matrix_index_t cols);

    // Set to identity
    void unit(matrix_index_t rows, matrix_index_t cols);

    // Is finite?
    auto is_finite() const -> bool;
    // Has NaN?
    auto has_nan() const -> bool;

    // Element access
    // Element access
    auto operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&;
    // Const element access
    auto operator() (matrix_index_t i, matrix_index_t j) const -> const Scalar_T&;

    // Operators
    // Add and assign
    auto operator+= (const eigen_matrix_wrapper& other) -> eigen_matrix_wrapper&;
    // Subtract and assign
    auto operator-= (const eigen_matrix_wrapper& other) -> eigen_matrix_wrapper&;
    // Multiply by scalar and assign
    auto operator*= (const Scalar_T& val) -> eigen_matrix_wrapper&;
    // Divide by scalar and assign
    auto operator/= (const Scalar_T& val) -> eigen_matrix_wrapper&;

    // Addition
    auto operator+ (const eigen_matrix_wrapper& other) const -> eigen_matrix_wrapper;
    // Subtraction
    auto operator- (const eigen_matrix_wrapper& other) const -> eigen_matrix_wrapper;

    // Matrix Multiplication
    auto operator* (const eigen_matrix_wrapper& other) const -> eigen_matrix_wrapper;

    // Unary -
    auto operator- () const -> eigen_matrix_wrapper;

    // Transpose
    auto t() const -> eigen_matrix_wrapper;

    // New Member Functions (formerly free functions)
    // Trace
    auto trace() const;
    // Eigenvalues
    auto eigenvalues() const -> std::vector<std::complex<double>>;
    // Infinity norm
    auto norm_inf() const;
    // Squared Frobenius norm
    auto norm_frob2() const;
    // Is NaN?
    auto isnan() const -> bool;
    // Is infinite?
    auto isinf() const -> bool;
    // Number of non-zeros
    auto nnz() const;

    // Inner product
    template< typename Result_Scalar_T, typename Other >
    auto inner(const Other& other) const -> Result_Scalar_T;

    // Kronecker matrix product
    auto kron(const eigen_matrix_wrapper& other) const -> eigen_matrix_wrapper;
    // Mixed Kronecker matrix product: Dense x Sparse -> Dense (wrapper)
    template< typename Other_Scalar_T >
    auto kron(const eigen_sparse_wrapper<Other_Scalar_T>& other) const -> eigen_matrix_wrapper<Other_Scalar_T>;

    // Left Kronecker quotient
    template< typename RHS_T >
    auto nork(const RHS_T& rhs, bool mono = true) const -> RHS_T;

    // Output to stream
    friend auto operator<< <>(std::ostream& os, const eigen_matrix_wrapper& m) -> std::ostream&;
  };

  // Mixed operations
  // Product of scalar and matrix wrapper
  template< typename Scalar_T >
  auto operator* (Scalar_T s, const eigen_matrix_wrapper<Scalar_T>& m) -> eigen_matrix_wrapper<Scalar_T>
  {
    return eigen_matrix_wrapper<Scalar_T>(s * m.m_mat);
  }
  // Product of matrix wrapper and scalar
  template< typename Scalar_T >
  auto operator* (const eigen_matrix_wrapper<Scalar_T>& m, Scalar_T s) -> eigen_matrix_wrapper<Scalar_T>
  {
    return eigen_matrix_wrapper<Scalar_T>(m.m_mat * s);
  }

  // =========================================================================
  // eigen_sparse_wrapper
  // =========================================================================

  // Output to stream (forward)
  template< typename Scalar_T >
  auto operator<< (std::ostream& os, const eigen_sparse_wrapper<Scalar_T>& m) -> std::ostream&;

  /// Wrapper for Eigen sparse matrix
  template< typename Scalar_T >
  class eigen_sparse_wrapper :
  public matrix_base<eigen_sparse_wrapper<Scalar_T>>
  {
  public:
    using MatrixType = Eigen::SparseMatrix<Scalar_T>;

    using value_type = Scalar_T;
    using size_type = typename MatrixType::Index;

    MatrixType m_mat;

    eigen_sparse_wrapper() = default;

    // Constructor from Eigen Sparse Matrix (e.g. expression result)
    explicit eigen_sparse_wrapper(const MatrixType& m);

    // Armadillo/uBLAS/Generator style constructor support
    eigen_sparse_wrapper(matrix_index_t rows, matrix_index_t cols, matrix_index_t estimated_nnz = 0);

    // Copy/Move similar to dense
    eigen_sparse_wrapper(const eigen_sparse_wrapper& other);

    // Move constructor
    eigen_sparse_wrapper(eigen_sparse_wrapper&& other) noexcept;

    // Copy assignment
    auto operator= (const eigen_sparse_wrapper& other) -> eigen_sparse_wrapper&;

    // Move assignment
    auto operator= (eigen_sparse_wrapper&& other) noexcept -> eigen_sparse_wrapper&;

    // Set size
    void set_size(matrix_index_t rows, matrix_index_t cols);

    // Make writable
    void resize(matrix_index_t rows, matrix_index_t cols, bool preserve = false);

    // Clear
    void clear();

    // Set to zero
    void zeros();

    // Set to zero
    void zeros(matrix_index_t rows, matrix_index_t cols);

    // Iterator support
    class const_iterator
    {
    public:
      using InnerIterator = typename MatrixType::InnerIterator;

      const MatrixType* mp_mat;
      int m_outer;
      InnerIterator m_inner;

      // Constructor for begin()
      const_iterator(const MatrixType* mat, bool start = true);

      // Advance iterator
      void advance();

      // Check if end
      auto is_end() const -> bool;
      // Prefix increment
      auto operator++ () -> const_iterator&;

      // Inequality comparison
      auto operator!= (const const_iterator& other) const -> bool;

      // Row index
      auto row() const -> matrix_index_t;
      // Column index
      auto col() const -> matrix_index_t;
      // Dereference
      auto operator* () const -> Scalar_T;
    };

    // Iterator support
    // Begin iterator
    auto begin() const -> const_iterator;
    // End iterator
    auto end() const -> const_iterator;

    // Number of rows
    auto nbr_rows() const -> matrix_index_t;
    // Number of columns
    auto nbr_cols() const -> matrix_index_t;

    // Const element access
    auto operator() (matrix_index_t i, matrix_index_t j) const -> Scalar_T;
    // Element access
    auto operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&;

    // Add and assign
    auto operator+= (const eigen_sparse_wrapper& other) -> eigen_sparse_wrapper&;
    // Subtract and assign
    auto operator-= (const eigen_sparse_wrapper& other) -> eigen_sparse_wrapper&;

    // Multiply by sparse wrapper
    auto operator* (const eigen_sparse_wrapper& other) const -> eigen_sparse_wrapper;

    // Multiply by scalar and assign
    auto operator*= (const Scalar_T& val) -> eigen_sparse_wrapper&;

    // New Member Functions
    // Trace
    auto trace() const; // Trace of sparse?
    // Eigenvalues
    auto eigenvalues() const -> std::vector<std::complex<double>>;
    // Infinity norm
    auto norm_inf() const;
    // Squared Frobenius norm
    auto norm_frob2() const;
    // Is NaN?
    auto isnan() const -> bool;
    // Is infinite?
    auto isinf() const -> bool;
    // Number of non-zeros
    auto nnz() const;

    // Inner product
    template< typename Result_Scalar_T, typename Other >
    auto inner(const Other& other) const -> Result_Scalar_T;

    // Mixed Kronecker matrix product: Sparse x Dense -> Dense (wrapper)
    template< typename Other_Scalar_T >
    auto kron(const eigen_matrix_wrapper<Other_Scalar_T>& other) const -> eigen_matrix_wrapper<Other_Scalar_T>;
    // Kronecker matrix product of sparse wrappers
    auto kron(const eigen_sparse_wrapper& other) const -> eigen_sparse_wrapper;

    // Left Kronecker quotient
    template< typename RHS_T >
    auto nork(const RHS_T& rhs, bool mono = true) const -> RHS_T;

    // Output to stream
    friend auto operator<< <>(std::ostream& os, const eigen_sparse_wrapper& m) -> std::ostream&;
  };

  // Product of sparse wrapper and scalar
  template< typename Scalar_T >
  auto operator* (const eigen_sparse_wrapper<Scalar_T>& m, Scalar_T s) -> eigen_sparse_wrapper<Scalar_T>
  {
    eigen_sparse_wrapper<Scalar_T> res(m);
    res *= s;
    return res;
  }

  // Product of scalar and sparse wrapper
  template< typename Scalar_T >
  auto operator* (Scalar_T s, const eigen_sparse_wrapper<Scalar_T>& m) -> eigen_sparse_wrapper<Scalar_T>
  {
    return m * s;
  }

  // Sum of sparse wrappers
  template< typename Scalar_T >
  auto operator+ (const eigen_sparse_wrapper<Scalar_T>& lhs, const eigen_sparse_wrapper<Scalar_T>& rhs) -> eigen_sparse_wrapper<Scalar_T>
  {
    eigen_sparse_wrapper<Scalar_T> res(lhs);
    res += rhs;
    return res;
  }

  // Difference of sparse wrappers
  template< typename Scalar_T >
  auto operator- (const eigen_sparse_wrapper<Scalar_T>& lhs, const eigen_sparse_wrapper<Scalar_T>& rhs) -> eigen_sparse_wrapper<Scalar_T>
  {
    eigen_sparse_wrapper<Scalar_T> res(lhs);
    res -= rhs;
    return res;
  }
} }

#endif  // _GLUCAT_MATRIX_EIGEN_H
