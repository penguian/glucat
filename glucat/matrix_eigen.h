#ifndef _GLUCAT_MATRIX_EIGEN_H
#define _GLUCAT_MATRIX_EIGEN_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_eigen.h : Declare Eigen matrix wrappers
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

#if defined(_GLUCAT_USE_QD)
# include "glucat/qd.h"
#endif

#ifdef _GLUCAT_USE_GCC_PRAGMAS
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
# if defined(__clang__)
#  pragma GCC diagnostic ignored "-Wunused-but-set-variable"
# endif
#endif

#if defined(_GLUCAT_USE_QD)
#include <Eigen/Core>
#endif

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>
#ifdef _GLUCAT_USE_GCC_PRAGMAS
#pragma GCC diagnostic pop
#endif

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
  std::ostream& operator<< (std::ostream& os, const eigen_matrix_wrapper<Scalar_T>& m);

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

    // Sized constructor (rows, cols)
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
    eigen_matrix_wrapper(const eigen_matrix_wrapper& other);

    // Move constructor
    eigen_matrix_wrapper(eigen_matrix_wrapper&& other) noexcept;

    // Assignment
    // Copy assignment
    eigen_matrix_wrapper& operator= (const eigen_matrix_wrapper& other);
    // Move assignment
    eigen_matrix_wrapper& operator= (eigen_matrix_wrapper&& other) noexcept;

    // Generic Interop Assignment
    template< typename Other_Matrix_T >
    eigen_matrix_wrapper& operator= (const Other_Matrix_T& other);

    // Constructor from Eigen
    eigen_matrix_wrapper(const MatrixType& m);
    // Constructor from Eigen (move)
    eigen_matrix_wrapper(MatrixType&& m);

    // Set size
    void set_size(matrix_index_t rows, matrix_index_t cols);



    // Helpers
    // Number of rows
    matrix_index_t nbr_rows() const;
    // Number of columns
    matrix_index_t nbr_cols() const;

    // Clear
    void clear();

    // Set to zero
    void zeros();
    // Set size then set to zero
    void zeros(matrix_index_t rows, matrix_index_t cols);

    // Set to identity
    void unit(matrix_index_t rows, matrix_index_t cols);

    // Is finite?
    bool is_finite() const;
    // Has NaN?
    bool has_nan() const;

    // Element access
    Scalar_T& operator() (matrix_index_t i, matrix_index_t j);
    // Const element access
    const Scalar_T& operator() (matrix_index_t i, matrix_index_t j) const;

    // Operators
    // Add and assign
    eigen_matrix_wrapper& operator+= (const eigen_matrix_wrapper& other);
    // Subtract and assign
    eigen_matrix_wrapper& operator-= (const eigen_matrix_wrapper& other);
    // Multiply by scalar and assign
    eigen_matrix_wrapper& operator*= (const Scalar_T& val);
    // Divide by scalar and assign
    eigen_matrix_wrapper& operator/= (const Scalar_T& val);

    // Addition
    eigen_matrix_wrapper operator+ (const eigen_matrix_wrapper& other) const;
    // Subtraction
    eigen_matrix_wrapper operator- (const eigen_matrix_wrapper& other) const;

    // Matrix Multiplication
    eigen_matrix_wrapper operator* (const eigen_matrix_wrapper& other) const;

    // Unary -
    eigen_matrix_wrapper operator- () const;

    // Transpose
    eigen_matrix_wrapper t() const;

    // New Member Functions (formerly free functions)
    // Trace
    Scalar_T trace() const;
    // Eigenvalues
    std::vector<std::complex<double>> eigenvalues() const;
    // Infinity norm
    typename Eigen::NumTraits<Scalar_T>::Real norm_inf() const;
    // Squared Frobenius norm
    typename Eigen::NumTraits<Scalar_T>::Real norm_frob2() const;
    // Is NaN?
    bool isnan() const;
    // Is infinite?
    bool isinf() const;
    // Number of non-zeros
    matrix_index_t nnz() const;

    // Inner product
    template< typename Result_Scalar_T, typename Other >
    Result_Scalar_T inner(const Other& other) const;

    // Kronecker matrix product
    eigen_matrix_wrapper kron(const eigen_matrix_wrapper& other) const;
    // Mixed Kronecker matrix product: Dense x Sparse -> Dense (wrapper)
    template< typename Other_Scalar_T >
    eigen_matrix_wrapper<Other_Scalar_T> kron(const eigen_sparse_wrapper<Other_Scalar_T>& other) const;

    // Left Kronecker quotient
    template< typename RHS_T >
    RHS_T nork(const RHS_T& rhs, bool mono = true) const;

    // Output to stream
    friend std::ostream& operator<< <>(std::ostream& os, const eigen_matrix_wrapper& m);
  };

  // Mixed operations
  // Product of scalar and matrix wrapper
  template< typename Scalar_T >
  eigen_matrix_wrapper<Scalar_T> operator* (Scalar_T s, const eigen_matrix_wrapper<Scalar_T>& m)
  {
    return eigen_matrix_wrapper<Scalar_T>(s * m.m_mat);
  }
  // Product of matrix wrapper and scalar
  template< typename Scalar_T >
  eigen_matrix_wrapper<Scalar_T> operator* (const eigen_matrix_wrapper<Scalar_T>& m, Scalar_T s)
  {
    return eigen_matrix_wrapper<Scalar_T>(m.m_mat * s);
  }

  // =========================================================================
  // eigen_sparse_wrapper
  // =========================================================================

  // Output to stream (forward)
  template< typename Scalar_T >
  std::ostream& operator<< (std::ostream& os, const eigen_sparse_wrapper<Scalar_T>& m);

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
    eigen_sparse_wrapper& operator= (const eigen_sparse_wrapper& other);

    // Move assignment
    eigen_sparse_wrapper& operator= (eigen_sparse_wrapper&& other) noexcept;

    // Set size
    void set_size(matrix_index_t rows, matrix_index_t cols);



    // Clear
    void clear();

    // Set to zero
    void zeros();

    // Set size then set to zero
    void zeros(matrix_index_t rows, matrix_index_t cols);

    // Set to identity
    void unit(matrix_index_t rows, matrix_index_t cols);

    // Iterator support
    class const_iterator
    {
    public:
      using InnerIterator = typename MatrixType::InnerIterator;

      const MatrixType* mp_mat;
      matrix_index_t m_outer;
      InnerIterator m_inner;

      // Constructor for begin()
      const_iterator(const MatrixType* mat, bool start = true);

      // Advance iterator
      void advance();

      // Check if end
      bool is_end() const;
      // Prefix increment
      const_iterator& operator++ ();

      // Inequality comparison
      bool operator!= (const const_iterator& other) const;

      // Row index
      matrix_index_t row() const;
      // Column index
      matrix_index_t col() const;
      // Dereference
      Scalar_T operator* () const;
    };

    // Iterator support
    // Begin iterator
    const_iterator begin() const;
    // End iterator
    const_iterator end() const;

    // Number of rows
    matrix_index_t nbr_rows() const;
    // Number of columns
    matrix_index_t nbr_cols() const;

    // Const element access
    Scalar_T operator() (matrix_index_t i, matrix_index_t j) const;
    // Element access
    Scalar_T& operator() (matrix_index_t i, matrix_index_t j);

    // Add and assign
    eigen_sparse_wrapper& operator+= (const eigen_sparse_wrapper& other);
    // Subtract and assign
    eigen_sparse_wrapper& operator-= (const eigen_sparse_wrapper& other);

    // Multiply by sparse wrapper
    eigen_sparse_wrapper operator* (const eigen_sparse_wrapper& other) const;

    // Multiply by scalar and assign
    eigen_sparse_wrapper& operator*= (const Scalar_T& val);

    // New Member Functions
    // Trace
    Scalar_T trace() const;
    // Eigenvalues
    std::vector<std::complex<double>> eigenvalues() const;
    // Infinity norm
    typename Eigen::NumTraits<Scalar_T>::Real norm_inf() const;
    // Squared Frobenius norm
    typename Eigen::NumTraits<Scalar_T>::Real norm_frob2() const;
    // Is NaN?
    bool isnan() const;
    // Is infinite?
    bool isinf() const;
    // Number of non-zeros
    matrix_index_t nnz() const;

    // Inner product
    template< typename Result_Scalar_T, typename Other >
    Result_Scalar_T inner(const Other& other) const;

    // Mixed Kronecker matrix product: Sparse x Dense -> Dense (wrapper)
    template< typename Other_Scalar_T >
    eigen_matrix_wrapper<Other_Scalar_T> kron(const eigen_matrix_wrapper<Other_Scalar_T>& other) const;
    // Kronecker matrix product of sparse wrappers
    eigen_sparse_wrapper kron(const eigen_sparse_wrapper& other) const;

    // Left Kronecker quotient
    template< typename RHS_T >
    RHS_T nork(const RHS_T& rhs, bool mono = true) const;

    // Output to stream
    friend std::ostream& operator<< <>(std::ostream& os, const eigen_sparse_wrapper& m);
  };

  // Product of sparse wrapper and scalar
  template< typename Scalar_T >
  eigen_sparse_wrapper<Scalar_T> operator* (const eigen_sparse_wrapper<Scalar_T>& m, Scalar_T s);

  // Product of scalar and sparse wrapper
  template< typename Scalar_T >
  eigen_sparse_wrapper<Scalar_T> operator* (Scalar_T s, const eigen_sparse_wrapper<Scalar_T>& m);

  // Sum of sparse wrappers
  template< typename Scalar_T >
  eigen_sparse_wrapper<Scalar_T> operator+ (const eigen_sparse_wrapper<Scalar_T>& lhs, const eigen_sparse_wrapper<Scalar_T>& rhs);

  // Difference of sparse wrappers
  template< typename Scalar_T >
  eigen_sparse_wrapper<Scalar_T> operator- (const eigen_sparse_wrapper<Scalar_T>& lhs, const eigen_sparse_wrapper<Scalar_T>& rhs);
} }

#endif  // _GLUCAT_MATRIX_EIGEN_H
