#ifndef _GLUCAT_MATRIX_ARMA_H
#define _GLUCAT_MATRIX_ARMA_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_arma.h : Declare common matrix functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2026 by Paul C. Leopardi
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

#if defined(_GLUCAT_USE_ARMADILLO)

#include "glucat/matrix_base.h"

#include <armadillo>

#include <type_traits>
#include <complex>
#include <vector>
#include <iostream>

namespace glucat { namespace matrix
{
  /// Matrix wrapper for Armadillo (forward)
  template< typename Scalar_T > class arma_matrix_wrapper;
  /// Sparse matrix wrapper for Armadillo (forward)
  template< typename Scalar_T > class arma_sparse_wrapper;

  // Output to stream (forward)
  template< typename Scalar_T >
  std::ostream& operator<< (std::ostream& os, const arma_matrix_wrapper<Scalar_T>& m);

  /// Wrapper for Armadillo matrix
  template< typename Scalar_T >
  class arma_matrix_wrapper :
  public matrix_base<arma_matrix_wrapper<Scalar_T>>
  {
  public:
    using MatrixType = arma::Mat<Scalar_T>;

    using value_type = Scalar_T;
    using size_type = matrix_index_t;

    MatrixType m_mat;

    // Constructors
    /// Default constructor
    arma_matrix_wrapper() = default;

    // Constructor with size
    arma_matrix_wrapper(matrix_index_t rows, matrix_index_t cols);

    /// Constructor from other matrix type
    template< typename Other_Matrix_T >
    explicit arma_matrix_wrapper(const Other_Matrix_T& other);

    // Constructor from eigen_matrix_wrapper
    template< typename Other_Scalar_T >
    explicit arma_matrix_wrapper(const eigen_matrix_wrapper<Other_Scalar_T>& other);

    // Constructor from eigen_sparse_wrapper
    template< typename Other_Scalar_T >
    explicit arma_matrix_wrapper(const eigen_sparse_wrapper<Other_Scalar_T>& other);

    // Constructor from arma_sparse_wrapper
    template< typename Other_Scalar_T >
    explicit arma_matrix_wrapper(const arma_sparse_wrapper<Other_Scalar_T>& other);

    // Copy/Move
    // Copy constructor
    arma_matrix_wrapper(const arma_matrix_wrapper& other);
    // Move constructor
    arma_matrix_wrapper(arma_matrix_wrapper&& other) noexcept(std::is_nothrow_move_constructible_v<Scalar_T>);

    // Copy assignment
    arma_matrix_wrapper& operator= (const arma_matrix_wrapper& other);
    // Move assignment
    arma_matrix_wrapper& operator= (arma_matrix_wrapper&& other) noexcept(std::is_nothrow_move_assignable_v<Scalar_T>);
    // Assignment from sparse wrapper
    arma_matrix_wrapper& operator= (const arma_sparse_wrapper<Scalar_T>& other);

    // Conversion to Arma Mat (implicit or explicit)
    // Conversion to const MatrixType reference
    operator const MatrixType&() const;
    // Conversion to MatrixType reference
    operator MatrixType&();

    // Attributes updated automatically by m_mat operations, accessors delegate directly
    // Set size
    void set_size(matrix_index_t rows, matrix_index_t cols);


    // Number of rows
    matrix_index_t nbr_rows() const;
    // Number of columns
    matrix_index_t nbr_cols() const;

    // Clear
    void clear();
    // Set to zero
    void zeros(matrix_index_t rows, matrix_index_t cols);
    // Set size then set to zero
    void zeros();
    // Set to identity
    void unit(matrix_index_t rows, matrix_index_t cols);

    // Element access
    // Element access
    Scalar_T& operator() (matrix_index_t i, matrix_index_t j);
    // Const element access
    const Scalar_T& operator() (matrix_index_t i, matrix_index_t j) const;

    // Operators
    // Add and assign
    arma_matrix_wrapper& operator+= (const arma_matrix_wrapper& other);
    // Subtract and assign
    arma_matrix_wrapper& operator-= (const arma_matrix_wrapper& other);
    // Multiply by scalar and assign
    arma_matrix_wrapper& operator*= (const Scalar_T& val);
    // Divide by scalar and assign
    arma_matrix_wrapper& operator/= (const Scalar_T& val);
    // Addition
    arma_matrix_wrapper operator+ (const arma_matrix_wrapper& other) const;
    // Subtraction
    arma_matrix_wrapper operator- (const arma_matrix_wrapper& other) const;
    // Multiplication
    arma_matrix_wrapper operator* (const arma_matrix_wrapper& other) const;
    // Unary negation
    arma_matrix_wrapper operator- () const;
    // Transpose
    arma_matrix_wrapper t() const;

    // Iterator support
    // Begin iterator
    auto begin();
    // End iterator
    auto end();
    // Begin const iterator
    auto begin() const;
    // End const iterator
    auto end() const;

    // Has infinity?
    bool has_inf() const;
    // Has NaN?
    bool has_nan() const;
    // Is finite?
    bool is_finite() const;

    // Trace
    Scalar_T trace() const;
    // Eigenvalues
    std::vector<std::complex<double>> eigenvalues() const;
    // Infinity norm
    Scalar_T norm_inf() const;
    // Squared Frobenius norm
    Scalar_T norm_frob2() const;
    // Is NaN?
    bool isnan() const;
    // Is infinite?
    bool isinf() const;
    // Number of non-zeros
    matrix_index_t nnz() const;

    // Inner product
    template< typename Result_Scalar_T, typename Other >
    Result_Scalar_T inner(const Other& other) const;

    // Trace of product: Trace(A*B) / Dim
    template< typename Other >
    Scalar_T trace_product(const Other& other) const;

    // Unary multivector mappings
    void similarity_transform(const std::vector<matrix_index_t>& perm, const std::vector<Scalar_T>& signs);
    void transpose_similarity_transform(const std::vector<matrix_index_t>& perm, const std::vector<Scalar_T>& signs);
    void involute();
    void reverse(index_t p, index_t q);

    // Kronecker matrix product
    arma_matrix_wrapper kron(const arma_matrix_wrapper& other) const;
    arma_matrix_wrapper mono_kron(const arma_matrix_wrapper& other) const;
    // Monomial matrix product
    arma_matrix_wrapper mono_prod(const arma_matrix_wrapper& other) const;
    // Mixed Kronecker matrix product: Dense x Sparse -> Dense (wrapper)
    template< typename Other_Scalar_T >
    arma_matrix_wrapper<Other_Scalar_T> kron(const arma_sparse_wrapper<Other_Scalar_T>& other) const;

    // Left Kronecker quotient
    template< typename RHS_T >
    RHS_T nork(const RHS_T& rhs, bool mono = true) const;

    // Output to stream
    friend std::ostream& operator<< <>(std::ostream& os, const arma_matrix_wrapper& m);

  private:
    // Helper to construct from raw arma mat
    arma_matrix_wrapper(const MatrixType& m);
    // Helper to construct from raw arma mat
    arma_matrix_wrapper(MatrixType&& m);
  };

  // Mixed op
  // Product of scalar and matrix wrapper
  template< typename Scalar_T >
  arma_matrix_wrapper<Scalar_T> operator* (Scalar_T s, const arma_matrix_wrapper<Scalar_T>& m);

  // Product of matrix wrapper and scalar
  template< typename Scalar_T >
  arma_matrix_wrapper<Scalar_T> operator* (const arma_matrix_wrapper<Scalar_T>& m, Scalar_T s);

  // Output to stream
  template< typename Scalar_T >
  std::ostream& operator<< (std::ostream& os, const arma_sparse_wrapper<Scalar_T>& m);

  /// Wrapper for Armadillo sparse matrix
  template< typename Scalar_T >
  class arma_sparse_wrapper :
  public matrix_base<arma_sparse_wrapper<Scalar_T>>
  {
  public:
    using MatrixType = arma::SpMat<Scalar_T>;

    using value_type = Scalar_T;
    using size_type = matrix_index_t;

    MatrixType m_mat;

    /// Default constructor
    arma_sparse_wrapper() = default;

    // Constructor from Armadillo SpMat
    explicit arma_sparse_wrapper(const MatrixType& m);

    // Constructor with size
    arma_sparse_wrapper(matrix_index_t rows, matrix_index_t cols);

    // Copy constructor
    arma_sparse_wrapper(const arma_sparse_wrapper& other);
    // Move constructor
    arma_sparse_wrapper(arma_sparse_wrapper&& other) noexcept(std::is_nothrow_move_constructible_v<Scalar_T>);

    // Copy assignment
    arma_sparse_wrapper& operator= (const arma_sparse_wrapper& other);
    // Move assignment
    arma_sparse_wrapper& operator= (arma_sparse_wrapper&& other) noexcept(std::is_nothrow_move_assignable_v<Scalar_T>);

    // Set size
    void set_size(matrix_index_t rows, matrix_index_t cols);



    // Clear
    void clear();
    // Set to zero
    void zeros(matrix_index_t rows, matrix_index_t cols);
    // Set to zero
    void zeros();

    // Set to identity
    void unit(matrix_index_t rows, matrix_index_t cols);

    using const_iterator = typename MatrixType::const_iterator;

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
    decltype(auto) operator() (matrix_index_t i, matrix_index_t j);

    // Add and assign
    arma_sparse_wrapper& operator+= (const arma_sparse_wrapper& other);
    // Multiply by sparse wrapper
    arma_sparse_wrapper operator* (const arma_sparse_wrapper& other) const;
    // Multiply by scalar and assign
    arma_sparse_wrapper& operator*= (const Scalar_T& val);
    // Mixed Kronecker matrix product: Sparse x Dense -> Dense (wrapper)
    template< typename Other_Scalar_T >
    arma_matrix_wrapper<Other_Scalar_T> kron(const arma_matrix_wrapper<Other_Scalar_T>& other) const;
    // Kronecker matrix product of sparse wrappers
    arma_sparse_wrapper kron(const arma_sparse_wrapper& other) const;
    arma_sparse_wrapper mono_kron(const arma_sparse_wrapper& other) const;
    // Monomial matrix product of sparse wrappers
    arma_sparse_wrapper mono_prod(const arma_sparse_wrapper& other) const;

    // Left Kronecker quotient
    template< typename RHS_T >
    RHS_T nork(const RHS_T& rhs, bool mono = true) const;

    // Output to stream
    friend std::ostream& operator<< <>(std::ostream& os, const arma_sparse_wrapper& m);

    // New Member Functions
    // Trace
    Scalar_T trace() const;
    // Eigenvalues
    std::vector<std::complex<double>> eigenvalues() const;
    // Infinity norm
    Scalar_T norm_inf() const;
    // Squared Frobenius norm
    Scalar_T norm_frob2() const;
    // Is NaN?
    bool isnan() const;
    // Is infinite?
    bool isinf() const;
    // Number of non-zeros
    matrix_index_t nnz() const;

    // Inner product
    template< typename Result_Scalar_T, typename Other >
    Result_Scalar_T inner(const Other& other) const;

    // Trace of product: Trace(A*B) / Dim
    template< typename Other >
    Scalar_T trace_product(const Other& other) const;

    // Unary multivector mappings
    void similarity_transform(const std::vector<matrix_index_t>& perm, const std::vector<Scalar_T>& signs);
    void transpose_similarity_transform(const std::vector<matrix_index_t>& perm, const std::vector<Scalar_T>& signs);
    void involute();
    void reverse(index_t p, index_t q);
  };

  // Product of scalar and sparse wrapper
  template< typename Scalar_T >
  arma_sparse_wrapper<Scalar_T> operator* (Scalar_T s, const arma_sparse_wrapper<Scalar_T>& m);

  // Product of sparse wrapper and scalar
  template< typename Scalar_T >
  arma_sparse_wrapper<Scalar_T> operator* (const arma_sparse_wrapper<Scalar_T>& m, Scalar_T s);

  // Sum of sparse wrappers
  template< typename Scalar_T >
  arma_sparse_wrapper<Scalar_T> operator+ (const arma_sparse_wrapper<Scalar_T>& lhs, const arma_sparse_wrapper<Scalar_T>& rhs);

  // Difference of sparse wrappers
  template< typename Scalar_T >
  arma_sparse_wrapper<Scalar_T> operator- (const arma_sparse_wrapper<Scalar_T>& lhs, const arma_sparse_wrapper<Scalar_T>& rhs);

} }
#endif // _GLUCAT_USE_ARMADILLO
#endif // _GLUCAT_MATRIX_ARMA_H
