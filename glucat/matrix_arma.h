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
#include "glucat/glucat_config.h"
#include "glucat/matrix.h"

#include <armadillo>

#include <type_traits>
#include <complex>
#include <vector>
#include <iostream>

namespace glucat { namespace matrix
{
  /// Matrix wrapper for Armadillo
  template< typename Scalar_T > class arma_matrix_wrapper; // Forward
  /// Sparse matrix wrapper for Armadillo
  template< typename Scalar_T > class arma_sparse_wrapper; // Forward

  /// Output to stream
  template< typename Scalar_T >
  auto operator<< (std::ostream& os, const arma_matrix_wrapper<Scalar_T>& m) -> std::ostream&;

  /// Wrapper for Armadillo matrix
  template< typename Scalar_T >
  class arma_matrix_wrapper :
  public matrix_impl_base<arma_matrix_wrapper<Scalar_T>>
  {
  public:
    using MatrixType = arma::Mat<Scalar_T>;

    using value_type = Scalar_T;
    using size_type = matrix_index_t;

    MatrixType m_mat;

    // Constructors
    /// Default constructor
    arma_matrix_wrapper() = default;

    /// Constructor with size
    arma_matrix_wrapper(matrix_index_t rows, matrix_index_t cols);

    /// Constructor from other matrix type
    template< typename Other_Matrix_T >
    explicit arma_matrix_wrapper(const Other_Matrix_T& other);

    /// Constructor from eigen_matrix_wrapper
    template< typename Other_Scalar_T >
    explicit arma_matrix_wrapper(const eigen_matrix_wrapper<Other_Scalar_T>& other);

    /// Constructor from eigen_sparse_wrapper
    template< typename Other_Scalar_T >
    explicit arma_matrix_wrapper(const eigen_sparse_wrapper<Other_Scalar_T>& other);

    /// Constructor from arma_sparse_wrapper
    template< typename Other_Scalar_T >
    explicit arma_matrix_wrapper(const arma_sparse_wrapper<Other_Scalar_T>& other);

    // Copy/Move
    /// Copy constructor
    arma_matrix_wrapper(const arma_matrix_wrapper& other);
    /// Move constructor
    arma_matrix_wrapper(arma_matrix_wrapper&& other) noexcept;

    /// Copy assignment
    auto operator= (const arma_matrix_wrapper& other) -> arma_matrix_wrapper&;
    /// Move assignment
    auto operator= (arma_matrix_wrapper&& other) noexcept -> arma_matrix_wrapper&;
    /// Assignment from sparse wrapper
    auto operator= (const arma_sparse_wrapper<Scalar_T>& other) -> arma_matrix_wrapper&;

    // Conversion to Arma Mat (implicit or explicit)
    /// Conversion to const MatrixType reference
    operator const MatrixType&() const;
    /// Conversion to MatrixType reference
    operator MatrixType&();

    // Attributes updated automatically by m_mat operations, accessors delegate directly
    /// Set size
    void set_size(matrix_index_t rows, matrix_index_t cols);
    /// Resize
    void resize(matrix_index_t rows, matrix_index_t cols, bool preserve = false);

    /// Number of rows
    auto nbr_rows() const -> matrix_index_t;
    /// Number of columns
    auto nbr_cols() const -> matrix_index_t;

    /// Clear
    void clear();
    /// Set to zero
    void zeros(matrix_index_t rows, matrix_index_t cols);
    /// Set size then set to zero
    void zeros();
    /// Set to identity
    void unit(matrix_index_t rows, matrix_index_t cols);

    // Element access
    /// Element access
    auto operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&;
    /// Const element access
    auto operator() (matrix_index_t i, matrix_index_t j) const -> const Scalar_T&;

    // Operators
    /// Add and assign
    auto operator+= (const arma_matrix_wrapper& other) -> arma_matrix_wrapper&;
    /// Subtract and assign
    auto operator-= (const arma_matrix_wrapper& other) -> arma_matrix_wrapper&;
    /// Multiply by scalar and assign
    auto operator*= (const Scalar_T& val) -> arma_matrix_wrapper&;
    /// Divide by scalar and assign
    auto operator/= (const Scalar_T& val) -> arma_matrix_wrapper&;
    /// Addition
    auto operator+ (const arma_matrix_wrapper& other) const -> arma_matrix_wrapper;
    /// Subtraction
    auto operator- (const arma_matrix_wrapper& other) const -> arma_matrix_wrapper;
    /// Multiplication
    auto operator* (const arma_matrix_wrapper& other) const -> arma_matrix_wrapper;
    /// Unary negation
    auto operator- () const -> arma_matrix_wrapper;
    /// Transpose
    auto t() const -> arma_matrix_wrapper;

    // Iterator support
    /// Begin iterator
    auto begin() { return m_mat.begin(); }
    /// End iterator
    auto end() { return m_mat.end(); }
    /// Begin const iterator
    auto begin() const { return m_mat.begin(); }
    /// End const iterator
    auto end() const { return m_mat.end(); }

    /// Has infinity?
    auto has_inf() const -> bool { return m_mat.has_inf(); }
    /// Has NaN?
    auto has_nan() const -> bool { return m_mat.has_nan(); }
    /// Is finite?
    auto is_finite() const -> bool { return m_mat.is_finite(); }

    /// Trace
    auto trace() const -> Scalar_T;
    /// Eigenvalues
    auto eigenvalues() const -> std::vector<std::complex<double>>;
    /// Infinity norm
    auto norm_inf() const;
    /// Squared Frobenius norm
    auto norm_frob2() const;
    /// Is NaN?
    auto isnan() const -> bool;
    /// Is infinite?
    auto isinf() const -> bool;
    /// Number of non-zeros
    auto nnz() const;

    /// Inner product
    template< typename Result_Scalar_T, typename Other >
    auto inner(const Other& other) const -> Result_Scalar_T;

    /// Kronecker matrix product
    auto kron(const arma_matrix_wrapper& other) const -> arma_matrix_wrapper;
    /// Mixed Kronecker matrix product: Dense x Sparse -> Dense (wrapper)
    template< typename Other_Scalar_T >
    auto kron(const arma_sparse_wrapper<Other_Scalar_T>& other) const -> arma_matrix_wrapper<Other_Scalar_T>;

    /// Left Kronecker quotient
    template< typename RHS_T >
    auto nork(const RHS_T& rhs, bool mono = true) const -> RHS_T;

    friend auto operator<< <>(std::ostream& os, const arma_matrix_wrapper& m) -> std::ostream&;

  private:
    /// Helper to construct from raw arma mat
    arma_matrix_wrapper(const MatrixType& m);
    /// Helper to construct from raw arma mat
    arma_matrix_wrapper(MatrixType&& m);
  };

  // Mixed op
  /// Product of scalar and matrix wrapper
  template< typename Scalar_T >
  auto operator* (Scalar_T s, const arma_matrix_wrapper<Scalar_T>& m) -> arma_matrix_wrapper<Scalar_T>
  {
    arma_matrix_wrapper<Scalar_T> result;
    result.m_mat = s * m.m_mat;
    return result;
  }

  /// Product of matrix wrapper and scalar
  template< typename Scalar_T >
  auto operator* (const arma_matrix_wrapper<Scalar_T>& m, Scalar_T s) -> arma_matrix_wrapper<Scalar_T>
  {  return s * m; }

  /// Output to stream
  template< typename Scalar_T >
  auto operator<< (std::ostream& os, const arma_sparse_wrapper<Scalar_T>& m) -> std::ostream&;

  /// Wrapper for Armadillo sparse matrix
  template< typename Scalar_T >
  class arma_sparse_wrapper :
  public matrix_impl_base<arma_sparse_wrapper<Scalar_T>>
  {
  public:
    using MatrixType = arma::SpMat<Scalar_T>;

    using value_type = Scalar_T;
    using size_type = matrix_index_t;

    MatrixType m_mat;

    /// Default constructor
    arma_sparse_wrapper() = default;

    /// Constructor from Armadillo SpMat
    explicit arma_sparse_wrapper(const MatrixType& m);

    /// Constructor with size
    arma_sparse_wrapper(matrix_index_t rows, matrix_index_t cols);

    /// Copy constructor
    arma_sparse_wrapper(const arma_sparse_wrapper& other);
    /// Move constructor
    arma_sparse_wrapper(arma_sparse_wrapper&& other) noexcept;

    /// Copy assignment
    auto operator= (const arma_sparse_wrapper& other) -> arma_sparse_wrapper&;
    /// Move assignment
    auto operator= (arma_sparse_wrapper&& other) noexcept -> arma_sparse_wrapper&;

    /// Set size
    void set_size(matrix_index_t rows, matrix_index_t cols);

    /// Resize
    void resize(matrix_index_t rows, matrix_index_t cols, bool preserve = false);

    /// Clear
    void clear();
    /// Set to zero
    void zeros(matrix_index_t rows, matrix_index_t cols);
    /// Set to zero
    void zeros();

    using const_iterator = typename MatrixType::const_iterator;

    /// Begin iterator
    auto begin() const -> const_iterator;
    /// End iterator
    auto end() const -> const_iterator;

    /// Number of rows
    auto nbr_rows() const -> matrix_index_t;
    /// Number of columns
    auto nbr_cols() const -> matrix_index_t;

    /// Const element access
    auto operator() (matrix_index_t i, matrix_index_t j) const -> Scalar_T;
    /// Element access
    auto operator() (matrix_index_t i, matrix_index_t j) -> auto;

    /// Add and assign
    auto operator+= (const arma_sparse_wrapper& other) -> arma_sparse_wrapper&;
    /// Multiply by sparse wrapper
    auto operator* (const arma_sparse_wrapper& other) const -> arma_sparse_wrapper;
    /// Multiply by scalar and assign
    auto operator*= (const Scalar_T& val) -> arma_sparse_wrapper&;
    /// Mixed Kronecker matrix product: Sparse x Dense -> Dense (wrapper)
    template< typename Other_Scalar_T >
    auto kron(const arma_matrix_wrapper<Other_Scalar_T>& other) const -> arma_matrix_wrapper<Other_Scalar_T>;

    /// Left Kronecker quotient
    template< typename RHS_T >
    auto nork(const RHS_T& rhs, bool mono = true) const -> RHS_T;

    friend auto operator<< <>(std::ostream& os, const arma_sparse_wrapper& m) -> std::ostream&;

    // New Member Functions
    /// Trace
    auto trace() const;
    /// Eigenvalues
    auto eigenvalues() const { throw std::runtime_error("Not implemented for sparse"); }
    /// Infinity norm
    auto norm_inf() const;
    /// Squared Frobenius norm
    auto norm_frob2() const;
    /// Is NaN?
    auto isnan() const -> bool;
    /// Is infinite?
    auto isinf() const -> bool;
    /// Number of non-zeros
    auto nnz() const;

    /// Inner product
    template< typename Result_Scalar_T, typename Other >
    auto inner(const Other& other) const -> Result_Scalar_T;
  };

  /// Product of scalar and sparse wrapper
  template< typename Scalar_T >
  auto operator* (Scalar_T s, const arma_sparse_wrapper<Scalar_T>& m) -> arma_sparse_wrapper<Scalar_T>
  {
    arma_sparse_wrapper<Scalar_T> result(m);
    result *= s;
    return result;
  }

  /// Product of sparse wrapper and scalar
  template< typename Scalar_T >
  auto operator* (const arma_sparse_wrapper<Scalar_T>& m, Scalar_T s) -> arma_sparse_wrapper<Scalar_T>
  { return s * m; }

  /// Sum of sparse wrappers
  template< typename Scalar_T >
  auto operator+ (const arma_sparse_wrapper<Scalar_T>& lhs, const arma_sparse_wrapper<Scalar_T>& rhs) -> arma_sparse_wrapper<Scalar_T>
  {
    arma_sparse_wrapper<Scalar_T> result(lhs);
    result += rhs;
    return result;
  }

  /// Difference of sparse wrappers
  template< typename Scalar_T >
  auto operator- (const arma_sparse_wrapper<Scalar_T>& lhs, const arma_sparse_wrapper<Scalar_T>& rhs) -> arma_sparse_wrapper<Scalar_T>
  {
    arma_sparse_wrapper<Scalar_T> result(lhs);
    result.m_mat -= rhs.m_mat; // Armadillo supports -=
    return result;
  }



  /// Armadillo support for float
  template<> struct is_arma_supported<float> : std::true_type {};
  /// Armadillo support for double
  template<> struct is_arma_supported<double> : std::true_type {};
  /// Armadillo support for complex float
  template<> struct is_arma_supported<std::complex<float>> : std::true_type {};
  /// Armadillo support for complex double
  template<> struct is_arma_supported<std::complex<double>> : std::true_type {};

  /// Matrix type selector specialization for Armadillo
  template< typename Scalar_T >
  struct matrix_type_selector<Scalar_T, true>
  {
    using type = arma_matrix_wrapper<Scalar_T>;
  };

  /// Sparse matrix type selector specialization for Armadillo
  template< typename Scalar_T >
  struct sparse_matrix_type_selector<Scalar_T, true>
  {
    using type = arma_sparse_wrapper<Scalar_T>;
  };
} }
#endif // _GLUCAT_USE_ARMADILLO
#endif // _GLUCAT_MATRIX_ARMA_H
