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
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <set>
#include <unsupported/Eigen/KroneckerProduct>

#if defined(_GLUCAT_USE_QD)
  #include <qd/dd_real.h>
  #include <qd/qd_real.h>
#endif

namespace glucat
{
  using matrix_index_t = std::size_t;

  namespace matrix
  {
    // =========================================================================
    // Traits
    // =========================================================================
    /// Helper trait for complex check
    template< typename Scalar_T > struct is_complex_t : std::false_type {};
    template< typename Scalar_T > struct is_complex_t<std::complex<Scalar_T>> : std::true_type {};

    /// Type traits to detect wrappers (Primary Templates): is_eigen_sparse
    template< typename Matrix_T > struct is_eigen_sparse : std::false_type {};
    /// Type traits to detect wrappers (Primary Templates): is_eigen_dense
    template< typename Matrix_T > struct is_eigen_dense : std::false_type {};

    /// Matrix wrapper for Armadillo
    template< typename Scalar_T > class arma_matrix_wrapper; // Forward
    /// Sparse matrix wrapper for Eigen
    template< typename Scalar_T > class eigen_sparse_wrapper; // Forward
#if defined(_GLUCAT_USE_ARMADILLO)
    /// Sparse matrix wrapper for Armadillo
    template< typename Scalar_T > class arma_sparse_wrapper; // Forward
#endif

    // =========================================================================
    // matrix_impl_base (CRTP Pattern)
    // Base class providing member functions that delegate to the derived implementation
    // =========================================================================
    /// Base class providing member functions that delegate to the derived implementation (CRTP Pattern)
    template< typename Derived >
    class matrix_impl_base
    {
    public:
      /// Return const reference to derived class
      auto derived() const -> const Derived&;
      /// Return reference to derived class
      auto derived() -> Derived&;

      // Member functions delegating to namespace matrix implementation
      // Defined in matrix_imp.h to resolve circular dependency

      /// Generic classify_eigenvalues relies on eigenvalues() member
      auto classify_eigenvalues() const;

      /// Inner product of matrix and other
      template< typename Scalar_T, typename Other >
      auto inner(const Other& other) const;
    };

#if defined(_GLUCAT_USE_ARMADILLO)
    /// Forward declaration of arma_sparse_wrapper
    template< typename Scalar_T > class arma_sparse_wrapper;
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
      arma_matrix_wrapper<Scalar_T> res;
      res.m_mat = s * m.m_mat;
      return res;
    }

    /// Product of matrix wrapper and scalar
    template< typename Scalar_T >
    auto operator* (const arma_matrix_wrapper<Scalar_T>& m, Scalar_T s) -> arma_matrix_wrapper<Scalar_T>
    {
      return s * m;
    }

    /// Kronecker matrix product for arma_matrix_wrapper
    template< typename Scalar_T >
    auto kron(const arma_matrix_wrapper<Scalar_T>& lhs, const arma_matrix_wrapper<Scalar_T>& rhs) -> arma_matrix_wrapper<Scalar_T>
    {
      arma_matrix_wrapper<Scalar_T> res;
      res.m_mat = arma::kron(lhs.m_mat, rhs.m_mat);
      return res;
    }

    /// Mixed Kronecker matrix product: Sparse x Dense -> Dense (wrapper)
    template< typename LHS_Scalar_T, typename RHS_Scalar_T >
    auto kron(const arma_sparse_wrapper<LHS_Scalar_T>& lhs, const arma_matrix_wrapper<RHS_Scalar_T>& rhs) -> arma_matrix_wrapper<RHS_Scalar_T>;

    /// Mixed Kronecker matrix product: Dense x Sparse -> Dense (wrapper)
    template< typename LHS_Scalar_T, typename RHS_Scalar_T >
    auto kron(const arma_matrix_wrapper<LHS_Scalar_T>& lhs, const arma_sparse_wrapper<RHS_Scalar_T>& rhs) -> arma_matrix_wrapper<RHS_Scalar_T>;
#endif

    // =========================================================================
    // eigen_matrix_wrapper
    // =========================================================================
    template< typename Scalar_T > class eigen_matrix_wrapper; // Forward
    /// Output to stream
    template< typename Scalar_T >
    auto operator<< (std::ostream& os, const eigen_matrix_wrapper<Scalar_T>& m) -> std::ostream&;

    /// Wrapper for Eigen matrix
    template< typename Scalar_T >
    class eigen_matrix_wrapper :
    public matrix_impl_base<eigen_matrix_wrapper<Scalar_T>>
    {
    public:
      using MatrixType = Eigen::Matrix<Scalar_T, Eigen::Dynamic, Eigen::Dynamic>;

      using value_type = Scalar_T;
      using size_type = typename MatrixType::Index;

      MatrixType m_mat;

      // Constructors
      /// Default constructor
      eigen_matrix_wrapper() = default;

      /// Armadillo constructor (rows, cols)
      eigen_matrix_wrapper(matrix_index_t rows, matrix_index_t cols);

      /// Constructor from Eigen expressions (e.g. m * s)
      template< typename Derived >
      eigen_matrix_wrapper(const Eigen::MatrixBase<Derived>& other);

      /// Generic Interop Constructor (e.g. from Armadillo matrix)
      template< typename Other_Matrix_T >
      explicit eigen_matrix_wrapper(const Other_Matrix_T& other);

      /// Constructor from eigen_sparse_wrapper
      template< typename Other_Scalar_T >
      explicit eigen_matrix_wrapper(const eigen_sparse_wrapper<Other_Scalar_T>& other);

      // Copy constructor
      eigen_matrix_wrapper(const eigen_matrix_wrapper& other);

      // Move constructor
      eigen_matrix_wrapper(eigen_matrix_wrapper&& other) noexcept;

      // Assignment
      /// Copy assignment
      auto operator= (const eigen_matrix_wrapper& other) -> eigen_matrix_wrapper&;
      /// Move assignment
      auto operator= (eigen_matrix_wrapper&& other) noexcept -> eigen_matrix_wrapper&;

      /// Generic Interop Assignment
      template< typename Other_Matrix_T >
      auto operator= (const Other_Matrix_T& other) -> eigen_matrix_wrapper&;

      /// Conversion to Armadillo (if enabled)
#if defined(_GLUCAT_USE_ARMADILLO)
      /// Conversion to Armadillo Mat
      operator arma::Mat<Scalar_T>() const;
#endif

      /// Constructor from Eigen
      eigen_matrix_wrapper(const MatrixType& m);
      /// Constructor from Eigen (move)
      eigen_matrix_wrapper(MatrixType&& m);

      /// Set size
      void set_size(matrix_index_t rows, matrix_index_t cols);

      /// Resize
      void resize(matrix_index_t rows, matrix_index_t cols, bool preserve = false);

      // Helpers
      /// Number of rows
      auto nbr_rows() const -> matrix_index_t;
      /// Number of columns
      auto nbr_cols() const -> matrix_index_t;

      /// Clear
      void clear();

      /// Set to zero
      void zeros();
      /// Set size then set to zero
      void zeros(matrix_index_t rows, matrix_index_t cols);

      /// Set to identity
      void unit(matrix_index_t rows, matrix_index_t cols);

      /// Is finite?
      auto is_finite() const -> bool;
      /// Has NaN?
      auto has_nan() const -> bool;

      // Element access
      /// Element access
      auto operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&;
      /// Const element access
      auto operator() (matrix_index_t i, matrix_index_t j) const -> const Scalar_T&;

      // Operators
      /// Add and assign
      auto operator+= (const eigen_matrix_wrapper& other) -> eigen_matrix_wrapper&;
      /// Subtract and assign
      auto operator-= (const eigen_matrix_wrapper& other) -> eigen_matrix_wrapper&;
      /// Multiply by scalar and assign
      auto operator*= (const Scalar_T& val) -> eigen_matrix_wrapper&;
      /// Divide by scalar and assign
      auto operator/= (const Scalar_T& val) -> eigen_matrix_wrapper&;

      /// Addition
      auto operator+ (const eigen_matrix_wrapper& other) const -> eigen_matrix_wrapper;
      /// Subtraction
      auto operator- (const eigen_matrix_wrapper& other) const -> eigen_matrix_wrapper;

      /// Matrix Multiplication
      auto operator* (const eigen_matrix_wrapper& other) const -> eigen_matrix_wrapper;

      /// Unary -
      auto operator- () const -> eigen_matrix_wrapper;

      /// Transpose
      auto t() const -> eigen_matrix_wrapper;

      // New Member Functions (formerly free functions)
      /// Trace
      auto trace() const;
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

      friend auto operator<< <>(std::ostream& os, const eigen_matrix_wrapper& m) -> std::ostream&;
    };

    // Mixed op
    /// Product of scalar and matrix wrapper
    template< typename Scalar_T >
    auto operator* (Scalar_T s, const eigen_matrix_wrapper<Scalar_T>& m) -> eigen_matrix_wrapper<Scalar_T>
    {
      return eigen_matrix_wrapper<Scalar_T>(s * m.m_mat);
    }
    /// Product of matrix wrapper and scalar
    template< typename Scalar_T >
    auto operator* (const eigen_matrix_wrapper<Scalar_T>& m, Scalar_T s) -> eigen_matrix_wrapper<Scalar_T>
    {
      return eigen_matrix_wrapper<Scalar_T>(m.m_mat * s);
    }

    // =========================================================================
    // eigen_sparse_wrapper
    // =========================================================================
    /// Wrapper for Eigen sparse matrix
    template< typename Scalar_T > class eigen_sparse_wrapper;
    /// Output to stream
    template< typename Scalar_T >
    auto operator<< (std::ostream& os, const eigen_sparse_wrapper<Scalar_T>& m) -> std::ostream&;

    /// Wrapper for Eigen sparse matrix
    template< typename Scalar_T >
    class eigen_sparse_wrapper :
    public matrix_impl_base<eigen_sparse_wrapper<Scalar_T>>
    {
    public:
      using MatrixType = Eigen::SparseMatrix<Scalar_T>;

      using value_type = Scalar_T;
      using size_type = typename MatrixType::Index;

      MatrixType m_mat;

      eigen_sparse_wrapper() = default;

      /// Constructor from Eigen Sparse Matrix (e.g. expression result)
      explicit eigen_sparse_wrapper(const MatrixType& m);

      /// Armadillo/uBLAS/Generator style constructor support
      eigen_sparse_wrapper(matrix_index_t rows, matrix_index_t cols, matrix_index_t estimated_nnz = 0);

      /// Copy/Move similar to dense
      eigen_sparse_wrapper(const eigen_sparse_wrapper& other);

      eigen_sparse_wrapper(eigen_sparse_wrapper&& other) noexcept;

      auto operator= (const eigen_sparse_wrapper& other) -> eigen_sparse_wrapper&;

      auto operator= (eigen_sparse_wrapper&& other) noexcept -> eigen_sparse_wrapper&;

        void set_size(matrix_index_t rows, matrix_index_t cols);

      /// Make writable
        void resize(matrix_index_t rows, matrix_index_t cols, bool preserve = false);

      void clear();

      void zeros();
        void zeros(matrix_index_t rows, matrix_index_t cols);

      /// Iterator support
      class const_iterator
      {
      public:
        using InnerIterator = typename MatrixType::InnerIterator;

        const MatrixType* mp_mat;
        int m_outer;
        InnerIterator m_inner;

        /// Constructor for begin()
        const_iterator(const MatrixType* mat, bool start = true);

        /// Advance iterator
        void advance();

        /// Check if end
        auto is_end() const -> bool;
        /// Prefix increment
        auto operator++ () -> const_iterator&;

        /// Inequality comparison
        auto operator!= (const const_iterator& other) const -> bool;

        /// Row index
        auto row() const -> matrix_index_t;
        /// Column index
        auto col() const -> matrix_index_t;
        /// Dereference
        auto operator* () const -> Scalar_T;
      };

      // Iterator support
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
      auto operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&;

      /// Add and assign
      auto operator+= (const eigen_sparse_wrapper& other) -> eigen_sparse_wrapper&;
      /// Subtract and assign
      auto operator-= (const eigen_sparse_wrapper& other) -> eigen_sparse_wrapper&;

      /// Multiply by sparse wrapper
      auto operator* (const eigen_sparse_wrapper& other) const -> eigen_sparse_wrapper;

      /// Multiply by scalar and assign
      auto operator*= (const Scalar_T& val) -> eigen_sparse_wrapper&;

      // New Member Functions
      /// Trace
      auto trace() const; // Trace of sparse?
      /// Eigenvalues
      auto eigenvalues() const -> std::vector<std::complex<double>> { throw std::runtime_error("Not implemented for sparse"); } // Usually not computed directly on sparse
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

      friend auto operator<< <>(std::ostream& os, const eigen_sparse_wrapper& m) -> std::ostream&;
    };

    /// Product of sparse wrapper and scalar
    template< typename Scalar_T >
    auto operator* (const eigen_sparse_wrapper<Scalar_T>& m, Scalar_T s) -> eigen_sparse_wrapper<Scalar_T>
    {
      eigen_sparse_wrapper<Scalar_T> res(m);
      res *= s;
      return res;
    }

    /// Product of scalar and sparse wrapper
    template< typename Scalar_T >
    auto operator* (Scalar_T s, const eigen_sparse_wrapper<Scalar_T>& m) -> eigen_sparse_wrapper<Scalar_T>
    {
      return m * s;
    }

    /// Sum of sparse wrappers
    template< typename Scalar_T >
    auto operator+ (const eigen_sparse_wrapper<Scalar_T>& lhs, const eigen_sparse_wrapper<Scalar_T>& rhs) -> eigen_sparse_wrapper<Scalar_T>
    {
      eigen_sparse_wrapper<Scalar_T> res(lhs);
      res += rhs;
      return res;
    }

    /// Difference of sparse wrappers
    template< typename Scalar_T >
    auto operator- (const eigen_sparse_wrapper<Scalar_T>& lhs, const eigen_sparse_wrapper<Scalar_T>& rhs) -> eigen_sparse_wrapper<Scalar_T>
    {
      eigen_sparse_wrapper<Scalar_T> res(lhs);
      res -= rhs;
      return res;
    }

    /// Kronecker matrix product
    template< typename Scalar_T >
    auto kron(const eigen_matrix_wrapper<Scalar_T>& lhs, const eigen_matrix_wrapper<Scalar_T>& rhs) -> eigen_matrix_wrapper<Scalar_T>;

    /// Mixed Kronecker matrix product: Sparse x Dense -> Dense (wrapper)
    template< typename LHS_Scalar_T, typename RHS_Scalar_T >
    auto kron(const eigen_sparse_wrapper<LHS_Scalar_T>& lhs, const eigen_matrix_wrapper<RHS_Scalar_T>& rhs) -> eigen_matrix_wrapper<RHS_Scalar_T>;

    /// Mixed Kronecker matrix product: Dense x Sparse -> Dense (wrapper)
    template< typename LHS_Scalar_T, typename RHS_Scalar_T >
    auto kron(const eigen_matrix_wrapper<LHS_Scalar_T>& lhs, const eigen_sparse_wrapper<RHS_Scalar_T>& rhs) -> eigen_matrix_wrapper<RHS_Scalar_T>;

    /// Kronecker matrix product of sparse wrappers
    template< typename Scalar_T >
    auto kron(const eigen_sparse_wrapper<Scalar_T>& lhs, const eigen_sparse_wrapper<Scalar_T>& rhs) -> eigen_sparse_wrapper<Scalar_T>;

    // Helper Free Functions for Member Implementation


#if defined(_GLUCAT_USE_ARMADILLO)
    // =========================================================================
    // arma_sparse_wrapper (moved from matrix_imp.h)
    // =========================================================================
    /// Wrapper for Armadillo sparse matrix
    template< typename Scalar_T > class arma_sparse_wrapper;
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
    };

    /// Product of scalar and sparse wrapper
    template< typename Scalar_T >
    auto operator* (Scalar_T s, const arma_sparse_wrapper<Scalar_T>& m) -> arma_sparse_wrapper<Scalar_T>
    {
      arma_sparse_wrapper<Scalar_T> res(m);
      res *= s;
      return res;
    }

    /// Product of sparse wrapper and scalar
    template< typename Scalar_T >
    auto operator* (const arma_sparse_wrapper<Scalar_T>& m, Scalar_T s) -> arma_sparse_wrapper<Scalar_T>
    {
      return s * m;
    }

    /// Sum of sparse wrappers
    template< typename Scalar_T >
    auto operator+ (const arma_sparse_wrapper<Scalar_T>& lhs, const arma_sparse_wrapper<Scalar_T>& rhs) -> arma_sparse_wrapper<Scalar_T>
    {
      arma_sparse_wrapper<Scalar_T> res(lhs);
      res += rhs;
      return res;
    }

    /// Difference of sparse wrappers
    template< typename Scalar_T >
    auto operator- (const arma_sparse_wrapper<Scalar_T>& lhs, const arma_sparse_wrapper<Scalar_T>& rhs) -> arma_sparse_wrapper<Scalar_T>
    {
      arma_sparse_wrapper<Scalar_T> res(lhs);
      res.m_mat -= rhs.m_mat; // Armadillo supports -=
      return res;
    }
#endif

    // Traits Specializations
    /// Type trait specialization: is_eigen_dense for eigen_matrix_wrapper
    template< typename Matrix_T > struct is_eigen_dense<eigen_matrix_wrapper<Matrix_T>> : std::true_type {};
    /// Type trait specialization: is_eigen_sparse for eigen_sparse_wrapper
    template< typename Matrix_T > struct is_eigen_sparse<eigen_sparse_wrapper<Matrix_T>> : std::true_type {};

#if defined(_GLUCAT_USE_ARMADILLO)
    /// Type trait specialization: is_eigen_sparse for arma_sparse_wrapper
    template< typename Matrix_T > struct is_eigen_sparse<arma_sparse_wrapper<Matrix_T>> : std::true_type {};
#endif

    // Forward declarations of Wrappers (already defined but to match pattern if needed)

    /// Trait to determine if T is natively supported by Armadillo
    template< typename Scalar_T >
    struct is_arma_supported : std::false_type {};

#if defined(_GLUCAT_USE_ARMADILLO)
    /// Armadillo support for float
    template<> struct is_arma_supported<float> : std::true_type {};
    /// Armadillo support for double
    template<> struct is_arma_supported<double> : std::true_type {};
    /// Armadillo support for complex float
    template<> struct is_arma_supported<std::complex<float>> : std::true_type {};
    /// Armadillo support for complex double
    template<> struct is_arma_supported<std::complex<double>> : std::true_type {};

    // Explicitly NOTE: long double, dd_real, qd_real use the default false_type
    // and thus be dispatched to the Eigen wrapper.
#endif

    // Dense Selector
    /// Dense Matrix Type Selector
    template< typename Scalar_T, bool UseArma = is_arma_supported<Scalar_T>::value >
    struct matrix_type_selector
    {
      using type = eigen_matrix_wrapper<Scalar_T>;
    };

#if defined(_GLUCAT_USE_ARMADILLO)
    /// Matrix type selector specialization for Armadillo
    template< typename Scalar_T >
    struct matrix_type_selector<Scalar_T, true>
    {
      using type = arma_matrix_wrapper<Scalar_T>;
    };
#endif

    /// Matrix type selector
    template< typename Scalar_T >
    using matrix_t = typename matrix_type_selector<Scalar_T>::type;

    /// Sparse Matrix Selector
    template< typename Scalar_T, bool UseArma = is_arma_supported<Scalar_T>::value >
    struct sparse_matrix_type_selector
    {
      using type = eigen_sparse_wrapper<Scalar_T>;
    };

#if defined(_GLUCAT_USE_ARMADILLO)
    /// Sparse matrix type selector specialization for Armadillo
    template< typename Scalar_T >
    struct sparse_matrix_type_selector<Scalar_T, true>
    {
      using type = arma_sparse_wrapper<Scalar_T>;
    };
#endif

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

    // Core Operations as Free Functions
    /// Normalization of rotation K
    template< typename LHS_T, typename RHS_T >
    auto nork(const LHS_T& lhs, const RHS_T& rhs, const bool mono = true) -> const RHS_T;

    /// Signed permutation of normalization of rotation K
    template< typename LHS_T, typename RHS_T >
    auto signed_perm_nork(const LHS_T& lhs, const RHS_T& rhs) -> const RHS_T;

    /// Identity matrix
    template< typename Matrix_T >
    auto unit(const matrix_index_t dim) -> const Matrix_T;

    /// Inner product
    template< typename Scalar_T, typename LHS_T, typename RHS_T >
    auto inner(const LHS_T& lhs, const RHS_T& rhs) -> Scalar_T;

    // Legacy Nbr Rows/Cols - kept as free functions for generic templates

    /// Classification of eigenvalues of a matrix
    using eig_case_t = enum
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
      /// Is the matrix singular?
      bool m_is_singular = false;
      /// What kind of eigenvalues does the matrix contain?
      eig_case_t m_eig_case = safe_eigs;
      /// Argument such that exp(pi-m_safe_arg) lies between arguments of eigenvalues
      Scalar_T   m_safe_arg = Scalar_T(0);
    };
  }
}

#endif  // _GLUCAT_MATRIX_H
