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

namespace glucat {

  // =========================================================================
  // Traits
  // =========================================================================
  // Helper trait for complex check
  template<typename T> struct is_complex_t : std::false_type {};
  template<typename T> struct is_complex_t<std::complex<T>> : std::true_type {};

  // Type traits to detect wrappers (Primary Templates)
  template<typename T> struct is_eigen_sparse : std::false_type {};
  template<typename T> struct is_eigen_dense : std::false_type {};

  template<typename Scalar_T> class arma_matrix_wrapper; // Forward

  // =========================================================================
  // matrix_impl_base (CRTP Pattern)
  // Base class providing member functions that delegate to the derived implementation
  // =========================================================================
  template<typename Derived>
  class matrix_impl_base {
  public:
      const Derived& derived() const;
      Derived& derived();

      // Member functions delegating to namespace matrix implementation
      // Defined in matrix_imp.h to resolve circular dependency
      auto trace() const;
      auto eigenvalues() const;
      auto classify_eigenvalues() const;
      auto norm_inf() const;
      auto norm_frob2() const;
      auto isnan() const;
      auto isinf() const;
      auto nnz() const;

      template<typename Scalar_T, typename Other>
      auto inner(const Other& other) const;
  };

#if defined(_GLUCAT_USE_ARMADILLO)
  template<typename Scalar_T> std::ostream& operator<<(std::ostream& os, const arma_matrix_wrapper<Scalar_T>& m);

  template<typename Scalar_T>
  class arma_matrix_wrapper : public matrix_impl_base<arma_matrix_wrapper<Scalar_T>> {
  public:
    using MatrixType = arma::Mat<Scalar_T>;
    using elem_type = Scalar_T;
    using value_type = Scalar_T;
    using uword = arma::uword;
    using size_type = arma::uword;

    MatrixType m_mat;

    // Attributes
    uword n_rows = 0;
    uword n_cols = 0;
    uword n_elem = 0;

    // Constructors
    arma_matrix_wrapper() = default;

    arma_matrix_wrapper(uword rows, uword cols);

    template<typename Other_Matrix_T>
    explicit arma_matrix_wrapper(const Other_Matrix_T& other);

    // Copy/Move
    arma_matrix_wrapper(const arma_matrix_wrapper& other);
    arma_matrix_wrapper(arma_matrix_wrapper&& other) noexcept;

    arma_matrix_wrapper& operator=(const arma_matrix_wrapper& other);
    arma_matrix_wrapper& operator=(arma_matrix_wrapper&& other) noexcept;

    // Conversion to Arma Mat (implicit or explicit)
    operator const MatrixType&() const;
    operator MatrixType&();

    void update_attributes();

    void set_size(uword rows, uword cols);

    void resize(uword rows, uword cols, bool preserve = false);

    auto rows() const;
    auto cols() const;
    uword size1() const;
    uword size2() const;

    void clear();
    void zeros(uword rows, uword cols);
    void zeros();
    void eye(uword rows, uword cols);

    // Element access
    Scalar_T& operator()(uword i, uword j);
    const Scalar_T& operator()(uword i, uword j) const;

    // Operators
    arma_matrix_wrapper& operator+=(const arma_matrix_wrapper& other);
    arma_matrix_wrapper& operator-=(const arma_matrix_wrapper& other);
    arma_matrix_wrapper& operator*=(const Scalar_T& val);
    arma_matrix_wrapper& operator/=(const Scalar_T& val);

    arma_matrix_wrapper operator+(const arma_matrix_wrapper& other) const;
    arma_matrix_wrapper operator-(const arma_matrix_wrapper& other) const;
    arma_matrix_wrapper operator*(const arma_matrix_wrapper& other) const;
    arma_matrix_wrapper operator-() const;

    arma_matrix_wrapper t() const;

    friend std::ostream& operator<< <>(std::ostream& os, const arma_matrix_wrapper& m);

  private:
    // Helper to construct from raw arma mat
    arma_matrix_wrapper(const MatrixType& m);
    arma_matrix_wrapper(MatrixType&& m);
  };

  // Mixed op
  template<typename Scalar_T>
  arma_matrix_wrapper<Scalar_T> operator*(Scalar_T s, const arma_matrix_wrapper<Scalar_T>& m) {
      arma_matrix_wrapper<Scalar_T> res;
      res.m_mat = s * m.m_mat;
      res.update_attributes();
      return res;
  }
  template<typename Scalar_T>
  arma_matrix_wrapper<Scalar_T> operator*(const arma_matrix_wrapper<Scalar_T>& m, Scalar_T s) {
      return s * m;
  }

  // Kron for arma_matrix_wrapper
  template<typename T>
  arma_matrix_wrapper<T> kron(const arma_matrix_wrapper<T>& A, const arma_matrix_wrapper<T>& B) {
      arma_matrix_wrapper<T> res;
      res.m_mat = arma::kron(A.m_mat, B.m_mat);
      res.update_attributes();
      return res;
  }
#endif

  // =========================================================================
  // eigen_matrix_wrapper
  // =========================================================================
  template<typename Scalar_T> class eigen_matrix_wrapper; // Forward
  template<typename Scalar_T> std::ostream& operator<<(std::ostream& os, const eigen_matrix_wrapper<Scalar_T>& m);

  template<typename Scalar_T>
  class eigen_matrix_wrapper : public matrix_impl_base<eigen_matrix_wrapper<Scalar_T>> {
  public:
    using MatrixType = Eigen::Matrix<Scalar_T, Eigen::Dynamic, Eigen::Dynamic>;
    using elem_type = Scalar_T;
    using value_type = Scalar_T;
    using uword = std::size_t;
    using size_type = typename MatrixType::Index;

    MatrixType m_mat;

    // Attributes to match Armadillo
    uword n_rows = 0;
    uword n_cols = 0;
    uword n_elem = 0;

    // Constructors
    eigen_matrix_wrapper() = default;

    // Armadillo constructor (rows, cols)
    eigen_matrix_wrapper(uword rows, uword cols);

    // Constructor from Eigen expressions (e.g. m * s)
    template<typename Derived>
    eigen_matrix_wrapper(const Eigen::MatrixBase<Derived>& other);

    // Generic Interop Constructor (e.g. from Armadillo matrix)
    template<typename Other_Matrix_T>
    explicit eigen_matrix_wrapper(const Other_Matrix_T& other);

    // Copy constructor
    eigen_matrix_wrapper(const eigen_matrix_wrapper& other);

    // Move constructor
    eigen_matrix_wrapper(eigen_matrix_wrapper&& other) noexcept;

    // Assignment
    eigen_matrix_wrapper& operator=(const eigen_matrix_wrapper& other);
    eigen_matrix_wrapper& operator=(eigen_matrix_wrapper&& other) noexcept;

    // Generic Interop Assignment
    template<typename Other_Matrix_T>
    eigen_matrix_wrapper& operator=(const Other_Matrix_T& other);

    // Conversion to Armadillo (if enabled)
    #if defined(_GLUCAT_USE_ARMADILLO)
    operator arma::Mat<Scalar_T>() const;
    #endif

    // Constructor from Eigen
    eigen_matrix_wrapper(const MatrixType& m);
    eigen_matrix_wrapper(MatrixType&& m);

    // Sync attributes after resize
    void update_attributes();

    void set_size(uword rows, uword cols);

    void resize(uword rows, uword cols, bool preserve = false);

    // Helpers
    auto rows() const;
    auto cols() const;
    uword size1() const; // ublas compat
    uword size2() const; // ublas compat

    void clear();

    void zeros();
    void zeros(uword rows, uword cols);

    void eye(uword rows, uword cols);

    bool is_finite() const;
    bool has_nan() const;

    // Element access
    Scalar_T& operator()(uword i, uword j);
    const Scalar_T& operator()(uword i, uword j) const;

    // Operators
    eigen_matrix_wrapper& operator+=(const eigen_matrix_wrapper& other);
    eigen_matrix_wrapper& operator-=(const eigen_matrix_wrapper& other);
    eigen_matrix_wrapper& operator*=(const Scalar_T& val);
    eigen_matrix_wrapper& operator/=(const Scalar_T& val);

    eigen_matrix_wrapper operator+(const eigen_matrix_wrapper& other) const;
    eigen_matrix_wrapper operator-(const eigen_matrix_wrapper& other) const;

    // Matrix Multiplication
    eigen_matrix_wrapper operator*(const eigen_matrix_wrapper& other) const;

    // Unary -
    eigen_matrix_wrapper operator-() const;

    // Transpose
    eigen_matrix_wrapper t() const;

    friend std::ostream& operator<< <>(std::ostream& os, const eigen_matrix_wrapper& m);
  };

  // Mixed op
  template<typename Scalar_T>
  eigen_matrix_wrapper<Scalar_T> operator*(Scalar_T s, const eigen_matrix_wrapper<Scalar_T>& m) {
      return eigen_matrix_wrapper<Scalar_T>(s * m.m_mat);
  }
  template<typename Scalar_T>
  eigen_matrix_wrapper<Scalar_T> operator*(const eigen_matrix_wrapper<Scalar_T>& m, Scalar_T s) {
      return eigen_matrix_wrapper<Scalar_T>(m.m_mat * s);
  }

  // =========================================================================
  // eigen_sparse_wrapper
  // =========================================================================
  template<typename Scalar_T> class eigen_sparse_wrapper;
  template<typename Scalar_T> std::ostream& operator<<(std::ostream& os, const eigen_sparse_wrapper<Scalar_T>& m);

  template<typename Scalar_T>
  class eigen_sparse_wrapper : public matrix_impl_base<eigen_sparse_wrapper<Scalar_T>> {
  public:
    using MatrixType = Eigen::SparseMatrix<Scalar_T>;
    using elem_type = Scalar_T;
    using value_type = Scalar_T;
    using uword = std::size_t;
    using size_type = typename MatrixType::Index;

    MatrixType m_mat;

    uword n_rows = 0;
    uword n_cols = 0;
    uword n_nonzero = 0;

    eigen_sparse_wrapper() = default;

    // Armadillo/uBLAS/Generator style constructor support
    eigen_sparse_wrapper(uword rows, uword cols, uword estimated_nnz = 0);

    // Copy/Move similar to dense
    eigen_sparse_wrapper(const eigen_sparse_wrapper& other);

    eigen_sparse_wrapper(eigen_sparse_wrapper&& other) noexcept;

    eigen_sparse_wrapper& operator=(const eigen_sparse_wrapper& other);

    eigen_sparse_wrapper& operator=(eigen_sparse_wrapper&& other) noexcept;

    void set_size(uword rows, uword cols);

    // Make writable
    void resize(uword rows, uword cols, bool preserve = false);

    void clear();

    void zeros();
    void zeros(uword rows, uword cols);

    void update_attributes();

    // Iterator support
    class const_iterator {
    public:
        using InnerIterator = typename MatrixType::InnerIterator;

        const MatrixType* mp_mat;
        int m_outer;
        InnerIterator m_inner;

        // Constructor for begin()
        const_iterator(const MatrixType* mat, bool start = true);

        void advance();

        bool is_end() const;
        const_iterator& operator++();

        bool operator!=(const const_iterator& other) const;

        uword row() const;
        uword col() const;
        Scalar_T operator*() const;
    };

    const_iterator begin() const;
    const_iterator end() const;

    auto rows() const;
    auto cols() const;
    uword size1() const;
    uword size2() const;

    Scalar_T operator()(uword i, uword j) const;
    Scalar_T& operator()(uword i, uword j);

    eigen_sparse_wrapper& operator+=(const eigen_sparse_wrapper& other);

    eigen_sparse_wrapper operator*(const eigen_sparse_wrapper& other) const;

    eigen_sparse_wrapper& operator*=(const Scalar_T& val);

    friend std::ostream& operator<< <>(std::ostream& os, const eigen_sparse_wrapper& m);
  };

#if defined(_GLUCAT_USE_ARMADILLO)
  // =========================================================================
  // arma_sparse_wrapper (moved from matrix_imp.h)
  // =========================================================================
  template<typename Scalar_T> class arma_sparse_wrapper;
  template<typename Scalar_T> std::ostream& operator<<(std::ostream& os, const arma_sparse_wrapper<Scalar_T>& m);

  template<typename Scalar_T>
  class arma_sparse_wrapper : public matrix_impl_base<arma_sparse_wrapper<Scalar_T>> {
  public:
    using MatrixType = arma::SpMat<Scalar_T>;
    using elem_type = Scalar_T;
    using value_type = Scalar_T;
    using uword = arma::uword;
    using size_type = arma::uword;

    MatrixType m_mat;
    uword n_rows = 0;
    uword n_cols = 0;
    uword n_nonzero = 0;

    arma_sparse_wrapper() = default;

    arma_sparse_wrapper(uword rows, uword cols);

    // Copy
    arma_sparse_wrapper(const arma_sparse_wrapper& other);
    // Move
    arma_sparse_wrapper(arma_sparse_wrapper&& other) noexcept;

    arma_sparse_wrapper& operator=(const arma_sparse_wrapper& other);
    arma_sparse_wrapper& operator=(arma_sparse_wrapper&& other) noexcept;

    void set_size(uword rows, uword cols);

    void resize(uword rows, uword cols, bool preserve = false);

    void clear();
    void zeros(uword rows, uword cols);
    void zeros();

    void update_attributes();

    using const_iterator = typename MatrixType::const_iterator;

    const_iterator begin() const;
    const_iterator end() const;

    auto rows() const;
    auto cols() const;
    uword size1() const;
    uword size2() const;

    Scalar_T operator()(uword i, uword j) const;
    Scalar_T& operator()(uword i, uword j);

    arma_sparse_wrapper& operator+=(const arma_sparse_wrapper& other);
    arma_sparse_wrapper operator*(const arma_sparse_wrapper& other) const;
    arma_sparse_wrapper& operator*=(const Scalar_T& val);
    friend std::ostream& operator<< <>(std::ostream& os, const arma_sparse_wrapper& m);
  };
#endif

  // Traits Specializations
  template<typename T> struct is_eigen_dense<eigen_matrix_wrapper<T>> : std::true_type {};
  template<typename T> struct is_eigen_sparse<eigen_sparse_wrapper<T>> : std::true_type {};

  #if defined(_GLUCAT_USE_ARMADILLO)
  template<typename T> struct is_eigen_dense<arma_matrix_wrapper<T>> : std::true_type {};
  template<typename T> struct is_eigen_sparse<arma_sparse_wrapper<T>> : std::true_type {};
  #endif

  // Forward declarations of Wrappers (already defined but to match pattern if needed)

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

  // =========================================================================
  // Matrix Template Classes (Facade)
  // Named dense_matrix to avoid collision with namespace matrix (legacy)
  // =========================================================================

  template<typename Scalar_T>
  class dense_matrix : public matrix_type_selector<Scalar_T>::type
  {
  public:
      using Base = typename matrix_type_selector<Scalar_T>::type;
      using Base::Base; // Inherit constructors
      using Base::operator=;

      dense_matrix() = default;
      dense_matrix(const dense_matrix&) = default;
      dense_matrix(dense_matrix&&) = default;
      dense_matrix& operator=(const dense_matrix&) = default;
      dense_matrix& operator=(dense_matrix&&) = default;

      template<typename T>
      dense_matrix(const T& other) : Base(other) {}
  };

  template<typename Scalar_T>
  class sparse_matrix : public sparse_matrix_type_selector<Scalar_T>::type
  {
  public:
      using Base = typename sparse_matrix_type_selector<Scalar_T>::type;
      using Base::Base; // Inherit constructors
      using Base::operator=;

      sparse_matrix() = default;
      sparse_matrix(const sparse_matrix&) = default;
      sparse_matrix(sparse_matrix&&) = default;
      sparse_matrix& operator=(const sparse_matrix&) = default;
      sparse_matrix& operator=(sparse_matrix&&) = default;

      template<typename T>
      sparse_matrix(const T& other) : Base(other) {}
  };

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
    nnz(const Matrix_T& m); // size_type check?

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


    /// Matrix infinity norm (max row sum)
    template< typename Matrix_T >
    auto
    norm_inf(const Matrix_T& val) -> typename Matrix_T::elem_type;

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

#endif  // _GLUCAT_MATRIX_H
