#ifndef _GLUCAT_MATRIX_IMP_H
#define _GLUCAT_MATRIX_IMP_H

#include "glucat/errors.h"
#include "glucat/scalar.h"
#include "glucat/matrix.h"

#include <set>
#include <vector>
#include <complex>
#include <type_traits>
#include <cstdio>
#include <unsupported/Eigen/KroneckerProduct>
#include <iostream>
#include <algorithm>

namespace glucat
{
  namespace matrix
  {

    // =========================================================================
    // matrix_impl_base Member Definitions
    // =========================================================================

    template< typename Derived >
    auto
    matrix_impl_base<Derived>::
    derived() const -> const Derived&
    { return static_cast<const Derived&>(*this); }

    template< typename Derived >
    auto
    matrix_impl_base<Derived>::
    derived() -> Derived&
    { return static_cast<Derived&>(*this); }

    // =========================================================================
    // Functions for Wrappers (to mimic Arma)
    // =========================================================================

    /// Kron
    template< typename T >
    auto
    kron(const eigen_matrix_wrapper<T>& A, const eigen_matrix_wrapper<T>& B) -> eigen_matrix_wrapper<T>
    {
      using namespace Eigen;
      return eigen_matrix_wrapper<T>(kroneckerProduct(A.m_mat, B.m_mat).eval());
    }

    /// Solve
    template< typename T >
    auto
    solve(eigen_matrix_wrapper<T>& X, const eigen_matrix_wrapper<T>& A, const eigen_matrix_wrapper<T>& B, int opts = 0) -> bool
    {
      // Solve A*X = B
      // The matrix representation of a real Clifford algebra is always a real square matrix.
      if (A.nbr_rows() != A.nbr_cols())
        return false;

      auto lu = A.m_mat.fullPivLu();
      if (lu.isInvertible())
      {
        X.m_mat = lu.solve(B.m_mat);
        return true;
      }
      return false;
    }

#if defined(_GLUCAT_USE_ARMADILLO)

    /// Solve for arma_matrix_wrapper
    template< typename T >
    auto
    solve(arma_matrix_wrapper<T>& X, const arma_matrix_wrapper<T>& A, const arma_matrix_wrapper<T>& B, int opts = 0) -> bool
    {
      if (A.nbr_rows() != A.nbr_cols())
        return false;

      return arma::solve(X.m_mat, A.m_mat, B.m_mat, arma::solve_opts::no_approx);
    }
#endif

    /// Mixed kron: Sparse x Dense -> Dense (wrapper)
    template< typename T1, typename T2 >
    auto
    kron(const eigen_sparse_wrapper<T1>& A, const eigen_matrix_wrapper<T2>& B) -> eigen_matrix_wrapper<T2>
    {
      // Convert A to compatible type (Dense)
      eigen_matrix_wrapper<T2> A_dense(A.nbr_rows(), A.nbr_cols());
      for (matrix_index_t i = 0; i < A.nbr_rows(); ++i)
        for (matrix_index_t j = 0; j < A.nbr_cols(); ++j)
          A_dense(i, j) = static_cast<T2>(A(i, j));

      return kron(A_dense, B);
    }

    /// Dense x Sparse -> Dense (wrapper)
    template< typename T1, typename T2 >
    auto
    kron(const eigen_matrix_wrapper<T1>& A, const eigen_sparse_wrapper<T2>& B) -> eigen_matrix_wrapper<T2>
    {
      // Convert B to compatible type (Dense)
      eigen_matrix_wrapper<T1> B_dense(B.nbr_rows(), B.nbr_cols());
      for (matrix_index_t i = 0; i < B.nbr_rows(); ++i)
        for (matrix_index_t j = 0; j < B.nbr_cols(); ++j)
          B_dense(i, j) = static_cast<T1>(B(i, j));

      return kron(A, B_dense);
    }

    /// Sparse x Sparse -> Sparse
    template< typename T >
    auto
    kron(const eigen_sparse_wrapper<T>& A, const eigen_sparse_wrapper<T>& B) -> eigen_sparse_wrapper<T>
    {
      eigen_sparse_wrapper<T> result(A.nbr_rows() * B.nbr_rows(), A.nbr_cols() * B.nbr_cols());
      std::vector<Eigen::Triplet<T>> triplets;

      // Iterate A
      for (int k = 0; k < A.m_mat.outerSize(); ++k)
      {
        for (typename eigen_sparse_wrapper<T>::MatrixType::InnerIterator itA(A.m_mat, k); itA; ++itA)
        {
          auto rA = itA.row();
          auto cA = itA.col();
          auto vA = itA.value();

          // Iterate B
          for (int l = 0; l < B.m_mat.outerSize(); ++l)
            for (typename eigen_sparse_wrapper<T>::MatrixType::InnerIterator itB(B.m_mat, l); itB; ++itB)
              triplets.emplace_back(rA * B.nbr_rows() + itB.row(), cA * B.nbr_cols() + itB.col(), vA * itB.value());
        }
      }
      result.m_mat.setFromTriplets(triplets.begin(), triplets.end());

      return result;
    }

#if defined(_GLUCAT_USE_ARMADILLO)
    /// Mixed kron: Sparse (wrapper) x Dense (Armadillo) -> Dense (Armadillo)
    template< typename T1, typename T2 >
    auto
    kron(const eigen_sparse_wrapper<T1>& A, const arma::Mat<T2>& B) -> arma::Mat<T2>
    {
      // Convert A to arma::Mat (dense)
      arma_matrix_wrapper<T2> wrapper(A);
      return arma::kron(wrapper.m_mat, B);
    }

    /// Mixed kron: Sparse (wrapper) x Dense (Armadillo Wrapper) -> Dense (Armadillo Wrapper)
    template< typename T1, typename T2 >
    auto
    kron(const eigen_sparse_wrapper<T1>& A, const arma_matrix_wrapper<T2>& B) -> arma_matrix_wrapper<T2>
    {
      // Convert A to compatible type (Dense Armadillo Wrapper) efficiently
      arma_matrix_wrapper<T2> A_dense(A);
      return kron(A_dense, B);
    }

    /// Mixed kron: Sparse (arma wrapper) x Dense (arma wrapper) -> Dense (arma wrapper)
    /// This handles the case causing verify_glucat_kron.cpp to fail.
    template< typename T1, typename T2 >
    auto
    kron(const arma_sparse_wrapper<T1>& A, const arma_matrix_wrapper<T2>& B) -> arma_matrix_wrapper<T2>
    {
      // Convert A to compatible type (Dense Wrapper)
      arma_matrix_wrapper<T2> A_dense(A); // Will use constructor converting sparse to dense
      return kron(A_dense, B);
    }

    /// Mixed kron: Dense (arma wrapper) x Sparse (arma wrapper) -> Dense (arma wrapper)
    template< typename T1, typename T2 >
    auto
    kron(const arma_matrix_wrapper<T1>& A, const arma_sparse_wrapper<T2>& B) -> arma_matrix_wrapper<T2>
    {
      // Convert B to compatible type (Dense Wrapper)
      arma_matrix_wrapper<T1> B_dense(B);
      return kron(A, B_dense);
    }
#endif

    // =========================================================================
    // matrix_impl_base Member Definitions
    // =========================================================================

    /// classify_eigenvalues implementation moved from matrix namespace
    template< typename Derived >
    auto
    matrix_impl_base<Derived>::
    classify_eigenvalues() const
    {
      using Scalar_T = typename Derived::value_type;
      eig_genus<Derived> result;

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

    template< typename Derived >
    template< typename Scalar_T, typename Other >
    auto
    matrix_impl_base<Derived>::
    inner(const Other& other) const
    { return ::glucat::matrix::inner<Scalar_T>(derived(), other); }

    // =========================================================================
    // eigen_matrix_wrapper Member Definitions
    // =========================================================================

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    nbr_rows() const -> matrix_index_t
    { return static_cast<matrix_index_t>(m_mat.rows()); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    nbr_cols() const -> matrix_index_t
    { return static_cast<matrix_index_t>(m_mat.cols()); }

    template< typename Scalar_T >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.setZero();
    }

    template< typename Scalar_T >
    template< typename Derived >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(const Eigen::MatrixBase<Derived>& other)
    { m_mat = other; }

    template< typename Scalar_T >
    template< typename Other_Matrix_T >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(const Other_Matrix_T& other)
    {
      if constexpr (requires { other.nbr_rows(); })
      {
        set_size(other.nbr_rows(), other.nbr_cols());
        for (matrix_index_t i = 0; i < nbr_rows(); ++i)
          for (matrix_index_t j = 0; j < nbr_cols(); ++j)
            (*this)(i, j) = static_cast<Scalar_T>(other(i, j));
      }
      else
      {
        // Assume Eigen-compatible
        m_mat = other;
      }
    }

    template< typename Scalar_T >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(const eigen_matrix_wrapper<Scalar_T>& other)
    : m_mat(other.m_mat)
    { }

    template< typename Scalar_T >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(eigen_matrix_wrapper<Scalar_T>&& other) noexcept
    : m_mat(std::move(other.m_mat))
    { }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator= (const eigen_matrix_wrapper<Scalar_T>& other) -> eigen_matrix_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = other.m_mat;
      return *this;
    }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator= (eigen_matrix_wrapper<Scalar_T>&& other) noexcept -> eigen_matrix_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = std::move(other.m_mat);
      return *this;
    }

    template< typename Scalar_T >
    template< typename Other_Scalar_T >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(const eigen_sparse_wrapper<Other_Scalar_T>& other)
    {
      set_size(other.nbr_rows(), other.nbr_cols());
      m_mat.setZero();
      for (auto it = other.begin(); it != other.end(); ++it)
        (*this)(it.row(), it.col()) = static_cast<Scalar_T>(*it);
    }

    template< typename Scalar_T >
    template< typename Other_Matrix_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator= (const Other_Matrix_T& other) -> eigen_matrix_wrapper<Scalar_T>&
    {
      if constexpr (requires { other.nbr_rows(); })
      {
        set_size(other.nbr_rows(), other.nbr_cols());
        for (matrix_index_t i = 0; i < nbr_rows(); ++i)
          for (matrix_index_t j = 0; j < nbr_cols(); ++j)
            (*this)(i, j) = static_cast<Scalar_T>(other(i, j));
      }
      else
        m_mat = other;
      return *this;
    }

#if defined(_GLUCAT_USE_ARMADILLO)
    template< typename Scalar_T >
    eigen_matrix_wrapper<Scalar_T>::
    operator arma::Mat<Scalar_T>() const
    {
      arma::Mat<Scalar_T> result(nbr_rows(), nbr_cols());
      for (matrix_index_t i = 0; i < nbr_rows(); ++i)
        for (matrix_index_t j = 0; j < nbr_cols(); ++j)
          result(i, j) = (*this)(i, j);
      return result;
    }
#endif

    template< typename Scalar_T >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(const MatrixType& m)
    : m_mat(m)
    { }

    template< typename Scalar_T >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(MatrixType&& m)
    : m_mat(std::move(m))
    { }

    template< typename Scalar_T >
    void
    eigen_matrix_wrapper<Scalar_T>::
    set_size(matrix_index_t rows, matrix_index_t cols)
    { m_mat.resize(rows, cols); }

    template< typename Scalar_T >
    void
    eigen_matrix_wrapper<Scalar_T>::
    resize(matrix_index_t rows, matrix_index_t cols, bool preserve)
    {
      if (preserve)
        m_mat.conservativeResize(rows, cols);
      else
        m_mat.resize(rows, cols);
    }

    template< typename Scalar_T >
    void
    eigen_matrix_wrapper<Scalar_T>::
    clear()
    { m_mat.setZero(); }

    template< typename Scalar_T >
    void
    eigen_matrix_wrapper<Scalar_T>::
    zeros()
    { m_mat.setZero(); }

    template< typename Scalar_T >
    void
    eigen_matrix_wrapper<Scalar_T>::
    zeros(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.setZero();
    }

    template< typename Scalar_T >
    void
    eigen_matrix_wrapper<Scalar_T>::
    unit(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.setIdentity();
    }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    is_finite() const -> bool
    { return m_mat.allFinite(); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    has_nan() const -> bool
    { return m_mat.hasNaN(); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&
    { return m_mat(i, j); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) const -> const Scalar_T&
    { return m_mat(i, j); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator+= (const eigen_matrix_wrapper<Scalar_T>& other) -> eigen_matrix_wrapper<Scalar_T>&
    {
      m_mat += other.m_mat;
      return *this;
    }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator-= (const eigen_matrix_wrapper<Scalar_T>& other) -> eigen_matrix_wrapper<Scalar_T>&
    {
      m_mat -= other.m_mat;
      return *this;
    }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator*= (const Scalar_T& val) -> eigen_matrix_wrapper<Scalar_T>&
    {
      m_mat *= val;
      return *this;
    }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator/= (const Scalar_T& val) -> eigen_matrix_wrapper<Scalar_T>&
    {
      m_mat /= val;
      return *this;
    }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator+ (const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
    { return eigen_matrix_wrapper<Scalar_T>(m_mat + other.m_mat); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator- (const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
    { return eigen_matrix_wrapper<Scalar_T>(m_mat - other.m_mat); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator* (const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
    { return eigen_matrix_wrapper<Scalar_T>(m_mat * other.m_mat); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator- () const -> eigen_matrix_wrapper<Scalar_T>
    { return eigen_matrix_wrapper<Scalar_T>(-m_mat); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    t() const -> eigen_matrix_wrapper<Scalar_T>
    { return eigen_matrix_wrapper<Scalar_T>(m_mat.transpose()); }

    template< typename Scalar_T >
    auto
    operator<< (std::ostream& os, const eigen_matrix_wrapper<Scalar_T>& m) -> std::ostream&
    { return os << m.m_mat; }

    // New Member Implementations
    // ========================

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    trace() const
    { return m_mat.trace(); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    norm_inf() const
    { return m_mat.cwiseAbs().rowwise().sum().maxCoeff(); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    norm_frob2() const
    { return m_mat.squaredNorm(); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    nnz() const
    { return (m_mat.array() != 0).count(); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    isnan() const -> bool
    { return m_mat.hasNaN(); }

    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    isinf() const -> bool
    { return !m_mat.allFinite() && !m_mat.hasNaN(); }

    template< typename T >
    auto
    eigen_matrix_wrapper<T>::
    eigenvalues() const -> std::vector<std::complex<double>>
    {
      // If T is real
      if constexpr (std::is_arithmetic_v<T> || std::is_same_v<T, double> || std::is_same_v<T, float> || std::is_same_v<T, long double>)
      {
        Eigen::EigenSolver<typename eigen_matrix_wrapper<T>::MatrixType> es(m_mat);
        const auto& E = es.eigenvalues();
        std::vector<std::complex<double>> result(E.size());
        for (int i = 0; i < E.size(); ++i)
          result[i] = std::complex<double>(E[i].real(), E[i].imag());

        return result;
      }
      // TBD !! If T is complex
      else if constexpr (is_complex_t<T>::value)
      {
        Eigen::ComplexEigenSolver<typename eigen_matrix_wrapper<T>::MatrixType> es(m_mat);
        const auto& E = es.eigenvalues();
        std::vector<std::complex<double>> result(E.size());
        for (int i = 0; i < E.size(); ++i)
        {
          // complex cast to double
          result[i] = std::complex<double>(std::real(E[i]), std::imag(E[i]));
        }

        return result;
      }
      else
      {
        Eigen::MatrixXcd dmat(nbr_rows(), nbr_cols());
        for (matrix_index_t i = 0; i < nbr_rows(); ++i)
          for (matrix_index_t j = 0; j < nbr_cols(); ++j)
            dmat(i, j) = std::complex<double>(numeric_traits<T>::to_double((*this)(i, j)), 0.0);

        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(dmat);
        const auto& E = es.eigenvalues();
        std::vector<std::complex<double>> result(E.size());
        for (int i = 0; i < E.size(); ++i)
          result[i] = E[i];

        return result;
      }
    }

#if defined(_GLUCAT_USE_ARMADILLO)
    // =========================================================================
    // arma_matrix_wrapper Member Definitions
    // =========================================================================

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    nbr_rows() const -> matrix_index_t
    { return m_mat.n_rows; }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    nbr_cols() const -> matrix_index_t
    { return m_mat.n_cols; }

    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.zeros();
    }

    template< typename Scalar_T >
    template< typename Other_Matrix_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(const Other_Matrix_T& other)
    {
      if constexpr (requires { other.m_mat; })
      {
        // Wrapper to wrapper
        m_mat = arma::conv_to<MatrixType>::from(other.m_mat);
      }
      else
      {
        // Direct
        m_mat = arma::conv_to<MatrixType>::from(other);
      }
    }

    template< typename Scalar_T >
    template< typename Other_Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(const eigen_sparse_wrapper<Other_Scalar_T>& other)
    {
      set_size(other.nbr_rows(), other.nbr_cols());
      m_mat.zeros();
      for (auto it = other.begin(); it != other.end(); ++it)
        m_mat(it.row(), it.col()) = static_cast<Scalar_T>(*it);
    }

    template< typename Scalar_T >
    template< typename Other_Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(const arma_sparse_wrapper<Other_Scalar_T>& other)
    {
      // Use efficient Armadillo conversion if possible
      m_mat = arma::conv_to<MatrixType>::from(other.m_mat);
    }

    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(const arma_matrix_wrapper<Scalar_T>& other)
    : m_mat(other.m_mat)
    {
      if (nbr_rows() == 0 && other.nbr_rows() != 0) // Only warn if source was NOT zero but dest IS zero
        std::fprintf(stderr, "DEBUG: arma_matrix_wrapper COPY: Source %lux%lu -> Dest %lux%lu\n", (unsigned long)other.nbr_rows(), (unsigned long)other.nbr_cols(), (unsigned long)nbr_rows(), (unsigned long)nbr_cols());
      else if (nbr_rows() == 0) // Warn on any 0x0 copy?
        std::fprintf(stderr, "DEBUG: arma_matrix_wrapper COPY: Copying 0x0 matrix.\n");
    }

    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(arma_matrix_wrapper<Scalar_T>&& other) noexcept
    : m_mat(std::move(other.m_mat))
    { }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator= (const arma_matrix_wrapper<Scalar_T>& other) -> arma_matrix_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = other.m_mat;
      return *this;
    }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator= (arma_matrix_wrapper<Scalar_T>&& other) noexcept -> arma_matrix_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = std::move(other.m_mat);
      return *this;
    }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator= (const arma_sparse_wrapper<Scalar_T>& other) -> arma_matrix_wrapper<Scalar_T>&
    {
      m_mat = other.m_mat;
      return *this;
    }

    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    operator const MatrixType&() const
    { return m_mat; }

    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    operator MatrixType&()
    { return m_mat; }

    template< typename Scalar_T >
    void
    arma_matrix_wrapper<Scalar_T>::
    set_size(matrix_index_t rows, matrix_index_t cols)
    { m_mat.set_size(rows, cols); }

    template< typename Scalar_T >
    void
    arma_matrix_wrapper<Scalar_T>::
    resize(matrix_index_t rows, matrix_index_t cols, bool preserve)
    {
      if (preserve)
        m_mat.resize(rows, cols); // Arma resize preserves data
      else
        m_mat.set_size(rows, cols); // set_size does not preserve (faster)
    }

    template< typename Scalar_T >
    void
    arma_matrix_wrapper<Scalar_T>::
    clear()
    { m_mat.zeros(); }

    template< typename Scalar_T >
    void
    arma_matrix_wrapper<Scalar_T>::
    zeros(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.zeros();
    }

    template< typename Scalar_T >
    void
    arma_matrix_wrapper<Scalar_T>::
    zeros()
    { m_mat.zeros(); }

    template< typename Scalar_T >
    void
    arma_matrix_wrapper<Scalar_T>::
    unit(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.eye();
    }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&
    { return m_mat(i, j); }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) const -> const Scalar_T&
    { return m_mat(i, j); }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator+= (const arma_matrix_wrapper<Scalar_T>& other) -> arma_matrix_wrapper<Scalar_T>&
    {
      m_mat += other.m_mat;
      return *this;
    }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator-= (const arma_matrix_wrapper<Scalar_T>& other) -> arma_matrix_wrapper<Scalar_T>&
    {
      m_mat -= other.m_mat;
      return *this;
    }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator*= (const Scalar_T& val) -> arma_matrix_wrapper<Scalar_T>&
    {
      m_mat *= val;
      return *this;
    }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator/= (const Scalar_T& val) -> arma_matrix_wrapper<Scalar_T>&
    {
      m_mat /= val;
      return *this;
    }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator+ (const arma_matrix_wrapper<Scalar_T>& other) const -> arma_matrix_wrapper<Scalar_T>
    { return arma_matrix_wrapper(MatrixType(m_mat + other.m_mat)); }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator- (const arma_matrix_wrapper<Scalar_T>& other) const -> arma_matrix_wrapper<Scalar_T>
    { return arma_matrix_wrapper(MatrixType(m_mat - other.m_mat)); }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator* (const arma_matrix_wrapper<Scalar_T>& other) const -> arma_matrix_wrapper<Scalar_T>
    {
      // Force evaluation to MatrixType (arma::Mat) to avoid resolving to generic template constructor with Glue
      MatrixType res_arma = m_mat * other.m_mat;
      return arma_matrix_wrapper(std::move(res_arma));
    }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator- () const -> arma_matrix_wrapper<Scalar_T>
    { return arma_matrix_wrapper(MatrixType(-m_mat)); }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    t() const -> arma_matrix_wrapper<Scalar_T>
    { return arma_matrix_wrapper(MatrixType(m_mat.t())); }

    template< typename Scalar_T >
    auto
    operator<< (std::ostream& os, const arma_matrix_wrapper<Scalar_T>& m) -> std::ostream&
    { return os << m.m_mat; }

    // New Member Implementations (moved from free functions)
    // ====================================================

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    trace() const -> Scalar_T
    { return arma::trace(m_mat); }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    eigenvalues() const -> std::vector<std::complex<double>>
    {
      arma::cx_vec eigval;
      arma::eig_gen(eigval, m_mat);
      std::vector<std::complex<double>> result(eigval.n_elem);
      for (matrix_index_t i = 0; i < eigval.n_elem; ++i)
        result[i] = std::complex<double>(eigval[i].real(), eigval[i].imag());

      return result;
    }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    isnan() const -> bool
    { return m_mat.has_nan(); }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    isinf() const -> bool
    { return m_mat.has_inf(); }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    norm_inf() const
    { return arma::norm(m_mat, "inf"); }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    norm_frob2() const
    {
      // TBD !!
      if constexpr (is_complex_t<Scalar_T>::value)
      {
        auto n = arma::norm(m_mat, "fro");
        return n * n;
      }
      else
        return arma::accu(arma::square(m_mat));
    }

    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    nnz() const
    { return arma::accu(m_mat != 0); }

    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(const MatrixType& m)
    : m_mat(m)
    { }

    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(MatrixType&& m)
    : m_mat(std::move(m))
    { }
#endif

    // =========================================================================
    // eigen_sparse_wrapper Member Definitions
    // =========================================================================

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    nbr_rows() const -> matrix_index_t
    { return static_cast<matrix_index_t>(m_mat.rows()); }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    nbr_cols() const -> matrix_index_t
    { return static_cast<matrix_index_t>(m_mat.cols()); }

    template< typename Scalar_T >
    eigen_sparse_wrapper<Scalar_T>::
    eigen_sparse_wrapper(const MatrixType& m)
    : m_mat(m)
    { }

    template< typename Scalar_T >
    eigen_sparse_wrapper<Scalar_T>::
    eigen_sparse_wrapper(matrix_index_t rows, matrix_index_t cols, matrix_index_t estimated_nnz)
    {
      set_size(rows, cols);
      if (estimated_nnz > 0) m_mat.reserve(estimated_nnz);
    }

    template< typename Scalar_T >
    eigen_sparse_wrapper<Scalar_T>::
    eigen_sparse_wrapper(const eigen_sparse_wrapper<Scalar_T>& other)
    : m_mat(other.m_mat)
    { }

    template< typename Scalar_T >
    eigen_sparse_wrapper<Scalar_T>::
    eigen_sparse_wrapper(eigen_sparse_wrapper<Scalar_T>&& other) noexcept
    : m_mat(std::move(other.m_mat))
    { }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator= (const eigen_sparse_wrapper<Scalar_T>& other) -> eigen_sparse_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = other.m_mat;
      return *this;
    }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator= (eigen_sparse_wrapper<Scalar_T>&& other) noexcept -> eigen_sparse_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = std::move(other.m_mat);
      return *this;
    }

    template< typename Scalar_T >
    void
    eigen_sparse_wrapper<Scalar_T>::
    set_size(matrix_index_t rows, matrix_index_t cols)
    { m_mat.resize(rows, cols); }

    template< typename Scalar_T >
    void
    eigen_sparse_wrapper<Scalar_T>::
    resize(matrix_index_t rows, matrix_index_t cols, bool preserve)
    {
      // preserve not directly supported in simple resize
      m_mat.resize(rows, cols);
    }

    template< typename Scalar_T >
    void
    eigen_sparse_wrapper<Scalar_T>::
    clear()
    { m_mat.setZero(); }

    template< typename Scalar_T >
    void
    eigen_sparse_wrapper<Scalar_T>::
    zeros()
    { m_mat.setZero(); }

    template< typename Scalar_T >
    void
    eigen_sparse_wrapper<Scalar_T>::
    zeros(matrix_index_t rows, matrix_index_t cols)
    {
      // TBD !!
      set_size(rows, cols);
      m_mat.setZero();
    }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    begin() const -> const_iterator
    { return const_iterator(&m_mat, true); }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    end() const -> const_iterator
    { return const_iterator(&m_mat, false); }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) const -> Scalar_T
    { return m_mat.coeff(i, j); }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&
    { return m_mat.coeffRef(i, j); }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator+= (const eigen_sparse_wrapper<Scalar_T>& other) -> eigen_sparse_wrapper<Scalar_T>&
    {
      m_mat += other.m_mat;
      return *this;
    }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator-= (const eigen_sparse_wrapper<Scalar_T>& other) -> eigen_sparse_wrapper<Scalar_T>&
    {
      m_mat -= other.m_mat;
      return *this;
    }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator* (const eigen_sparse_wrapper<Scalar_T>& other) const -> eigen_sparse_wrapper<Scalar_T>
    { return eigen_sparse_wrapper(m_mat * other.m_mat); }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator*= (const Scalar_T& val) -> eigen_sparse_wrapper<Scalar_T>&
    {
      m_mat *= val;
      return *this;
    }

    template< typename Scalar_T >
    auto
    operator<< (std::ostream& os, const eigen_sparse_wrapper<Scalar_T>& m) -> std::ostream&
    { return os << m.m_mat; }

    // const_iterator implementation
    template< typename Scalar_T >
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    const_iterator(const MatrixType* mat, bool start)
    : mp_mat(mat),
      m_outer(0),
      m_inner(*mat, 0)
    {
      if (start)
      {
        if (mp_mat->outerSize() == 0)
        {
          m_outer = 0;
          return;
        }
        m_inner = InnerIterator(*mp_mat, 0);
        if (!m_inner) advance();
      }
      else
        m_outer = mp_mat->outerSize();
    }

    template< typename Scalar_T >
    void
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    advance()
    {
      if (m_inner)
        ++m_inner;
      while (!m_inner && m_outer < mp_mat->outerSize())
      {
        m_outer++;
        if (m_outer < mp_mat->outerSize())
          m_inner = InnerIterator(*mp_mat, m_outer);
      }
    }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    is_end() const -> bool
    { return m_outer >= mp_mat->outerSize(); }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    operator++ () -> const_iterator&
    {
      advance();
      return *this;
    }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    operator!= (const const_iterator& other) const -> bool
    {
      if (m_outer != other.m_outer) return true;
      if (m_outer >= mp_mat->outerSize()) return false;
      return m_inner != other.m_inner;
    }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    row() const -> matrix_index_t
    { return m_inner.row(); }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    col() const -> matrix_index_t
    { return m_inner.col(); }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    operator* () const -> Scalar_T
    { return m_inner.value(); }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    isinf() const -> bool
    {
      // TBD!!
      return false; // approximation
    }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    isnan() const -> bool
    {
      // TBD!!
      // Use generic check strictly or return false
      return false;
    }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    trace() const
    {
      Scalar_T sum = 0;
      for (int k = 0; k < m_mat.outerSize(); ++k)
        for (typename eigen_sparse_wrapper<Scalar_T>::MatrixType::InnerIterator it(m_mat, k); it; ++it)
          if (it.row() == it.col())
            sum += it.value();
      return sum;
    }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    norm_inf() const
    {
      Eigen::Vector<typename numeric_traits<Scalar_T>::real_t, Eigen::Dynamic> row_sums(nbr_rows());
      row_sums.setZero();
      for (int k = 0; k < m_mat.outerSize(); ++k)
        for (typename eigen_sparse_wrapper<Scalar_T>::MatrixType::InnerIterator it(m_mat, k); it; ++it)
          row_sums(it.row()) += numeric_traits<Scalar_T>::abs(it.value());
      return row_sums.maxCoeff();
    }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    norm_frob2() const
    { return m_mat.squaredNorm(); }

    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    nnz() const
    { return m_mat.nonZeros(); }

    template< typename Matrix_T >
    auto
    unit(const matrix_index_t dim) -> const Matrix_T
    {
      Matrix_T result(dim, dim);
      // Set to identity
      if constexpr (requires { result.unit(dim, dim); })
        result.unit(dim, dim);
      else if constexpr (requires { result.eye(dim, dim); })
        result.eye(dim, dim);
      else if constexpr (requires { result.setIdentity(); })
        result.setIdentity();
#if defined(_GLUCAT_USE_ARMADILLO)
      else if constexpr (requires { result.m_mat.eye(); })
        result.m_mat.eye();
#endif
      else if constexpr (requires { result.m_mat.setIdentity(); })
        result.m_mat.setIdentity();
      else
      {
        // Manual identity (may be slow for sparse if insertion not optimized)
        for (matrix_index_t i = 0; i < dim; ++i)
          result(i, i) = static_cast<typename Matrix_T::value_type>(1);
      }
      return result;
    }

    // nork / signed_perm_nork implementation
    template< typename LHS_T, typename RHS_T >
    auto
    signed_perm_nork(const LHS_T& lhs, const RHS_T& rhs) -> const RHS_T
    {
      matrix_index_t blk_rows = rhs.nbr_rows() / (std::max)(matrix_index_t(1), matrix_index_t(lhs.nbr_rows()));
      matrix_index_t blk_cols = rhs.nbr_cols() / (std::max)(matrix_index_t(1), matrix_index_t(lhs.nbr_cols()));

      if (lhs.nbr_rows() == 0 || lhs.nbr_cols() == 0)
      {
        if constexpr (requires { RHS_T(blk_rows, blk_cols); })
          return RHS_T(blk_rows, blk_cols);
        RHS_T result;
        result.set_size(blk_rows, blk_cols);
        return result;
      }

      RHS_T result(blk_rows, blk_cols);
      result.zeros();

      // Iterate over non-zero elements of LHS using iterator if available, or fallback
      if constexpr (requires { lhs.begin(); lhs.end(); })
      {
        // Sparse iterator optimization
        for (auto it = lhs.begin(); it != lhs.end(); ++it)
        {
          auto val = *it; // Value
          if (val != 0)
          {
            matrix_index_t r = it.row();
            matrix_index_t c = it.col();

            matrix_index_t start_row = r * blk_rows;
            matrix_index_t start_col = c * blk_cols;
            for (matrix_index_t i = 0; i < blk_rows; ++i)
              for (matrix_index_t j = 0; j < blk_cols; ++j)
                result(i, j) += static_cast<typename RHS_T::value_type>(val) * static_cast<typename RHS_T::value_type>(rhs(start_row + i, start_col + j));
          }
        }
      }
      else
      {
        // Fallback for types without sparse iterators (e.g. dense)
        for (matrix_index_t r = 0; r < lhs.nbr_rows(); ++r)
        {
          for (matrix_index_t c = 0; c < lhs.nbr_cols(); ++c)
          {
            auto val = lhs(r, c);
            if (val != 0)
            {
              matrix_index_t start_row = r * blk_rows;
              matrix_index_t start_col = c * blk_cols;
              for (matrix_index_t i = 0; i < blk_rows; ++i)
                for (matrix_index_t j = 0; j < blk_cols; ++j)
                  result(i, j) += static_cast<typename RHS_T::value_type>(val) * static_cast<typename RHS_T::value_type>(rhs(start_row + i, start_col + j));
            }
          }
        }
      }

      // Normalize by dividing by Frobenius norm squared of LHS.
      // For a signed permutation matrix of size N, norm_frob2 is N (assuming entries are +/- 1).
      // This matches legacy implementation which used lhs.size1().
      // Use to_scalar_t to safely convert n_rows (size_t/long) to Scalar (potentially qd_real)
      auto norm_sq = numeric_traits<typename RHS_T::value_type>::to_scalar_t(lhs.nbr_rows());
      if (norm_sq != numeric_traits<typename RHS_T::value_type>::to_scalar_t(1))
      {
        // If not 1, we must scale result
        if constexpr (requires { result(0, 0); })
          for (matrix_index_t i = 0; i < result.nbr_rows(); ++i)
            for (matrix_index_t j = 0; j < result.nbr_cols(); ++j)
              result(i, j) /= norm_sq;
      }
      return result;
    }

    template< typename LHS_T, typename RHS_T >
    auto
    nork(const LHS_T& lhs, const RHS_T& rhs, const bool mono) -> const RHS_T
    { return signed_perm_nork(lhs, rhs); }

    template< typename Scalar_T, typename LHS_T, typename RHS_T >
    auto
    inner(const LHS_T& lhs, const RHS_T& rhs) -> Scalar_T
    {
      Scalar_T sum = Scalar_T(0);
      // Generic implementation assuming n_rows/n_cols and operator()
      // Note: For sparse LHS, this iteration is inefficient (O(N^2) instead of O(NNZ))
      // but it is correct and works for all wrappers.
      // Optimization TODO: Use iterators if LHS is sparse.
      if constexpr (requires { lhs.begin(); lhs.end(); })
      {
        // Sparse iterator optimization
        for (auto it = lhs.begin(); it != lhs.end(); ++it)
        {
          matrix_index_t r = it.row();
          matrix_index_t c = it.col();
          auto val = *it;

          // For inner product: sum( lhs(i,j) * rhs(i,j) )
          // We only need to visit non-zeros of LHS.
          // Assuming RHS has efficient random access or is dense.
          sum += static_cast<Scalar_T>(val) * static_cast<Scalar_T>(rhs(r, c));
        }
      }
      else
      {
        // Fallback
        for (size_t i = 0; i < lhs.nbr_rows(); ++i)
          for (size_t j = 0; j < lhs.nbr_cols(); ++j)
            sum += static_cast<Scalar_T>(lhs(i, j)) * static_cast<Scalar_T>(rhs(i, j));
      }
      if (lhs.nbr_rows() == 0)
        return Scalar_T(0);
      return sum / Scalar_T(double(lhs.nbr_rows()));
    }

#if defined(_GLUCAT_USE_ARMADILLO)
    // =========================================================================
    // arma_sparse_wrapper Member Definitions
    // =========================================================================

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    nbr_rows() const -> matrix_index_t
    { return m_mat.n_rows; }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    nbr_cols() const -> matrix_index_t
    { return m_mat.n_cols; }

    template< typename Scalar_T >
    arma_sparse_wrapper<Scalar_T>::
    arma_sparse_wrapper(const MatrixType& m)
    : m_mat(m)
    { }

    template< typename Scalar_T >
    arma_sparse_wrapper<Scalar_T>::
    arma_sparse_wrapper(matrix_index_t rows, matrix_index_t cols)
    { set_size(rows, cols); }

    template< typename Scalar_T >
    arma_sparse_wrapper<Scalar_T>::
    arma_sparse_wrapper(const arma_sparse_wrapper<Scalar_T>& other)
    : m_mat(other.m_mat)
    { }

    template< typename Scalar_T >
    arma_sparse_wrapper<Scalar_T>::
    arma_sparse_wrapper(arma_sparse_wrapper<Scalar_T>&& other) noexcept
    : m_mat(std::move(other.m_mat))
    { }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator= (const arma_sparse_wrapper<Scalar_T>& other) -> arma_sparse_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = other.m_mat;
      return *this;
    }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator= (arma_sparse_wrapper<Scalar_T>&& other) noexcept -> arma_sparse_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = std::move(other.m_mat);
      return *this;
    }

    template< typename Scalar_T >
    void
    arma_sparse_wrapper<Scalar_T>::
    set_size(matrix_index_t rows, matrix_index_t cols)
    { m_mat.set_size(rows, cols); }

    template< typename Scalar_T >
    void
    arma_sparse_wrapper<Scalar_T>::
    resize(matrix_index_t rows, matrix_index_t cols, bool preserve)
    { m_mat.resize(rows, cols); }

    template< typename Scalar_T >
    void
    arma_sparse_wrapper<Scalar_T>::
    clear()
    { m_mat.zeros(); }

    template< typename Scalar_T >
    void
    arma_sparse_wrapper<Scalar_T>::
    zeros(matrix_index_t rows, matrix_index_t cols)
    {
      // TBD !!
      set_size(rows, cols);
      m_mat.zeros();
    }

    template< typename Scalar_T >
    void
    arma_sparse_wrapper<Scalar_T>::
    zeros()
    {
      // TBD !!
      m_mat.zeros();
    }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    begin() const -> const_iterator
    { return m_mat.begin(); }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    end() const -> const_iterator
    { return m_mat.end(); }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) const -> Scalar_T
    { return m_mat(i, j); }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) -> auto
    { return m_mat(i, j); }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator+= (const arma_sparse_wrapper<Scalar_T>& other) -> arma_sparse_wrapper<Scalar_T>&
    {
      m_mat += other.m_mat;
      return *this;
    }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator* (const arma_sparse_wrapper<Scalar_T>& other) const -> arma_sparse_wrapper<Scalar_T>
    { return arma_sparse_wrapper(m_mat * other.m_mat); }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator*= (const Scalar_T& val) -> arma_sparse_wrapper<Scalar_T>&
    {
      m_mat *= val;
      return *this;
    }

    template< typename Scalar_T >
    auto
    operator<< (std::ostream& os, const arma_sparse_wrapper<Scalar_T>& m) -> std::ostream&
    { return os << m.m_mat; }

    // Armadillo Sparse Helper Member Functions
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    isinf() const -> bool
    { return m_mat.has_inf(); }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    isnan() const -> bool
    { return m_mat.has_nan(); }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    nnz() const
    { return m_mat.n_nonzero; }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    norm_inf() const
    { return arma::norm(m_mat, "inf"); }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    norm_frob2() const
    {
      // TBD !!
      if constexpr (is_complex_t<Scalar_T>::value)
      {
        auto n = arma::norm(m_mat, "fro");
        return n * n;
      }
      else
      {
        auto n = arma::norm(m_mat, "fro");
        return n * n;
      }
    }

    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    trace() const
    {
      // generic trace for sparse
      Scalar_T sum = 0;
      // Armadillo SpMat iterators
      for (auto it = m_mat.begin(); it != m_mat.end(); ++it)
        if (it.row() == it.col())
          sum += *it;
      return sum;
    }
#endif
  }
}

#endif // _GLUCAT_MATRIX_IMP_H
