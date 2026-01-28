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

    /// Return const reference to derived class
    template< typename Derived >
    auto
    matrix_impl_base<Derived>::
    derived() const -> const Derived&
    { return static_cast<const Derived&>(*this); }

    /// Return reference to derived class
    template< typename Derived >
    auto
    matrix_impl_base<Derived>::
    derived() -> Derived&
    { return static_cast<Derived&>(*this); }

    // =========================================================================
    // Functions for Wrappers (to mimic Arma)
    // =========================================================================

    /// Solve
    template< typename Scalar_T >
    auto
    solve(eigen_matrix_wrapper<Scalar_T>& X, const eigen_matrix_wrapper<Scalar_T>& lhs, const eigen_matrix_wrapper<Scalar_T>& rhs, int opts = 0) -> bool
    {
      // Solve lhs*X = rhs
      // The matrix representation of a real Clifford algebra is always a real square matrix.
      if (lhs.nbr_rows() != lhs.nbr_cols())
        return false;

      auto lu = lhs.m_mat.fullPivLu();
      if (lu.isInvertible())
      {
        X.m_mat = lu.solve(rhs.m_mat);
        return true;
      }
      return false;
    }

#if defined(_GLUCAT_USE_ARMADILLO)

    /// Solve for arma_matrix_wrapper
    template< typename Scalar_T >
    auto
    solve(arma_matrix_wrapper<Scalar_T>& X, const arma_matrix_wrapper<Scalar_T>& lhs, const arma_matrix_wrapper<Scalar_T>& rhs, int opts = 0) -> bool
    {
      if (lhs.nbr_rows() != lhs.nbr_cols())
        return false;

      return arma::solve(X.m_mat, lhs.m_mat, rhs.m_mat, arma::solve_opts::no_approx);
    }
#endif


    // =========================================================================
    // matrix_impl_base Member Definitions
    // =========================================================================

    /// Generic classify_eigenvalues relies on eigenvalues() member
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



    // =========================================================================
    // eigen_matrix_wrapper Member Definitions
    // =========================================================================

    /// Number of rows
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    nbr_rows() const -> matrix_index_t
    { return static_cast<matrix_index_t>(m_mat.rows()); }

    /// Number of columns
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    nbr_cols() const -> matrix_index_t
    { return static_cast<matrix_index_t>(m_mat.cols()); }

    /// Armadillo constructor (rows, cols)
    template< typename Scalar_T >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.setZero();
    }

    /// Constructor from Eigen expressions (e.g. m * s)
    template< typename Scalar_T >
    template< typename Derived >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(const Eigen::MatrixBase<Derived>& other)
    { m_mat = other; }

    /// Generic Interop Constructor (e.g. from Armadillo matrix)
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

    /// Copy constructor
    template< typename Scalar_T >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(const eigen_matrix_wrapper<Scalar_T>& other)
    : m_mat(other.m_mat)
    { }

    /// Move constructor
    template< typename Scalar_T >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(eigen_matrix_wrapper<Scalar_T>&& other) noexcept
    : m_mat(std::move(other.m_mat))
    { }

    /// Assignment
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator= (const eigen_matrix_wrapper<Scalar_T>& other) -> eigen_matrix_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = other.m_mat;
      return *this;
    }

    /// Move Assignment
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator= (eigen_matrix_wrapper<Scalar_T>&& other) noexcept -> eigen_matrix_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = std::move(other.m_mat);
      return *this;
    }

    /// Constructor from eigen_sparse_wrapper
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

    /// Generic Interop Assignment
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
    /// Conversion to Armadillo Mat
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

    /// Constructor from Eigen
    template< typename Scalar_T >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(const MatrixType& m)
    : m_mat(m)
    { }

    /// Constructor from Eigen (move)
    template< typename Scalar_T >
    eigen_matrix_wrapper<Scalar_T>::
    eigen_matrix_wrapper(MatrixType&& m)
    : m_mat(std::move(m))
    { }

    /// Set size
    template< typename Scalar_T >
    void
    eigen_matrix_wrapper<Scalar_T>::
    set_size(matrix_index_t rows, matrix_index_t cols)
    { m_mat.resize(rows, cols); }

    /// Resize
    template< typename Scalar_T >
    void
    eigen_matrix_wrapper<Scalar_T>::
    resize(matrix_index_t rows, matrix_index_t cols, bool preserve)
    {
      // Preserve support can be simplified or removed if unused
      if (preserve)
        m_mat.conservativeResize(rows, cols);
      else
        m_mat.resize(rows, cols);
    }

    /// Clear
    template< typename Scalar_T >
    void
    eigen_matrix_wrapper<Scalar_T>::
    clear()
    { m_mat.setZero(); }

    /// Set to zero
    template< typename Scalar_T >
    void
    eigen_matrix_wrapper<Scalar_T>::
    zeros()
    { m_mat.setZero(); }

    /// Set size then set to zero
    template< typename Scalar_T >
    void
    eigen_matrix_wrapper<Scalar_T>::
    zeros(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      zeros();
    }

    /// Set to identity
    template< typename Scalar_T >
    void
    eigen_matrix_wrapper<Scalar_T>::
    unit(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.setIdentity();
    }

    /// Is finite?
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    is_finite() const -> bool
    { return m_mat.allFinite(); }

    /// Has NaN?
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    has_nan() const -> bool
    { return m_mat.hasNaN(); }

    /// Element access
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&
    { return m_mat(i, j); }

    /// Const element access
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) const -> const Scalar_T&
    { return m_mat(i, j); }

    /// Add and assign
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator+= (const eigen_matrix_wrapper<Scalar_T>& other) -> eigen_matrix_wrapper<Scalar_T>&
    {
      m_mat += other.m_mat;
      return *this;
    }

    /// Subtract and assign
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator-= (const eigen_matrix_wrapper<Scalar_T>& other) -> eigen_matrix_wrapper<Scalar_T>&
    {
      m_mat -= other.m_mat;
      return *this;
    }

    /// Multiply by scalar and assign
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator*= (const Scalar_T& val) -> eigen_matrix_wrapper<Scalar_T>&
    {
      m_mat *= val;
      return *this;
    }

    /// Divide by scalar and assign
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator/= (const Scalar_T& val) -> eigen_matrix_wrapper<Scalar_T>&
    {
      m_mat /= val;
      return *this;
    }

    /// Addition
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator+ (const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
    { return eigen_matrix_wrapper<Scalar_T>(m_mat + other.m_mat); }

    /// Subtraction
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator- (const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
    { return eigen_matrix_wrapper<Scalar_T>(m_mat - other.m_mat); }

    /// Matrix Multiplication
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator* (const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
    { return eigen_matrix_wrapper<Scalar_T>(m_mat * other.m_mat); }

    /// Unary -
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    operator- () const -> eigen_matrix_wrapper<Scalar_T>
    { return eigen_matrix_wrapper<Scalar_T>(-m_mat); }

    /// Transpose
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    t() const -> eigen_matrix_wrapper<Scalar_T>
    { return eigen_matrix_wrapper<Scalar_T>(m_mat.transpose()); }

    /// Output to stream
    template< typename Scalar_T >
    auto
    operator<< (std::ostream& os, const eigen_matrix_wrapper<Scalar_T>& m) -> std::ostream&
    { return os << m.m_mat; }

    // New Member Implementations
    // ========================

    /// Kronecker matrix product
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    kron(const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
    {
      using namespace Eigen;
      return eigen_matrix_wrapper<Scalar_T>(kroneckerProduct(m_mat, other.m_mat).eval());
    }

    /// Mixed Kronecker matrix product: Dense x Sparse -> Dense (wrapper)
    template< typename Scalar_T >
    template< typename Other_Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    kron(const eigen_sparse_wrapper<Other_Scalar_T>& other) const -> eigen_matrix_wrapper<Other_Scalar_T>
    {
      // Convert lhs (this, dense) to Other_Scalar_T
      eigen_matrix_wrapper<Other_Scalar_T> lhs_conv(nbr_rows(), nbr_cols());
      for (matrix_index_t i = 0; i < nbr_rows(); ++i)
        for (matrix_index_t j = 0; j < nbr_cols(); ++j)
          lhs_conv(i, j) = static_cast<Other_Scalar_T>((*this)(i, j));

      // Convert rhs (sparse) to compatible type (Dense)
      eigen_matrix_wrapper<Other_Scalar_T> rhs_dense(other.nbr_rows(), other.nbr_cols());
      for (matrix_index_t i = 0; i < other.nbr_rows(); ++i)
        for (matrix_index_t j = 0; j < other.nbr_cols(); ++j)
          rhs_dense(i, j) = static_cast<Other_Scalar_T>(other(i, j));

      return lhs_conv.kron(rhs_dense);
    }

    /// Trace
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    trace() const
    { return m_mat.trace(); }

    /// Infinity norm
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    norm_inf() const
    { return m_mat.cwiseAbs().rowwise().sum().maxCoeff(); }

    /// Squared Frobenius norm
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    norm_frob2() const
    { return m_mat.squaredNorm(); }

    /// Number of non-zeros
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    nnz() const
    { return (m_mat.array() != 0).count(); }

    /// Is NaN?
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    isnan() const -> bool
    { return m_mat.hasNaN(); }

    /// Is infinite?
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    isinf() const -> bool
    { return !m_mat.allFinite() && !m_mat.hasNaN(); }

    /// Eigenvalues
    template< typename Scalar_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    eigenvalues() const -> std::vector<std::complex<double>>
    {
      // Optimized for real matrices as per matrix representation
      if constexpr (std::is_arithmetic_v<Scalar_T> || std::is_same_v<Scalar_T, double> || std::is_same_v<Scalar_T, float> || std::is_same_v<Scalar_T, long double>)
      {
         Eigen::EigenSolver<typename eigen_matrix_wrapper<Scalar_T>::MatrixType> es(m_mat);
         const auto& E = es.eigenvalues();
         std::vector<std::complex<double>> result(E.size());
         for (int i = 0; i < E.size(); ++i)
           result[i] = std::complex<double>(E[i].real(), E[i].imag());
         return result;
      }
      else
      {
         // Fallback for complex or custom scalar types
         Eigen::MatrixXcd dmat(nbr_rows(), nbr_cols());
         for (matrix_index_t i = 0; i < nbr_rows(); ++i)
           for (matrix_index_t j = 0; j < nbr_cols(); ++j)
             dmat(i, j) = std::complex<double>(numeric_traits<Scalar_T>::to_double((*this)(i, j)), 0.0);

         Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(dmat);
         const auto& E = es.eigenvalues();
         std::vector<std::complex<double>> result(E.size());
         for (int i = 0; i < E.size(); ++i)
           result[i] = E[i];
         return result;
      }
    }

    /// Inner product
    template< typename Scalar_T >
    template< typename Result_Scalar_T, typename Other >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    inner(const Other& other) const -> Result_Scalar_T
    {
      Result_Scalar_T sum = Result_Scalar_T(0);
      for (matrix_index_t i = 0; i < nbr_rows(); ++i)
        for (matrix_index_t j = 0; j < nbr_cols(); ++j)
          sum += static_cast<Result_Scalar_T>((*this)(i, j)) * static_cast<Result_Scalar_T>(other(i, j));

      if (nbr_rows() == 0) return Result_Scalar_T(0);
      return sum / Result_Scalar_T(static_cast<double>(nbr_rows()));
    }

    /// Normalization of rotation K
    template< typename Scalar_T >
    template< typename RHS_T >
    auto
    eigen_matrix_wrapper<Scalar_T>::
    nork(const RHS_T& rhs, bool mono) const -> RHS_T
    {
        // Dense implementation (inverse kron logic)
        matrix_index_t blk_rows = rhs.nbr_rows() / (std::max)(matrix_index_t(1), nbr_rows());
        matrix_index_t blk_cols = rhs.nbr_cols() / (std::max)(matrix_index_t(1), nbr_cols());

        if (nbr_rows() == 0 || nbr_cols() == 0)
        {
          if constexpr (requires { RHS_T(blk_rows, blk_cols); })
            return RHS_T(blk_rows, blk_cols);
          RHS_T result;
          result.set_size(blk_rows, blk_cols);
          return result;
        }

        RHS_T result(blk_rows, blk_cols);
        result.zeros();

        // Loop over LHS (this) elements using dense indexing
        for (matrix_index_t r = 0; r < nbr_rows(); ++r)
        {
          for (matrix_index_t c = 0; c < nbr_cols(); ++c)
          {
            auto val = (*this)(r, c);
            if (val != Scalar_T(0))
            {
              matrix_index_t start_row = r * blk_rows;
              matrix_index_t start_col = c * blk_cols;
              for (matrix_index_t i = 0; i < blk_rows; ++i)
                for (matrix_index_t j = 0; j < blk_cols; ++j)
                  result(i, j) += static_cast<typename RHS_T::value_type>(val) * static_cast<typename RHS_T::value_type>(rhs(start_row + i, start_col + j));
            }
          }
        }

        // Normalize
        auto norm_sq = numeric_traits<typename RHS_T::value_type>::to_scalar_t(nbr_rows());
        if (norm_sq != numeric_traits<typename RHS_T::value_type>::to_scalar_t(1))
        {
          // We can use generic iteration or dense access depending on RHS_T
          // Assuming RHS_T supports (i,j) access or we should optimize if sparse
          // But result is typically dense if RHS is dense-ish or small.
          // Usually RHS_T is same as LHS_T (wrapper).
           if constexpr (requires { result(0, 0); })
            for (matrix_index_t i = 0; i < result.nbr_rows(); ++i)
              for (matrix_index_t j = 0; j < result.nbr_cols(); ++j)
                result(i, j) /= norm_sq;
        }
        return result;
    }


#if defined(_GLUCAT_USE_ARMADILLO)
    // =========================================================================
    // arma_matrix_wrapper Member Definitions
    // =========================================================================

    /// Number of rows
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    nbr_rows() const -> matrix_index_t
    { return m_mat.n_rows; }

    /// Number of columns
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    nbr_cols() const -> matrix_index_t
    { return m_mat.n_cols; }

    /// Constructor with size
    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.zeros();
    }

    /// Constructor from other matrix type
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

    /// Constructor from eigen_sparse_wrapper
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

    /// Constructor from arma_sparse_wrapper
    template< typename Scalar_T >
    template< typename Other_Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(const arma_sparse_wrapper<Other_Scalar_T>& other)
    {
      // Use efficient Armadillo conversion if possible
      m_mat = arma::conv_to<MatrixType>::from(other.m_mat);
    }

    /// Copy constructor
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

    /// Move constructor
    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(arma_matrix_wrapper<Scalar_T>&& other) noexcept
    : m_mat(std::move(other.m_mat))
    { }

    /// Copy assignment
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator= (const arma_matrix_wrapper<Scalar_T>& other) -> arma_matrix_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = other.m_mat;
      return *this;
    }

    /// Move assignment
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator= (arma_matrix_wrapper<Scalar_T>&& other) noexcept -> arma_matrix_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = std::move(other.m_mat);
      return *this;
    }

    /// Assignment from sparse wrapper
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator= (const arma_sparse_wrapper<Scalar_T>& other) -> arma_matrix_wrapper<Scalar_T>&
    {
      m_mat = other.m_mat;
      return *this;
    }

    /// Conversion to const MatrixType reference
    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    operator const MatrixType&() const
    { return m_mat; }

    /// Conversion to MatrixType reference
    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    operator MatrixType&()
    { return m_mat; }

    /// Set size
    template< typename Scalar_T >
    void
    arma_matrix_wrapper<Scalar_T>::
    set_size(matrix_index_t rows, matrix_index_t cols)
    { m_mat.set_size(rows, cols); }

    /// Resize
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

    /// Clear
    template< typename Scalar_T >
    void
    arma_matrix_wrapper<Scalar_T>::
    clear()
    { m_mat.zeros(); }

    /// Set to zero
    template< typename Scalar_T >
    void
    arma_matrix_wrapper<Scalar_T>::
    zeros(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      zeros();
    }

    /// Set size then set to zero
    template< typename Scalar_T >
    void
    arma_matrix_wrapper<Scalar_T>::
    zeros()
    { m_mat.zeros(); }

    /// Set to identity
    template< typename Scalar_T >
    void
    arma_matrix_wrapper<Scalar_T>::
    unit(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.eye();
    }

    /// Element access
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&
    { return m_mat(i, j); }

    /// Const element access
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) const -> const Scalar_T&
    { return m_mat(i, j); }

    /// Add and assign
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator+= (const arma_matrix_wrapper<Scalar_T>& other) -> arma_matrix_wrapper<Scalar_T>&
    {
      m_mat += other.m_mat;
      return *this;
    }

    /// Subtract and assign
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator-= (const arma_matrix_wrapper<Scalar_T>& other) -> arma_matrix_wrapper<Scalar_T>&
    {
      m_mat -= other.m_mat;
      return *this;
    }

    /// Multiply by scalar and assign
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator*= (const Scalar_T& val) -> arma_matrix_wrapper<Scalar_T>&
    {
      m_mat *= val;
      return *this;
    }

    /// Divide by scalar and assign
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator/= (const Scalar_T& val) -> arma_matrix_wrapper<Scalar_T>&
    {
      m_mat /= val;
      return *this;
    }

    /// Addition
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator+ (const arma_matrix_wrapper<Scalar_T>& other) const -> arma_matrix_wrapper<Scalar_T>
    { return arma_matrix_wrapper(MatrixType(m_mat + other.m_mat)); }

    /// Subtraction
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator- (const arma_matrix_wrapper<Scalar_T>& other) const -> arma_matrix_wrapper<Scalar_T>
    { return arma_matrix_wrapper(MatrixType(m_mat - other.m_mat)); }

    /// Multiplication
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator* (const arma_matrix_wrapper<Scalar_T>& other) const -> arma_matrix_wrapper<Scalar_T>
    {
      // Force evaluation to MatrixType (arma::Mat) to avoid resolving to generic template constructor with Glue
      MatrixType res_arma = m_mat * other.m_mat;
      return arma_matrix_wrapper(std::move(res_arma));
    }

    /// Unary negation
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    operator- () const -> arma_matrix_wrapper<Scalar_T>
    { return arma_matrix_wrapper(MatrixType(-m_mat)); }

    /// Transpose
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    t() const -> arma_matrix_wrapper<Scalar_T>
    { return arma_matrix_wrapper(MatrixType(m_mat.t())); }

    /// Output to stream
    template< typename Scalar_T >
    auto
    operator<< (std::ostream& os, const arma_matrix_wrapper<Scalar_T>& m) -> std::ostream&
    { return os << m.m_mat; }

    // New Member Implementations (moved from free functions)
    // ====================================================

    /// Kronecker matrix product
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    kron(const arma_matrix_wrapper<Scalar_T>& other) const -> arma_matrix_wrapper<Scalar_T>
    {
      arma_matrix_wrapper<Scalar_T> res;
      res.m_mat = arma::kron(m_mat, other.m_mat);
      return res;
    }

    /// Mixed Kronecker matrix product: Dense x Sparse -> Dense (wrapper)
    template< typename Scalar_T >
    template< typename Other_Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    kron(const arma_sparse_wrapper<Other_Scalar_T>& other) const -> arma_matrix_wrapper<Other_Scalar_T>
    {
        // Convert lhs (dense) to Other_Scalar_T
        arma_matrix_wrapper<Other_Scalar_T> lhs_conv;
        if constexpr (std::is_same_v<Scalar_T, Other_Scalar_T>)
        {
             lhs_conv = *this; 
        }
        else
        {
             // conversion logic
             lhs_conv.m_mat = arma::conv_to<typename arma_matrix_wrapper<Other_Scalar_T>::MatrixType>::from(m_mat);
        }

        // Convert sparse RHS to dense RHS
        arma_matrix_wrapper<Other_Scalar_T> rhs_dense(other);
        
        return lhs_conv.kron(rhs_dense);
    }

    /// Trace
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    trace() const -> Scalar_T
    { return arma::trace(m_mat); }

    /// Eigenvalues
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

    /// Is NaN?
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    isnan() const -> bool
    { return m_mat.has_nan(); }

    /// Is infinite?
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    isinf() const -> bool
    { return m_mat.has_inf(); }

    /// Infinity norm
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    norm_inf() const
    { return arma::norm(m_mat, "inf"); }

    /// Squared Frobenius norm
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    norm_frob2() const
    { return arma::accu(arma::square(m_mat)); }

    /// Number of non-zeros
    template< typename Scalar_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    nnz() const
    { return arma::accu(m_mat != 0); }

    /// Helper to construct from raw arma mat
    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(const MatrixType& m)
    : m_mat(m)
    { }

    /// Helper to construct from raw arma mat
    template< typename Scalar_T >
    arma_matrix_wrapper<Scalar_T>::
    arma_matrix_wrapper(MatrixType&& m)
    : m_mat(std::move(m))
    { }

    /// Inner product
    template< typename Scalar_T >
    template< typename Result_Scalar_T, typename Other >
    auto
    arma_matrix_wrapper<Scalar_T>::
    inner(const Other& other) const -> Result_Scalar_T
    {
      Result_Scalar_T sum = Result_Scalar_T(0);
      for (matrix_index_t i = 0; i < nbr_rows(); ++i)
        for (matrix_index_t j = 0; j < nbr_cols(); ++j)
          sum += static_cast<Result_Scalar_T>(m_mat(i, j)) * static_cast<Result_Scalar_T>(other(i, j));

      if (nbr_rows() == 0) return Result_Scalar_T(0);
      return sum / Result_Scalar_T(static_cast<double>(nbr_rows()));
    }

    /// Normalization of rotation K
    template< typename Scalar_T >
    template< typename RHS_T >
    auto
    arma_matrix_wrapper<Scalar_T>::
    nork(const RHS_T& rhs, bool mono) const -> RHS_T
    {
        // Dense implementation (inverse kron logic)
        matrix_index_t blk_rows = rhs.nbr_rows() / (std::max)(matrix_index_t(1), nbr_rows());
        matrix_index_t blk_cols = rhs.nbr_cols() / (std::max)(matrix_index_t(1), nbr_cols());

        if (nbr_rows() == 0 || nbr_cols() == 0)
        {
          if constexpr (requires { RHS_T(blk_rows, blk_cols); })
            return RHS_T(blk_rows, blk_cols);
          RHS_T result;
          result.set_size(blk_rows, blk_cols);
          return result;
        }

        RHS_T result(blk_rows, blk_cols);
        result.zeros();

        // Loop over LHS (this) elements using dense indexing
        for (matrix_index_t r = 0; r < nbr_rows(); ++r)
        {
          for (matrix_index_t c = 0; c < nbr_cols(); ++c)
          {
            auto val = (*this)(r, c);
            if (val != Scalar_T(0))
            {
              matrix_index_t start_row = r * blk_rows;
              matrix_index_t start_col = c * blk_cols;
              for (matrix_index_t i = 0; i < blk_rows; ++i)
                for (matrix_index_t j = 0; j < blk_cols; ++j)
                  result(i, j) += static_cast<typename RHS_T::value_type>(val) * static_cast<typename RHS_T::value_type>(rhs(start_row + i, start_col + j));
            }
          }
        }

        // Normalize
        auto norm_sq = numeric_traits<typename RHS_T::value_type>::to_scalar_t(nbr_rows());
        if (norm_sq != numeric_traits<typename RHS_T::value_type>::to_scalar_t(1))
        {
          if constexpr (requires { result(0, 0); })
            for (matrix_index_t i = 0; i < result.nbr_rows(); ++i)
              for (matrix_index_t j = 0; j < result.nbr_cols(); ++j)
                result(i, j) /= norm_sq;
        }
        return result;
    }
#endif

    // =========================================================================
    // eigen_sparse_wrapper Member Definitions
    // =========================================================================

    /// Number of rows
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    nbr_rows() const -> matrix_index_t
    { return static_cast<matrix_index_t>(m_mat.rows()); }

    /// Number of columns
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    nbr_cols() const -> matrix_index_t
    { return static_cast<matrix_index_t>(m_mat.cols()); }

    /// Constructor from Eigen Sparse Matrix (e.g. expression result)
    template< typename Scalar_T >
    eigen_sparse_wrapper<Scalar_T>::
    eigen_sparse_wrapper(const MatrixType& m)
    : m_mat(m)
    { }

    /// Armadillo/uBLAS/Generator style constructor support
    template< typename Scalar_T >
    eigen_sparse_wrapper<Scalar_T>::
    eigen_sparse_wrapper(matrix_index_t rows, matrix_index_t cols, matrix_index_t estimated_nnz)
    {
      set_size(rows, cols);
      if (estimated_nnz > 0) m_mat.reserve(estimated_nnz);
    }

    /// Copy/Move similar to dense
    template< typename Scalar_T >
    eigen_sparse_wrapper<Scalar_T>::
    eigen_sparse_wrapper(const eigen_sparse_wrapper<Scalar_T>& other)
    : m_mat(other.m_mat)
    { }

    /// Move constructor
    template< typename Scalar_T >
    eigen_sparse_wrapper<Scalar_T>::
    eigen_sparse_wrapper(eigen_sparse_wrapper<Scalar_T>&& other) noexcept
    : m_mat(std::move(other.m_mat))
    { }

    /// Copy assignment
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator= (const eigen_sparse_wrapper<Scalar_T>& other) -> eigen_sparse_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = other.m_mat;
      return *this;
    }

    /// Move assignment
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator= (eigen_sparse_wrapper<Scalar_T>&& other) noexcept -> eigen_sparse_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = std::move(other.m_mat);
      return *this;
    }

    /// Set size
    template< typename Scalar_T >
    void
    eigen_sparse_wrapper<Scalar_T>::
    set_size(matrix_index_t rows, matrix_index_t cols)
    { m_mat.resize(rows, cols); }

    /// Make writable
    template< typename Scalar_T >
    void
    eigen_sparse_wrapper<Scalar_T>::
    resize(matrix_index_t rows, matrix_index_t cols, bool preserve)
    {
      // preserve not directly supported in simple resize
      m_mat.resize(rows, cols);
    }

    /// Clear
    template< typename Scalar_T >
    void
    eigen_sparse_wrapper<Scalar_T>::
    clear()
    { m_mat.setZero(); }

    /// Set to zero
    template< typename Scalar_T >
    void
    eigen_sparse_wrapper<Scalar_T>::
    zeros()
    { m_mat.setZero(); }

    /// Set to zero
    template< typename Scalar_T >
    void
    eigen_sparse_wrapper<Scalar_T>::
    zeros(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      zeros();
    }

    /// Begin iterator
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    begin() const -> const_iterator
    { return const_iterator(&m_mat, true); }

    /// End iterator
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    end() const -> const_iterator
    { return const_iterator(&m_mat, false); }

    /// Const element access
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) const -> Scalar_T
    { return m_mat.coeff(i, j); }

    /// Element access
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&
    { return m_mat.coeffRef(i, j); }

    /// Add and assign
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator+= (const eigen_sparse_wrapper<Scalar_T>& other) -> eigen_sparse_wrapper<Scalar_T>&
    {
      m_mat += other.m_mat;
      return *this;
    }

    /// Subtract and assign
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator-= (const eigen_sparse_wrapper<Scalar_T>& other) -> eigen_sparse_wrapper<Scalar_T>&
    {
      m_mat -= other.m_mat;
      return *this;
    }

    /// Multiply by sparse wrapper
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator* (const eigen_sparse_wrapper<Scalar_T>& other) const -> eigen_sparse_wrapper<Scalar_T>
    { return eigen_sparse_wrapper(m_mat * other.m_mat); }

    /// Multiply by scalar and assign
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    operator*= (const Scalar_T& val) -> eigen_sparse_wrapper<Scalar_T>&
    {
      m_mat *= val;
      return *this;
    }

    /// Output to stream
    template< typename Scalar_T >
    auto
    operator<< (std::ostream& os, const eigen_sparse_wrapper<Scalar_T>& m) -> std::ostream&
    { return os << m.m_mat; }

    // const_iterator implementation
    /// Iterator support
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

    /// Constructor for begin()
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

    /// Check if end
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    is_end() const -> bool
    { return m_outer >= mp_mat->outerSize(); }

    /// Prefix increment
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    operator++ () -> const_iterator&
    {
      advance();
      return *this;
    }

    /// Inequality comparison
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    operator!= (const const_iterator& other) const -> bool
    {
      if (m_outer != other.m_outer)
        return true;
      if (m_outer >= mp_mat->outerSize())
        return false;
      return m_inner != other.m_inner;
    }

    /// Row index
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    row() const -> matrix_index_t
    { return m_inner.row(); }

    /// Column index
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    col() const -> matrix_index_t
    { return m_inner.col(); }

    /// Dereference
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    const_iterator::
    operator* () const -> Scalar_T
    { return m_inner.value(); }

    /// Is infinite?
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    isinf() const -> bool
    {
      for (int k = 0; k < m_mat.outerSize(); ++k)
        for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
          if (numeric_traits<Scalar_T>::isinf(it.value()))
            return true;
      return false;
    }

    /// Is NaN?
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    isnan() const -> bool
    {
      for (int k = 0; k < m_mat.outerSize(); ++k)
        for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
          if (numeric_traits<Scalar_T>::isnan(it.value()))
            return true;
      return false;
    }

    /// Trace
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    trace() const
    {
      Scalar_T sum = 0;
      for (int k = 0; k < m_mat.outerSize(); ++k)
        for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
          if (it.row() == it.col())
            sum += it.value();
      return sum;
    }

    /// Mixed Kronecker matrix product: Sparse x Dense -> Dense (wrapper)
    template< typename Scalar_T >
    template< typename Other_Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    kron(const eigen_matrix_wrapper<Other_Scalar_T>& other) const -> eigen_matrix_wrapper<Other_Scalar_T>
    {
      // Convert lhs (sparse) to compatible type (Dense of Other_Scalar_T)
      eigen_matrix_wrapper<Other_Scalar_T> lhs_dense(nbr_rows(), nbr_cols());
      for (matrix_index_t i = 0; i < nbr_rows(); ++i)
        for (matrix_index_t j = 0; j < nbr_cols(); ++j)
          lhs_dense(i, j) = static_cast<Other_Scalar_T>((*this)(i, j)); // using operator() const

      return lhs_dense.kron(other);
    }

    /// Kronecker matrix product of sparse wrappers
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    kron(const eigen_sparse_wrapper<Scalar_T>& other) const -> eigen_sparse_wrapper<Scalar_T>
    {
      eigen_sparse_wrapper<Scalar_T> result(nbr_rows() * other.nbr_rows(), nbr_cols() * other.nbr_cols());
      std::vector<Eigen::Triplet<Scalar_T>> triplets;
      // Iterate lhs
      for (int k = 0; k < m_mat.outerSize(); ++k)
      {
        for (typename MatrixType::InnerIterator itA(m_mat, k); itA; ++itA)
        {
          auto rA = itA.row();
          auto cA = itA.col();
          auto vA = itA.value();

          // Iterate rhs
          for (int l = 0; l < other.m_mat.outerSize(); ++l)
            for (typename MatrixType::InnerIterator itB(other.m_mat, l); itB; ++itB)
              triplets.emplace_back(rA * other.nbr_rows() + itB.row(), cA * other.nbr_cols() + itB.col(), vA * itB.value());
        }
      }
      result.m_mat.setFromTriplets(triplets.begin(), triplets.end());

      return result;
    }

    /// Infinity norm
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    norm_inf() const
    {
      Eigen::Vector<typename numeric_traits<Scalar_T>::real_t, Eigen::Dynamic> row_sums(nbr_rows());
      row_sums.setZero();
      for (int k = 0; k < m_mat.outerSize(); ++k)
        for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
          row_sums(it.row()) += numeric_traits<Scalar_T>::abs(it.value());
      return row_sums.maxCoeff();
    }

    /// Squared Frobenius norm
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    norm_frob2() const
    { return m_mat.squaredNorm(); }

    /// Number of non-zeros
    template< typename Scalar_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    nnz() const
    { return m_mat.nonZeros(); }

    /// Inner product
    template< typename Scalar_T >
    template< typename Result_Scalar_T, typename Other >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    inner(const Other& other) const -> Result_Scalar_T
    {
       Result_Scalar_T sum = Result_Scalar_T(0);
       for (int k=0; k < m_mat.outerSize(); ++k)
         for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
         {
            sum += static_cast<Result_Scalar_T>(it.value()) * static_cast<Result_Scalar_T>(other(it.row(), it.col()));
         }
       if (nbr_rows() == 0) return Result_Scalar_T(0);
       return sum / Result_Scalar_T(static_cast<double>(nbr_rows()));
    }

    /// Normalization of rotation K
    template< typename Scalar_T >
    template< typename RHS_T >
    auto
    eigen_sparse_wrapper<Scalar_T>::
    nork(const RHS_T& rhs, bool mono) const -> RHS_T
    {
       // Sparse implementation (signed_perm_nork logic)
       matrix_index_t blk_rows = rhs.nbr_rows() / (std::max)(matrix_index_t(1), nbr_rows());
       matrix_index_t blk_cols = rhs.nbr_cols() / (std::max)(matrix_index_t(1), nbr_cols());

       if (nbr_rows() == 0 || nbr_cols() == 0)
       {
         if constexpr (requires { RHS_T(blk_rows, blk_cols); })
           return RHS_T(blk_rows, blk_cols);
         RHS_T result;
         result.set_size(blk_rows, blk_cols);
         return result;
       }

       RHS_T result(blk_rows, blk_cols);
       result.zeros();

       // Sparse iterator optimization
       for (auto it = begin(); it != end(); ++it)
       {
         auto val = *it; // Value
         if (val != Scalar_T(0))
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

       // Normalize
       auto norm_sq = numeric_traits<typename RHS_T::value_type>::to_scalar_t(nbr_rows());
       if (norm_sq != numeric_traits<typename RHS_T::value_type>::to_scalar_t(1))
       {
         if constexpr (requires { result(0, 0); })
           for (matrix_index_t i = 0; i < result.nbr_rows(); ++i)
             for (matrix_index_t j = 0; j < result.nbr_cols(); ++j)
               result(i, j) /= norm_sq;
       }
       return result;
    }

    /// Identity matrix
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




#if defined(_GLUCAT_USE_ARMADILLO)
    // =========================================================================
    // arma_sparse_wrapper Member Definitions
    // =========================================================================

    /// Number of rows
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    nbr_rows() const -> matrix_index_t
    { return m_mat.n_rows; }

    /// Number of columns
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    nbr_cols() const -> matrix_index_t
    { return m_mat.n_cols; }

    /// Constructor from Armadillo SpMat
    template< typename Scalar_T >
    arma_sparse_wrapper<Scalar_T>::
    arma_sparse_wrapper(const MatrixType& m)
    : m_mat(m)
    { }

    /// Constructor with size
    template< typename Scalar_T >
    arma_sparse_wrapper<Scalar_T>::
    arma_sparse_wrapper(matrix_index_t rows, matrix_index_t cols)
    { set_size(rows, cols); }

    /// Copy constructor
    template< typename Scalar_T >
    arma_sparse_wrapper<Scalar_T>::
    arma_sparse_wrapper(const arma_sparse_wrapper<Scalar_T>& other)
    : m_mat(other.m_mat)
    { }

    /// Move constructor
    template< typename Scalar_T >
    arma_sparse_wrapper<Scalar_T>::
    arma_sparse_wrapper(arma_sparse_wrapper<Scalar_T>&& other) noexcept
    : m_mat(std::move(other.m_mat))
    { }

    /// Copy assignment
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator= (const arma_sparse_wrapper<Scalar_T>& other) -> arma_sparse_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = other.m_mat;
      return *this;
    }

    /// Move assignment
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator= (arma_sparse_wrapper<Scalar_T>&& other) noexcept -> arma_sparse_wrapper<Scalar_T>&
    {
      if (this != &other)
        m_mat = std::move(other.m_mat);
      return *this;
    }

    /// Set size
    template< typename Scalar_T >
    void
    arma_sparse_wrapper<Scalar_T>::
    set_size(matrix_index_t rows, matrix_index_t cols)
    { m_mat.set_size(rows, cols); }

    /// Resize
    template< typename Scalar_T >
    void
    arma_sparse_wrapper<Scalar_T>::
    resize(matrix_index_t rows, matrix_index_t cols, bool preserve)
    { m_mat.resize(rows, cols); }

    /// Clear
    template< typename Scalar_T >
    void
    arma_sparse_wrapper<Scalar_T>::
    clear()
    { m_mat.zeros(); }

    /// Set to zero
    template< typename Scalar_T >
    void
    arma_sparse_wrapper<Scalar_T>::
    zeros(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      zeros();
    }

    /// Set to zero
    template< typename Scalar_T >
    void
    arma_sparse_wrapper<Scalar_T>::
    zeros()
    { m_mat.zeros(); }

    /// Begin iterator
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    begin() const -> const_iterator
    { return m_mat.begin(); }

    /// End iterator
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    end() const -> const_iterator
    { return m_mat.end(); }

    /// Const element access
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) const -> Scalar_T
    { return m_mat(i, j); }

    /// Element access
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator() (matrix_index_t i, matrix_index_t j) -> auto
    { return m_mat(i, j); }

    /// Add and assign
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator+= (const arma_sparse_wrapper<Scalar_T>& other) -> arma_sparse_wrapper<Scalar_T>&
    {
      m_mat += other.m_mat;
      return *this;
    }

    /// Multiply by sparse wrapper
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator* (const arma_sparse_wrapper<Scalar_T>& other) const -> arma_sparse_wrapper<Scalar_T>
    { return arma_sparse_wrapper(m_mat * other.m_mat); }

    /// Multiply by scalar and assign
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    operator*= (const Scalar_T& val) -> arma_sparse_wrapper<Scalar_T>&
    {
      m_mat *= val;
      return *this;
    }

    /// Output to stream
    template< typename Scalar_T >
    auto
    operator<< (std::ostream& os, const arma_sparse_wrapper<Scalar_T>& m) -> std::ostream&
    { return os << m.m_mat; }

    // Armadillo Sparse Helper Member Functions
    /// Is infinite?
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    isinf() const -> bool
    { return m_mat.has_inf(); }

    /// Is NaN?
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    isnan() const -> bool
    { return m_mat.has_nan(); }

    /// Number of non-zeros
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    nnz() const
    { return m_mat.n_nonzero; }

    /// Infinity norm
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    norm_inf() const
    { return arma::norm(m_mat, "inf"); }

    /// Squared Frobenius norm
    template< typename Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    norm_frob2() const
    { return arma::accu(arma::square(m_mat)); }

    /// Trace
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

    /// Inner product
    template< typename Scalar_T >
    template< typename Result_Scalar_T, typename Other >
    auto
    arma_sparse_wrapper<Scalar_T>::
    inner(const Other& other) const -> Result_Scalar_T
    {
      Result_Scalar_T sum = Result_Scalar_T(0);
      for (auto it = m_mat.begin(); it != m_mat.end(); ++it)
      {
         sum += static_cast<Result_Scalar_T>(*it) * static_cast<Result_Scalar_T>(other(it.row(), it.col()));
      }
      if (nbr_rows() == 0) return Result_Scalar_T(0);
      return sum / Result_Scalar_T(static_cast<double>(nbr_rows()));
    }
    /// Mixed Kronecker matrix product: Sparse x Dense -> Dense (wrapper)
    template< typename Scalar_T >
    template< typename Other_Scalar_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    kron(const arma_matrix_wrapper<Other_Scalar_T>& other) const -> arma_matrix_wrapper<Other_Scalar_T>
    {
      // Convert lhs to compatible type (Dense Wrapper)
      arma_matrix_wrapper<Other_Scalar_T> lhs_dense(*this); // Will use constructor converting sparse to dense
      return lhs_dense.kron(other);
    }

    /// Normalization of rotation K
    template< typename Scalar_T >
    template< typename RHS_T >
    auto
    arma_sparse_wrapper<Scalar_T>::
    nork(const RHS_T& rhs, bool mono) const -> RHS_T
    {
       // Sparse implementation
       matrix_index_t blk_rows = rhs.nbr_rows() / (std::max)(matrix_index_t(1), nbr_rows());
       matrix_index_t blk_cols = rhs.nbr_cols() / (std::max)(matrix_index_t(1), nbr_cols());

       if (nbr_rows() == 0 || nbr_cols() == 0)
       {
         if constexpr (requires { RHS_T(blk_rows, blk_cols); })
           return RHS_T(blk_rows, blk_cols);
         RHS_T result;
         result.set_size(blk_rows, blk_cols);
         return result;
       }

       RHS_T result(blk_rows, blk_cols);
       result.zeros();

       // Sparse iterator optimization
       for (auto it = begin(); it != end(); ++it)
       {
         auto val = *it; // Value
         if (val != Scalar_T(0))
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

       // Normalize
       auto norm_sq = numeric_traits<typename RHS_T::value_type>::to_scalar_t(nbr_rows());
       if (norm_sq != numeric_traits<typename RHS_T::value_type>::to_scalar_t(1))
       {
         if constexpr (requires { result(0, 0); })
           for (matrix_index_t i = 0; i < result.nbr_rows(); ++i)
             for (matrix_index_t j = 0; j < result.nbr_cols(); ++j)
               result(i, j) /= norm_sq;
       }
       return result;
    }
#endif
  }
}

#endif // _GLUCAT_MATRIX_IMP_H
