#ifndef _GLUCAT_MATRIX_EIGEN_IMP_H
#define _GLUCAT_MATRIX_EIGEN_IMP_H
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

namespace glucat { namespace matrix
{
  // =========================================================================
  // Functions for Wrappers (to mimic Armadillo)
  // =========================================================================

  /**
   * @brief Solve
   * @details
   */
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

  // =========================================================================
  // eigen_matrix_wrapper Member Definitions
  // =========================================================================

  /**
   * @brief Number of rows
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  nbr_rows() const -> matrix_index_t
  { return static_cast<matrix_index_t>(m_mat.rows()); }

  /**
   * @brief Number of columns
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  nbr_cols() const -> matrix_index_t
  { return static_cast<matrix_index_t>(m_mat.cols()); }

  /**
   * @brief Armadillo constructor (rows, cols)
   * @details
   */
  template< typename Scalar_T >
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(matrix_index_t rows, matrix_index_t cols)
  {
    set_size(rows, cols);
    m_mat.setZero();
  }

  /**
   * @brief Constructor from Eigen expressions (e.g. m * s)
   * @details
   */
  template< typename Scalar_T >
  template< typename Derived >
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(const Eigen::MatrixBase<Derived>& other)
  { m_mat = other; }

  /**
   * @brief Generic Interop Constructor (e.g. from Armadillo matrix)
   * @details
   */
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

  /**
   * @brief Copy constructor
   * @details
   */
  template< typename Scalar_T >
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(const eigen_matrix_wrapper<Scalar_T>& other)
  : m_mat(other.m_mat)
  { }

  /**
   * @brief Move constructor
   * @details
   */
  template< typename Scalar_T >
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(eigen_matrix_wrapper<Scalar_T>&& other) noexcept
  : m_mat(std::move(other.m_mat))
  { }

  /**
   * @brief Assignment
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  operator= (const eigen_matrix_wrapper<Scalar_T>& other) -> eigen_matrix_wrapper<Scalar_T>&
  {
    if (this != &other)
      m_mat = other.m_mat;
    return *this;
  }

  /**
   * @brief Move Assignment
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  operator= (eigen_matrix_wrapper<Scalar_T>&& other) noexcept -> eigen_matrix_wrapper<Scalar_T>&
  {
    if (this != &other)
      m_mat = std::move(other.m_mat);
    return *this;
  }

  /**
   * @brief Constructor from eigen_sparse_wrapper
   * @details
   */
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

  /**
   * @brief Generic Interop Assignment
   * @details
   */
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

  /**
   * @brief Constructor from Eigen
   * @details
   */
  template< typename Scalar_T >
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(const MatrixType& m)
  : m_mat(m)
  { }

  /**
   * @brief Constructor from Eigen (move)
   * @details
   */
  template< typename Scalar_T >
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(MatrixType&& m)
  : m_mat(std::move(m))
  { }

  /**
   * @brief Set size
   * @details
   */
  template< typename Scalar_T >
  void
  eigen_matrix_wrapper<Scalar_T>::
  set_size(matrix_index_t rows, matrix_index_t cols)
  { m_mat.resize(rows, cols); }

  /**
   * @brief Resize
   * @details
   */
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

  /**
   * @brief Clear
   * @details
   */
  template< typename Scalar_T >
  void
  eigen_matrix_wrapper<Scalar_T>::
  clear()
  { m_mat.setZero(); }

  /**
   * @brief Set to zero
   * @details
   */
  template< typename Scalar_T >
  void
  eigen_matrix_wrapper<Scalar_T>::
  zeros()
  { m_mat.setZero(); }

  /**
   * @brief Set size then set to zero
   * @details
   */
  template< typename Scalar_T >
  void
  eigen_matrix_wrapper<Scalar_T>::
  zeros(matrix_index_t rows, matrix_index_t cols)
  {
    set_size(rows, cols);
    zeros();
  }

  /**
   * @brief Set to identity
   * @details
   */
  template< typename Scalar_T >
  void
  eigen_matrix_wrapper<Scalar_T>::
  unit(matrix_index_t rows, matrix_index_t cols)
  {
    set_size(rows, cols);
    m_mat.setIdentity();
  }

  /**
   * @brief Is finite?
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  is_finite() const -> bool
  { return m_mat.allFinite(); }

  /**
   * @brief Has NaN?
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  has_nan() const -> bool
  { return m_mat.hasNaN(); }

  /**
   * @brief Element access
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&
  { return m_mat(i, j); }

  /**
   * @brief Const element access
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  operator() (matrix_index_t i, matrix_index_t j) const -> const Scalar_T&
  { return m_mat(i, j); }

  /**
   * @brief Add and assign
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  operator+= (const eigen_matrix_wrapper<Scalar_T>& other) -> eigen_matrix_wrapper<Scalar_T>&
  {
    m_mat += other.m_mat;
    return *this;
  }

  /**
   * @brief Subtract and assign
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  operator-= (const eigen_matrix_wrapper<Scalar_T>& other) -> eigen_matrix_wrapper<Scalar_T>&
  {
    m_mat -= other.m_mat;
    return *this;
  }

  /**
   * @brief Multiply by scalar and assign
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  operator*= (const Scalar_T& val) -> eigen_matrix_wrapper<Scalar_T>&
  {
    m_mat *= val;
    return *this;
  }

  /**
   * @brief Divide by scalar and assign
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  operator/= (const Scalar_T& val) -> eigen_matrix_wrapper<Scalar_T>&
  {
    m_mat /= val;
    return *this;
  }

  /**
   * @brief Addition
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  operator+ (const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
  { return eigen_matrix_wrapper<Scalar_T>(m_mat + other.m_mat); }

  /**
   * @brief Subtraction
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  operator- (const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
  { return eigen_matrix_wrapper<Scalar_T>(m_mat - other.m_mat); }

  /**
   * @brief Matrix Multiplication
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  operator* (const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
  { return eigen_matrix_wrapper<Scalar_T>(m_mat * other.m_mat); }

  /**
   * @brief Unary -
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  operator- () const -> eigen_matrix_wrapper<Scalar_T>
  { return eigen_matrix_wrapper<Scalar_T>(-m_mat); }

  /**
   * @brief Transpose
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  t() const -> eigen_matrix_wrapper<Scalar_T>
  { return eigen_matrix_wrapper<Scalar_T>(m_mat.transpose()); }

  /**
   * @brief Output to stream
   * @details
   */
  template< typename Scalar_T >
  auto
  operator<< (std::ostream& os, const eigen_matrix_wrapper<Scalar_T>& m) -> std::ostream&
  { return os << m.m_mat; }

  // New Member Implementations
  // ========================

  /**
   * @brief Kronecker matrix product
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  kron(const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
  {
    using namespace Eigen;
    return eigen_matrix_wrapper<Scalar_T>(kroneckerProduct(m_mat, other.m_mat).eval());
  }

  /**
   * @brief Mixed Kronecker matrix product: Dense x Sparse -> Dense (wrapper)
   * @details
   */
  template< typename Scalar_T >
  template< typename Other_Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  kron(const eigen_sparse_wrapper<Other_Scalar_T>& other) const -> eigen_matrix_wrapper<Other_Scalar_T>
  {
    eigen_matrix_wrapper<Other_Scalar_T> result(nbr_rows() * other.nbr_rows(), nbr_cols() * other.nbr_cols());
    result.zeros();

    for (matrix_index_t i = 0; i < nbr_rows(); ++i)
    {
      for (matrix_index_t j = 0; j < nbr_cols(); ++j)
      {
        auto val = static_cast<Other_Scalar_T>((*this)(i, j));
        if (val != Other_Scalar_T(0))
        {
          matrix_index_t r_offset = i * other.nbr_rows();
          matrix_index_t c_offset = j * other.nbr_cols();
          for (auto it = other.begin(); it != other.end(); ++it)
            result(r_offset + it.row(), c_offset + it.col()) = val * static_cast<Other_Scalar_T>(*it);
        }
      }
    }
    return result;
  }

  /**
   * @brief Trace
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  trace() const
  { return m_mat.trace(); }

  /**
   * @brief Infinity norm
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  norm_inf() const
  { return m_mat.cwiseAbs().rowwise().sum().maxCoeff(); }

  /**
   * @brief Squared Frobenius norm
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  norm_frob2() const
  { return m_mat.squaredNorm(); }

  /**
   * @brief Number of non-zeros
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  nnz() const
  { return (m_mat.array() != 0).count(); }

  /**
   * @brief Is NaN?
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  isnan() const -> bool
  { return m_mat.hasNaN(); }

  /**
   * @brief Is infinite?
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_matrix_wrapper<Scalar_T>::
  isinf() const -> bool
  { return !m_mat.allFinite() && !m_mat.hasNaN(); }

  /**
   * @brief Eigenvalues
   * @details
   */
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

  /**
   * @brief Inner product
   * @details
   */
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

  /**
   * @brief Normalization of rotation K
   * @details
   */
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

  // =========================================================================
  // eigen_sparse_wrapper member definitions
  // =========================================================================

  /**
   * @brief Number of rows
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  nbr_rows() const -> matrix_index_t
  { return static_cast<matrix_index_t>(m_mat.rows()); }

  /**
   * @brief Number of columns
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  nbr_cols() const -> matrix_index_t
  { return static_cast<matrix_index_t>(m_mat.cols()); }

  /**
   * @brief Constructor from Eigen Sparse Matrix (e.g. expression result)
   * @details
   */
  template< typename Scalar_T >
  eigen_sparse_wrapper<Scalar_T>::
  eigen_sparse_wrapper(const MatrixType& m)
  : m_mat(m)
  { }

  /**
   * @brief Armadillo/uBLAS/Generator style constructor support
   * @details
   */
  template< typename Scalar_T >
  eigen_sparse_wrapper<Scalar_T>::
  eigen_sparse_wrapper(matrix_index_t rows, matrix_index_t cols, matrix_index_t estimated_nnz)
  {
    set_size(rows, cols);
    if (estimated_nnz > 0) m_mat.reserve(estimated_nnz);
  }

  /**
   * @brief Copy/Move similar to dense
   * @details
   */
  template< typename Scalar_T >
  eigen_sparse_wrapper<Scalar_T>::
  eigen_sparse_wrapper(const eigen_sparse_wrapper<Scalar_T>& other)
  : m_mat(other.m_mat)
  { }

  /**
   * @brief Move constructor
   * @details
   */
  template< typename Scalar_T >
  eigen_sparse_wrapper<Scalar_T>::
  eigen_sparse_wrapper(eigen_sparse_wrapper<Scalar_T>&& other) noexcept
  : m_mat(std::move(other.m_mat))
  { }

  /**
   * @brief Copy assignment
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  operator= (const eigen_sparse_wrapper<Scalar_T>& other) -> eigen_sparse_wrapper<Scalar_T>&
  {
    if (this != &other)
      m_mat = other.m_mat;
    return *this;
  }

  /**
   * @brief Move assignment
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  operator= (eigen_sparse_wrapper<Scalar_T>&& other) noexcept -> eigen_sparse_wrapper<Scalar_T>&
  {
    if (this != &other)
      m_mat = std::move(other.m_mat);
    return *this;
  }

  /**
   * @brief Set size
   * @details
   */
  template< typename Scalar_T >
  void
  eigen_sparse_wrapper<Scalar_T>::
  set_size(matrix_index_t rows, matrix_index_t cols)
  { m_mat.resize(rows, cols); }

  /**
   * @brief Make writable
   * @details
   */
  template< typename Scalar_T >
  void
  eigen_sparse_wrapper<Scalar_T>::
  resize(matrix_index_t rows, matrix_index_t cols, bool preserve)
  {
    // preserve not directly supported in simple resize
    m_mat.resize(rows, cols);
  }

  /**
   * @brief Clear
   * @details
   */
  template< typename Scalar_T >
  void
  eigen_sparse_wrapper<Scalar_T>::
  clear()
  { m_mat.setZero(); }

  /**
   * @brief Set to zero
   * @details
   */
  template< typename Scalar_T >
  void
  eigen_sparse_wrapper<Scalar_T>::
  zeros()
  { m_mat.setZero(); }

  /**
   * @brief Set to zero
   * @details
   */
  template< typename Scalar_T >
  void
  eigen_sparse_wrapper<Scalar_T>::
  zeros(matrix_index_t rows, matrix_index_t cols)
  {
    set_size(rows, cols);
    zeros();
  }

  /**
   * @brief Begin iterator
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  begin() const -> const_iterator
  { return const_iterator(&m_mat, true); }

  /**
   * @brief End iterator
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  end() const -> const_iterator
  { return const_iterator(&m_mat, false); }

  /**
   * @brief Const element access
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  operator() (matrix_index_t i, matrix_index_t j) const -> Scalar_T
  { return m_mat.coeff(i, j); }

  /**
   * @brief Element access
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&
  { return m_mat.coeffRef(i, j); }

  /**
   * @brief Add and assign
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  operator+= (const eigen_sparse_wrapper<Scalar_T>& other) -> eigen_sparse_wrapper<Scalar_T>&
  {
    m_mat += other.m_mat;
    return *this;
  }

  /**
   * @brief Subtract and assign
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  operator-= (const eigen_sparse_wrapper<Scalar_T>& other) -> eigen_sparse_wrapper<Scalar_T>&
  {
    m_mat -= other.m_mat;
    return *this;
  }

  /**
   * @brief Multiply by sparse wrapper
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  operator* (const eigen_sparse_wrapper<Scalar_T>& other) const -> eigen_sparse_wrapper<Scalar_T>
  { return eigen_sparse_wrapper(m_mat * other.m_mat); }

  /**
   * @brief Multiply by scalar and assign
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  operator*= (const Scalar_T& val) -> eigen_sparse_wrapper<Scalar_T>&
  {
    m_mat *= val;
    return *this;
  }

  /**
   * @brief Output to stream
   * @details
   */
  template< typename Scalar_T >
  auto
  operator<< (std::ostream& os, const eigen_sparse_wrapper<Scalar_T>& m) -> std::ostream&
  { return os << m.m_mat; }

  // const_iterator implementation
  /**
   * @brief Iterator support
   * @details
   */
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

  /**
   * @brief Constructor for begin()
   * @details
   */
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

  /**
   * @brief Check if end
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  const_iterator::
  is_end() const -> bool
  { return m_outer >= mp_mat->outerSize(); }

  /**
   * @brief Prefix increment
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  const_iterator::
  operator++ () -> const_iterator&
  {
    advance();
    return *this;
  }

  /**
   * @brief Inequality comparison
   * @details
   */
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

  /**
   * @brief Row index
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  const_iterator::
  row() const -> matrix_index_t
  { return m_inner.row(); }

  /**
   * @brief Column index
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  const_iterator::
  col() const -> matrix_index_t
  { return m_inner.col(); }

  /**
   * @brief Dereference
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  const_iterator::
  operator* () const -> Scalar_T
  { return m_inner.value(); }

  /**
   * @brief Is infinite?
   * @details
   */
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

  /**
   * @brief Is NaN?
   * @details
   */
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

  /**
   * @brief Trace
   * @details
   */
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

  /**
   * @brief Mixed Kronecker matrix product: Sparse x Dense -> Dense (wrapper)
   * @details
   */
  template< typename Scalar_T >
  template< typename Other_Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  kron(const eigen_matrix_wrapper<Other_Scalar_T>& other) const -> eigen_matrix_wrapper<Other_Scalar_T>
  {
    eigen_matrix_wrapper<Other_Scalar_T> result(nbr_rows() * other.nbr_rows(), nbr_cols() * other.nbr_cols());
    result.zeros();

    for (auto it = begin(); it != end(); ++it)
    {
      auto val = static_cast<Other_Scalar_T>(*it);
      matrix_index_t r_offset = it.row() * other.nbr_rows();
      matrix_index_t c_offset = it.col() * other.nbr_cols();

      for (matrix_index_t i = 0; i < other.nbr_rows(); ++i)
        for (matrix_index_t j = 0; j < other.nbr_cols(); ++j)
           result(r_offset + i, c_offset + j) = val * static_cast<Other_Scalar_T>(other(i, j));
    }
    return result;
  }

  /**
   * @brief Kronecker matrix product of sparse wrappers
   * @details
   */
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

  /**
   * @brief Infinity norm
   * @details
   */
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

  /**
   * @brief Squared Frobenius norm
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  norm_frob2() const
  { return m_mat.squaredNorm(); }

  /**
   * @brief Number of non-zeros
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  nnz() const
  { return m_mat.nonZeros(); }

  /**
   * @brief Inner product
   * @details
   */
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

  /**
   * @brief Left Kronecker quotient
   * @details
   */
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

  /**
   * @brief Eigenvalues
   * @details
   */
  template< typename Scalar_T >
  auto
  eigen_sparse_wrapper<Scalar_T>::
  eigenvalues() const -> std::vector<std::complex<double>>
  { throw std::runtime_error("Not implemented for sparse"); } // Usually not computed directly on sparse

  // =========================================================================
  // unit_helper Specializations

  // =========================================================================

  // Specialization for Eigen Wrapper
  template< typename Scalar_T >
  struct unit_helper< eigen_matrix_wrapper<Scalar_T> >
  {
    static auto apply(matrix_index_t dim) -> const eigen_matrix_wrapper<Scalar_T>
    {
      eigen_matrix_wrapper<Scalar_T> result(dim, dim);
      result.m_mat.setIdentity();
      return result;
    }
  };

  // Specialization for Eigen Sparse Wrapper
  template< typename Scalar_T >
  struct unit_helper< eigen_sparse_wrapper<Scalar_T> >
  {
    static auto apply(matrix_index_t dim) -> const eigen_sparse_wrapper<Scalar_T>
    {
      eigen_sparse_wrapper<Scalar_T> result(dim, dim);
      result.m_mat.setIdentity();
      return result;
    }
  };

} }

#endif // _GLUCAT_MATRIX_EIGEN_IMP_H
