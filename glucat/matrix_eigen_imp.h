#ifndef _GLUCAT_MATRIX_EIGEN_IMP_H
#define _GLUCAT_MATRIX_EIGEN_IMP_H
/**************************************************************************
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
 ***************************************************************************
 ***************************************************************************/

#include "glucat/errors.h"
#include "glucat/scalar.h"
#include "glucat/matrix.h"

#include <set>
#include <vector>
#include <complex>
#include <type_traits>
#include <cstdio>
#include <iostream>
#include <algorithm>

namespace glucat { namespace matrix
{
  // =========================================================================
  // Functions for Wrappers (to mimic Armadillo)
  // =========================================================================

  /*
   * @brief Solve
   * @details
   * @tparam Scalar_T
   * @param X Value
   * @param lhs Left hand side
   * @param rhs Right hand side
   * @param opts Options
   */
  template< typename Scalar_T >
  inline bool
  solve(eigen_matrix_wrapper<Scalar_T>& X, const eigen_matrix_wrapper<Scalar_T>& lhs, const eigen_matrix_wrapper<Scalar_T>& rhs, int opts = 0)
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


  /*
   * @brief Armadillo constructor (rows, cols)
   * @details
   * @tparam Scalar_T
   * @param rows Number of rows
   * @param cols Number of columns
   */
  template< typename Scalar_T >
  inline
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(matrix_index_t rows, matrix_index_t cols)
  {
    set_size(rows, cols);
    m_mat.setZero();
  }

  /*
   * @brief Constructor from Eigen expressions (e.g. m * s)
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   */
  template< typename Scalar_T >
  template< typename Derived_T >
  inline
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(const Eigen::MatrixBase<Derived_T>& other)
  { m_mat = other; }

  /*
   * @brief Generic Interop Constructor (e.g. from Armadillo matrix)
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   */
  template< typename Scalar_T >
  template< typename Other_Matrix_T >
  inline
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

  /*
   * @brief Copy constructor
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   */
  template< typename Scalar_T >
  inline
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(const eigen_matrix_wrapper<Scalar_T>& other)
  : m_mat(other.m_mat)
  { }

  /*
   * @brief Move constructor
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   */
  template< typename Scalar_T >
  inline
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(eigen_matrix_wrapper<Scalar_T>&& other) noexcept(std::is_nothrow_move_constructible_v<Scalar_T>)
  : m_mat(std::move(other.m_mat))
  { }

  /*
   * @brief Copy assignment
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>&
  eigen_matrix_wrapper<Scalar_T>::
  operator= (const eigen_matrix_wrapper<Scalar_T>& other)
  {
    if (this != &other)
      m_mat = other.m_mat;
    return *this;
  }

  /*
   * @brief Move assignment
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   */
  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>&
  eigen_matrix_wrapper<Scalar_T>::
  operator= (eigen_matrix_wrapper<Scalar_T>&& other) noexcept(std::is_nothrow_move_assignable_v<Scalar_T>)
  {
    if (this != &other)
      m_mat = std::move(other.m_mat);
    return *this;
  }

  /*
   * @brief Constructor from eigen_sparse_wrapper
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   */
  template< typename Scalar_T >
  template< typename Other_Scalar_T >
  inline
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(const eigen_sparse_wrapper<Other_Scalar_T>& other)
  {
    set_size(other.nbr_rows(), other.nbr_cols());
    m_mat.setZero();
    for (auto it = other.begin(); it != other.end(); ++it)
      (*this)(it.row(), it.col()) = static_cast<Scalar_T>(*it);
  }

  /*
   * @brief Generic Interop Assignment
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   */
  template< typename Scalar_T >
  template< typename Other_Matrix_T >
  inline eigen_matrix_wrapper<Scalar_T>&
  eigen_matrix_wrapper<Scalar_T>::
  operator= (const Other_Matrix_T& other)
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

  /*
   * @brief Constructor from Eigen
   * @details
   * @tparam Scalar_T
   * @param m Matrix
   */
  template< typename Scalar_T >
  inline
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(const MatrixType& m)
  : m_mat(m)
  { }

  /*
   * @brief Constructor from Eigen (move)
   * @details
   * @tparam Scalar_T
   * @param m Matrix
   */
  template< typename Scalar_T >
  inline
  eigen_matrix_wrapper<Scalar_T>::
  eigen_matrix_wrapper(MatrixType&& m)
  : m_mat(std::move(m))
  { }

  /*
   * @brief Set size
   * @details
   * @tparam Scalar_T
   * @param rows Number of rows
   * @param cols Number of columns
   */
  template< typename Scalar_T >
  inline void
  eigen_matrix_wrapper<Scalar_T>::
  set_size(matrix_index_t rows, matrix_index_t cols)
  { m_mat.resize(rows, cols); }

  /*
   * @brief Resize
   * @details
   * @tparam Scalar_T
   * @param rows Number of rows
   * @param cols Number of columns
   * @param preserve Preserve data?
   */


  /*
   * @brief Number of rows
   * @details
   * @tparam Scalar_T
   * @return Number of rows
   */
  template< typename Scalar_T >
  inline matrix_index_t
  eigen_matrix_wrapper<Scalar_T>::
  nbr_rows() const
  { return static_cast<matrix_index_t>(m_mat.rows()); }

  /*
   * @brief Number of columns
   * @details
   * @tparam Scalar_T
   * @return Number of cols
   */
  template< typename Scalar_T >
  inline matrix_index_t
  eigen_matrix_wrapper<Scalar_T>::
  nbr_cols() const
  { return static_cast<matrix_index_t>(m_mat.cols()); }

  /*
   * @brief Clear
   * @details
   * @tparam Scalar_T
   */
  template< typename Scalar_T >
  inline void
  eigen_matrix_wrapper<Scalar_T>::
  clear()
  { m_mat.setZero(); }

  /*
   * @brief Set to zero
   * @details
   * @tparam Scalar_T
   */
  template< typename Scalar_T >
  inline void
  eigen_matrix_wrapper<Scalar_T>::
  zeros()
  { m_mat.setZero(); }

  /*
   * @brief Set size then set to zero
   * @details
   * @tparam Scalar_T
   * @param rows Number of rows
   * @param cols Number of columns
   */
  template< typename Scalar_T >
  inline void
  eigen_matrix_wrapper<Scalar_T>::
  zeros(matrix_index_t rows, matrix_index_t cols)
  {
    set_size(rows, cols);
    zeros();
  }

  /*
   * @brief Set to identity
   * @details
   * @tparam Scalar_T
   * @param rows Number of rows
   * @param cols Number of columns
   */
  template< typename Scalar_T >
  inline void
  eigen_matrix_wrapper<Scalar_T>::
  unit(matrix_index_t rows, matrix_index_t cols)
  {
    set_size(rows, cols);
    m_mat.setIdentity();
  }

  /*
   * @brief Is finite?
   * @details
   * @tparam Scalar_T
   * @return True if is finite
   */
  template< typename Scalar_T >
  inline bool
  eigen_matrix_wrapper<Scalar_T>::
  is_finite() const
  { return m_mat.allFinite(); }

  /*
   * @brief Has NaN?
   * @details
   * @tparam Scalar_T
   * @return True if has nan
   */
  template< typename Scalar_T >
  inline bool
  eigen_matrix_wrapper<Scalar_T>::
  has_nan() const
  { return m_mat.hasNaN(); }

  /*
   * @brief Element access
   * @details
   * @tparam Scalar_T
   * @return Element
   */
  template< typename Scalar_T >
  inline Scalar_T&
  eigen_matrix_wrapper<Scalar_T>::
  operator() (matrix_index_t i, matrix_index_t j)
  { return m_mat(i, j); }

  /*
   * @brief Const element access
   * @details
   * @tparam Scalar_T
   * @return Element
   */
  template< typename Scalar_T >
  inline const Scalar_T&
  eigen_matrix_wrapper<Scalar_T>::
  operator() (matrix_index_t i, matrix_index_t j) const
  { return m_mat(i, j); }

  /*
   * @brief Add and assign
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>&
  eigen_matrix_wrapper<Scalar_T>::
  operator+= (const eigen_matrix_wrapper<Scalar_T>& other)
  {
    m_mat += other.m_mat;
    return *this;
  }

  /*
   * @brief Subtract and assign
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>&
  eigen_matrix_wrapper<Scalar_T>::
  operator-= (const eigen_matrix_wrapper<Scalar_T>& other)
  {
    m_mat -= other.m_mat;
    return *this;
  }

  /*
   * @brief Multiply by scalar and assign
   * @details
   * @tparam Scalar_T
   * @param val Value
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>&
  eigen_matrix_wrapper<Scalar_T>::
  operator*= (const Scalar_T& val)
  {
    m_mat *= val;
    return *this;
  }

  /*
   * @brief Divide by scalar and assign
   * @details
   * @tparam Scalar_T
   * @param val Value
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>&
  eigen_matrix_wrapper<Scalar_T>::
  operator/= (const Scalar_T& val)
  {
    m_mat /= val;
    return *this;
  }

  /*
   * @brief Addition
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Sum
   */
  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>
  eigen_matrix_wrapper<Scalar_T>::
  operator+ (const eigen_matrix_wrapper<Scalar_T>& other) const
  { return eigen_matrix_wrapper<Scalar_T>(MatrixType(m_mat + other.m_mat)); }

  /*
   * @brief Subtraction
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Difference
   */
  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>
  eigen_matrix_wrapper<Scalar_T>::
  operator- (const eigen_matrix_wrapper<Scalar_T>& other) const
  { return eigen_matrix_wrapper<Scalar_T>(MatrixType(m_mat - other.m_mat)); }

  /*
   * @brief Matrix Multiplication
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Product
   */
  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>
  eigen_matrix_wrapper<Scalar_T>::
  operator* (const eigen_matrix_wrapper<Scalar_T>& other) const
  { return eigen_matrix_wrapper<Scalar_T>(MatrixType(m_mat * other.m_mat)); }

  /*
   * @brief Unary -
   * @details
   * @tparam Scalar_T
   * @return Unary minus
   */
  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>
  eigen_matrix_wrapper<Scalar_T>::
  operator- () const
  { return eigen_matrix_wrapper<Scalar_T>(MatrixType(-m_mat)); }

  /*
   * @brief Transpose
   * @details
   * @tparam Scalar_T
   * @return Transpose
   */
  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>
  eigen_matrix_wrapper<Scalar_T>::
  t() const
  { return eigen_matrix_wrapper<Scalar_T>(MatrixType(m_mat.transpose())); }

  /*
   * @brief Output to stream
   * @details
   * @tparam Scalar_T
   * @param os Output stream
   * @param m Matrix
   * @return Output stream
   */
  template< typename Scalar_T >
  inline std::ostream&
  operator<< (std::ostream& os, const eigen_matrix_wrapper<Scalar_T>& m)
  { return os << m.m_mat; }

  // New Member Implementations
  // ========================

  /*
   * @brief Trace
   * @details
   * @tparam Scalar_T
   * @return Trace
   */
  template< typename Scalar_T >
  inline Scalar_T
  eigen_matrix_wrapper<Scalar_T>::
  trace() const
  { return m_mat.trace(); }

  /*
   * @brief Eigenvalues
   * @details
   * @tparam Scalar_T
   * @return Eigenvalues
   */
  template< typename Scalar_T >
  inline std::vector<std::complex<double>>
  eigen_matrix_wrapper<Scalar_T>::
  eigenvalues() const
  {
    // Optimized for real matrices as per matrix representation
    if constexpr (std::is_arithmetic_v<Scalar_T> || std::is_same_v<Scalar_T, double> || std::is_same_v<Scalar_T, float> || std::is_same_v<Scalar_T, long double>)
    {
       Eigen::EigenSolver<typename eigen_matrix_wrapper<Scalar_T>::MatrixType> es(m_mat);
       const auto& E = es.eigenvalues();
       std::vector<std::complex<double>> result(E.size());
       for (matrix_index_t i = 0; i < static_cast<matrix_index_t>(E.size()); ++i)
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
       for (matrix_index_t i = 0; i < static_cast<matrix_index_t>(E.size()); ++i)
         result[i] = E[i];
       return result;
    }
  }

  template< typename Scalar_T >
  inline bool
  eigen_matrix_wrapper<Scalar_T>::
  operator==(const eigen_matrix_wrapper& other) const
  {
    return m_mat == other.m_mat;
  }

  /*
   * @brief Infinity norm
   * @details
   * @tparam Scalar_T
   * @return Inf norm
   */
  template< typename Scalar_T >
  inline Scalar_T
  eigen_matrix_wrapper<Scalar_T>::
  norm_inf() const
  { return static_cast<Scalar_T>(m_mat.cwiseAbs().rowwise().sum().maxCoeff()); }

  /*
   * @brief Squared Frobenius norm
   * @details
   * @tparam Scalar_T
   * @return Frob2 norm
   */
  template< typename Scalar_T >
  inline Scalar_T
  eigen_matrix_wrapper<Scalar_T>::
  norm_frob2() const
  { return static_cast<Scalar_T>(m_mat.squaredNorm()); }

  /*
   * @brief Is NaN?
   * @details
   * @tparam Scalar_T
   * @return True if successful or condition met
   */
  template< typename Scalar_T >
  inline bool
  eigen_matrix_wrapper<Scalar_T>::
  isnan() const
  { return m_mat.hasNaN(); }

  /*
   * @brief Is infinite?
   * @details
   * @tparam Scalar_T
   * @return True if successful or condition met
   */
  template< typename Scalar_T >
  inline bool
  eigen_matrix_wrapper<Scalar_T>::
  isinf() const
  { return !m_mat.allFinite() && !m_mat.hasNaN(); }

  /*
   * @brief Number of non-zeros
   * @details
   * @tparam Scalar_T
   * @return Number of non-zeros
   */
  template< typename Scalar_T >
  inline matrix_index_t
  eigen_matrix_wrapper<Scalar_T>::
  nnz() const
  { return (m_mat.array() != 0).count(); }

  template< typename Scalar_T >
  inline bool
  eigen_matrix_wrapper<Scalar_T>::
  is_zero() const
  { return (m_mat.array() == 0).all(); }

  /*
   * @brief Inner product
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Inner product
   */
  template< typename Scalar_T >
  template< typename Result_Scalar_T, typename Other >
  inline Result_Scalar_T
  eigen_matrix_wrapper<Scalar_T>::
  inner(const Other& other) const
  {
    Result_Scalar_T sum = Result_Scalar_T(0);
    for (matrix_index_t i = 0; i < nbr_rows(); ++i)
      for (matrix_index_t j = 0; j < nbr_cols(); ++j)
        sum += static_cast<Result_Scalar_T>((*this)(i, j)) * static_cast<Result_Scalar_T>(other(i, j));

    if (nbr_rows() == 0) return Result_Scalar_T(0);
    return sum / Result_Scalar_T(static_cast<double>(nbr_rows()));
  }

  /*
   * @brief Trace of product
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Trace(A*B) / Dim
   */
  template< typename Scalar_T >
  template< typename Other >
  inline Scalar_T
  eigen_matrix_wrapper<Scalar_T>::
  trace_product(const Other& other) const
  {
    using traits_t = numeric_traits<Scalar_T>;
    const auto dim = matrix_index_t(m_mat.rows());
    if (dim == 0) return Scalar_T(0);
    Scalar_T tr = (m_mat.cwiseProduct(other.m_mat.transpose())).sum();
    return tr / traits_t::to_scalar_t(dim);
  }

  template< typename Scalar_T >
  inline void
  eigen_matrix_wrapper<Scalar_T>::
  involute()
  {
    const auto rows = m_mat.rows();
    const auto cols = m_mat.cols();
    for (matrix_index_t r = 0; r < rows; ++r)
      for (matrix_index_t c = 0; c < cols; ++c)
        if (std::popcount((unsigned int)(r ^ c)) % 2 != 0)
          m_mat(r, c) *= -1.0;
  }

  template< typename Scalar_T >
  inline void
  eigen_matrix_wrapper<Scalar_T>::
  reverse(index_t p, index_t q)
  {
    // Reversion: M(x^\dagger) = H M(x)^T H^{-1}
    // For Euclidean algebras (q=0), H = I.
    if (q == 0)
    {
      m_mat.transposeInPlace();
      return;
    }
    
    // For non-Euclidean, use similarity transformation of transpose.
    // For now, implement as a general sign flip on transpose.
    // The sign depends on the signature.
    m_mat.transposeInPlace();
    const auto rows = m_mat.rows();
    const auto cols = m_mat.cols();
    for (matrix_index_t r = 0; r < rows; ++r)
      for (matrix_index_t c = 0; c < cols; ++c)
      {
        int k = std::popcount((unsigned int)(r ^ c));
        if (((k * (k - 1)) / 2) % 2 != 0)
          m_mat(r, c) *= -1.0;
      }
  }

  template< typename Scalar_T >
  inline void
  eigen_matrix_wrapper<Scalar_T>::
  similarity_transform(const std::vector<matrix_index_t>& perm, const std::vector<Scalar_T>& signs)
  {
    const matrix_index_t n = m_mat.rows();
    std::vector<matrix_index_t> inv_perm(n);
    for (matrix_index_t i = 0; i < n; ++i) inv_perm[perm[i]] = i;
    MatrixType tmp = MatrixType::Zero(n, n);
    for (matrix_index_t k = 0; k < n; ++k)
      for (matrix_index_t l = 0; l < n; ++l)
      {
        const Scalar_T val = m_mat(k, l);
        if (val == Scalar_T(0)) continue;
        const matrix_index_t i = inv_perm[k];
        const matrix_index_t j = inv_perm[l];
        tmp(i, j) = signs[i] * val * signs[j];
      }
    m_mat = tmp;
  }

  template< typename Scalar_T >
  inline void
  eigen_matrix_wrapper<Scalar_T>::
  transpose_similarity_transform(const std::vector<matrix_index_t>& perm, const std::vector<Scalar_T>& signs)
  {
    const matrix_index_t n = m_mat.rows();
    std::vector<matrix_index_t> inv_perm(n);
    for (matrix_index_t i = 0; i < n; ++i) inv_perm[perm[i]] = i;
    MatrixType tmp = MatrixType::Zero(n, n);
    for (matrix_index_t k = 0; k < n; ++k)
      for (matrix_index_t l = 0; l < n; ++l)
      {
        const Scalar_T val = m_mat(k, l);
        if (val == Scalar_T(0)) continue;
        const matrix_index_t i = inv_perm[l];
        const matrix_index_t j = inv_perm[k];
        tmp(i, j) = signs[i] * val * signs[j];
      }
    m_mat = tmp;
  }

  /*
   * @brief Kronecker matrix product
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Kronecker product
   */
  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>
  eigen_matrix_wrapper<Scalar_T>::
  kron(const eigen_matrix_wrapper<Scalar_T>& other) const
  {
    using namespace Eigen;
    if (nbr_rows() == 2 && nbr_cols() == 2)
    {
      eigen_matrix_wrapper result(2 * other.nbr_rows(), 2 * other.nbr_cols());
      matrix_index_t br = other.nbr_rows();
      matrix_index_t bc = other.nbr_cols();
      result.m_mat.block(0,  0,  br, bc) = (*this)(0,0) * other.m_mat;
      result.m_mat.block(0,  bc, br, bc) = (*this)(0,1) * other.m_mat;
      result.m_mat.block(br, 0,  br, bc) = (*this)(1,0) * other.m_mat;
      result.m_mat.block(br, bc, br, bc) = (*this)(1,1) * other.m_mat;
      return result;
    }
    return eigen_matrix_wrapper<Scalar_T>(kroneckerProduct(m_mat, other.m_mat).eval());
  }

  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>
  eigen_matrix_wrapper<Scalar_T>::
  mono_kron(const eigen_matrix_wrapper<Scalar_T>& other) const
  {
    return kron(other);
  }

  template< typename Scalar_T >
  inline eigen_matrix_wrapper<Scalar_T>
  eigen_matrix_wrapper<Scalar_T>::
  mono_prod(const eigen_matrix_wrapper<Scalar_T>& other) const
  {
    return eigen_matrix_wrapper<Scalar_T>(MatrixType(m_mat * other.m_mat));
  }

  /*
   * @brief Mixed Kronecker matrix product: Dense x Sparse -> Dense (wrapper)
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Kronecker product
   */
  template< typename Scalar_T >
  template< typename Other_Scalar_T >
  inline eigen_matrix_wrapper<Other_Scalar_T>
  eigen_matrix_wrapper<Scalar_T>::
  kron(const eigen_sparse_wrapper<Other_Scalar_T>& other) const
  {
    if (nbr_rows() == 2 && nbr_cols() == 2)
    {
      eigen_matrix_wrapper<Other_Scalar_T> result(2 * other.nbr_rows(), 2 * other.nbr_cols());
      matrix_index_t br = other.nbr_rows();
      matrix_index_t bc = other.nbr_cols();
      result.m_mat.block(0,  0,  br, bc) = static_cast<Other_Scalar_T>((*this)(0,0)) * other.m_mat;
      result.m_mat.block(0,  bc, br, bc) = static_cast<Other_Scalar_T>((*this)(0,1)) * other.m_mat;
      result.m_mat.block(br, 0,  br, bc) = static_cast<Other_Scalar_T>((*this)(1,0)) * other.m_mat;
      result.m_mat.block(br, bc, br, bc) = static_cast<Other_Scalar_T>((*this)(1,1)) * other.m_mat;
      return result;
    }

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

  /*
   * @brief Left Kronecker quotient
   * @details
   * @tparam Scalar_T
   * @param rhs Right hand side
   * @param mono Value
   * @return Left Kronecker quotient
   */
  template< typename Scalar_T >
  template< typename RHS_T >
  inline RHS_T
  eigen_matrix_wrapper<Scalar_T>::
  nork(const RHS_T& rhs, bool mono) const
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
            if constexpr (requires { result.m_mat; rhs.m_mat; })
            {
              result.m_mat += static_cast<typename RHS_T::value_type>(val) * rhs.m_mat.block(start_row, start_col, blk_rows, blk_cols);
            }
            else
            {
              for (matrix_index_t i = 0; i < blk_rows; ++i)
                for (matrix_index_t j = 0; j < blk_cols; ++j)
                  result(i, j) += static_cast<typename RHS_T::value_type>(val) * static_cast<typename RHS_T::value_type>(rhs(start_row + i, start_col + j));
            }
          }
        }
      }

      // Normalize
      auto norm_sq = numeric_traits<typename RHS_T::value_type>::to_scalar_t(nbr_rows());
      if (norm_sq != numeric_traits<typename RHS_T::value_type>::to_scalar_t(1))
      {
        if constexpr (requires { result.m_mat /= norm_sq; })
          result.m_mat /= norm_sq;
        else if constexpr (requires { result(0, 0); })
          for (matrix_index_t i = 0; i < result.nbr_rows(); ++i)
            for (matrix_index_t j = 0; j < result.nbr_cols(); ++j)
              result(i, j) /= norm_sq;
      }
      return result;
  }

  // =========================================================================
  // eigen_sparse_wrapper member definitions
  // =========================================================================

  /*
   * @brief Number of rows
   * @details
   * @tparam Scalar_T
   * @return Number of rows
   */
  template< typename Scalar_T >
  inline matrix_index_t
  eigen_sparse_wrapper<Scalar_T>::
  nbr_rows() const
  { return static_cast<matrix_index_t>(m_mat.rows()); }

  /*
   * @brief Number of columns
   * @details
   * @tparam Scalar_T
   */
  template< typename Scalar_T >
  inline matrix_index_t
  eigen_sparse_wrapper<Scalar_T>::
  nbr_cols() const
  { return static_cast<matrix_index_t>(m_mat.cols()); }

  /*
   * @brief Constructor from Eigen Sparse Matrix (e.g. expression result)
   * @details
   * @tparam Scalar_T
   * @param m Matrix
   */
  template< typename Scalar_T >
  inline
  eigen_sparse_wrapper<Scalar_T>::
  eigen_sparse_wrapper(const MatrixType& m)
  : m_mat(m)
  { }

  /*
   * @brief Armadillo/uBLAS/Generator style constructor support
   * @details
   * @tparam Scalar_T
   * @param rows Number of rows
   * @param cols Number of columns
   * @param estimated_nnz Estimated number of non-zeros
   */
  template< typename Scalar_T >
  inline
  eigen_sparse_wrapper<Scalar_T>::
  eigen_sparse_wrapper(matrix_index_t rows, matrix_index_t cols, matrix_index_t estimated_nnz)
  {
    set_size(rows, cols);
    if (estimated_nnz > 0) m_mat.reserve(estimated_nnz);
  }

  /*
   * @brief Copy/Move similar to dense
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   */
  template< typename Scalar_T >
  inline
  eigen_sparse_wrapper<Scalar_T>::
  eigen_sparse_wrapper(const eigen_sparse_wrapper<Scalar_T>& other)
  : m_mat(other.m_mat)
  { }

  /*
   * @brief Move constructor
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   */
  template< typename Scalar_T >
  inline
  eigen_sparse_wrapper<Scalar_T>::
  eigen_sparse_wrapper(eigen_sparse_wrapper<Scalar_T>&& other) noexcept(std::is_nothrow_move_constructible_v<Scalar_T>)
  : m_mat(std::move(other.m_mat))
  { }

  /*
   * @brief Copy assignment
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>&
  eigen_sparse_wrapper<Scalar_T>::
  operator= (const eigen_sparse_wrapper<Scalar_T>& other)
  {
    if (this != &other)
      m_mat = other.m_mat;
    return *this;
  }

  /*
   * @brief Move assignment
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>&
  eigen_sparse_wrapper<Scalar_T>::
  operator= (eigen_sparse_wrapper<Scalar_T>&& other) noexcept(std::is_nothrow_move_assignable_v<Scalar_T>)
  {
    if (this != &other)
      m_mat = std::move(other.m_mat);
    return *this;
  }

  /*
   * @brief Set size
   * @details
   * @tparam Scalar_T
   * @param rows Number of rows
   * @param cols Number of columns
   */
  template< typename Scalar_T >
  inline void
  eigen_sparse_wrapper<Scalar_T>::
  set_size(matrix_index_t rows, matrix_index_t cols)
  { m_mat.resize(rows, cols); }

  /*
   * @brief Make writable
   * @details
   * @tparam Scalar_T
   * @param rows Number of rows
   * @param cols Number of columns
   * @param preserve Preserve data?
   */


  /*
   * @brief Clear
   * @details
   * @tparam Scalar_T
   */
  template< typename Scalar_T >
  inline void
  eigen_sparse_wrapper<Scalar_T>::
  clear()
  { m_mat.setZero(); }

  /*
   * @brief Set to zero
   * @details
   * @tparam Scalar_T
   */
  template< typename Scalar_T >
  inline void
  eigen_sparse_wrapper<Scalar_T>::
  zeros()
  { m_mat.setZero(); }

  /*
   * @brief Set to zero
   * @details
   * @tparam Scalar_T
   * @param rows Number of rows
   * @param cols Number of columns
   */
  template< typename Scalar_T >
  inline void
  eigen_sparse_wrapper<Scalar_T>::
  zeros(matrix_index_t rows, matrix_index_t cols)
  {
    set_size(rows, cols);
    zeros();
  }

  /*
   * @brief Set to identity
   * @details
   * @tparam Scalar_T
   * @param rows Number of rows
   * @param cols Number of columns
   */
  template< typename Scalar_T >
  inline void
  eigen_sparse_wrapper<Scalar_T>::
  unit(matrix_index_t rows, matrix_index_t cols)
  {
    set_size(rows, cols);
    m_mat.setIdentity();
  }

  /*
   * @brief Begin iterator
   * @details
   * @tparam Scalar_T
   * @return Iterator
   */
  template< typename Scalar_T >
  inline typename eigen_sparse_wrapper<Scalar_T>::const_iterator
  eigen_sparse_wrapper<Scalar_T>::
  begin() const
  { return const_iterator(&m_mat, true); }

  /*
   * @brief End iterator
   * @details
   * @tparam Scalar_T
   * @return Iterator
   */
  template< typename Scalar_T >
  inline typename eigen_sparse_wrapper<Scalar_T>::const_iterator
  eigen_sparse_wrapper<Scalar_T>::
  end() const
  { return const_iterator(&m_mat, false); }

  /*
   * @brief Const element access
   * @details
   * @tparam Scalar_T
   * @return Element
   */
  template< typename Scalar_T >
  inline Scalar_T
  eigen_sparse_wrapper<Scalar_T>::
  operator() (matrix_index_t i, matrix_index_t j) const
  { return m_mat.coeff(i, j); }

  /*
   * @brief Element access
   * @details
   * @tparam Scalar_T
   * @return Element
   */
  template< typename Scalar_T >
  inline Scalar_T&
  eigen_sparse_wrapper<Scalar_T>::
  operator() (matrix_index_t i, matrix_index_t j)
  { return m_mat.coeffRef(i, j); }

  /*
   * @brief Add and assign
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>&
  eigen_sparse_wrapper<Scalar_T>::
  operator+= (const eigen_sparse_wrapper<Scalar_T>& other)
  {
    m_mat += other.m_mat;
    return *this;
  }

  /*
   * @brief Subtract and assign
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>&
  eigen_sparse_wrapper<Scalar_T>::
  operator-= (const eigen_sparse_wrapper<Scalar_T>& other)
  {
    m_mat -= other.m_mat;
    return *this;
  }

  /*
   * @brief Multiply by sparse wrapper
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Product
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>
  eigen_sparse_wrapper<Scalar_T>::
  operator* (const eigen_sparse_wrapper<Scalar_T>& other) const
  { return eigen_sparse_wrapper(m_mat * other.m_mat); }

  /*
   * @brief Multiply by scalar and assign
   * @details
   * @tparam Scalar_T
   * @param val Value
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>&
  eigen_sparse_wrapper<Scalar_T>::
  operator*= (const Scalar_T& val)
  {
    m_mat *= val;
    return *this;
  }

  /*
   * @brief Output to stream
   * @details
   * @tparam Scalar_T
   * @param os Output stream
   * @param m Matrix
   * @return Output stream
   */
  template< typename Scalar_T >
  inline std::ostream&
  operator<< (std::ostream& os, const eigen_sparse_wrapper<Scalar_T>& m)
  { return os << m.m_mat; }

  // const_iterator implementation
  /*
   * @brief Iterator support
   * @details
   * @tparam Scalar_T
   * @param mat Value
   * @param start Value
   */
  template< typename Scalar_T >
  inline
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
      m_outer = static_cast<matrix_index_t>(mp_mat->outerSize());
  }

  /*
   * @brief Constructor for begin()
   * @details
   * @tparam Scalar_T
   */
  template< typename Scalar_T >
  inline void
  eigen_sparse_wrapper<Scalar_T>::
  const_iterator::
  advance()
  {
    if (m_inner)
      ++m_inner;
    while (!m_inner && m_outer < static_cast<matrix_index_t>(mp_mat->outerSize()))
    {
      m_outer++;
      if (m_outer < static_cast<matrix_index_t>(mp_mat->outerSize()))
        m_inner = InnerIterator(*mp_mat, m_outer);
    }
  }

  /*
   * @brief Check if end
   * @details
   * @tparam Scalar_T
   * @return True if is end
   */
  template< typename Scalar_T >
  inline bool
  eigen_sparse_wrapper<Scalar_T>::const_iterator::
  is_end() const
  { return m_outer >= static_cast<matrix_index_t>(mp_mat->outerSize()); }

  /*
   * @brief Prefix increment
   * @details
   * @tparam Scalar_T
   * @return Sum
   */
  template< typename Scalar_T >
  inline typename eigen_sparse_wrapper<Scalar_T>::const_iterator&
  eigen_sparse_wrapper<Scalar_T>::const_iterator::
  operator++ ()
  {
    advance();
    return *this;
  }

  /*
   * @brief Inequality comparison
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return True if not equal
   */
  template< typename Scalar_T >
  inline bool
  eigen_sparse_wrapper<Scalar_T>::const_iterator::
  operator!= (const const_iterator& other) const
  {
    if (mp_mat != other.mp_mat)
      return true;
    if (m_outer != other.m_outer)
      return true;
    if (m_outer >= static_cast<matrix_index_t>(mp_mat->outerSize()))
      return false;
    return m_inner.index() != other.m_inner.index();
  }

  /*
   * @brief Equality comparison
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return True if equal
   */
  template< typename Scalar_T >
  inline bool
  eigen_sparse_wrapper<Scalar_T>::const_iterator::
  operator== (const const_iterator& other) const
  { return !(*this != other); }

  /*
   * @brief Row index
   * @details
   * @tparam Scalar_T
   * @return Result
   */
  template< typename Scalar_T >
  inline matrix_index_t
  eigen_sparse_wrapper<Scalar_T>::const_iterator::
  row() const
  { return m_inner.row(); }

  /*
   * @brief Column index
   * @details
   * @tparam Scalar_T
   * @return Result
   */
  template< typename Scalar_T >
  inline matrix_index_t
  eigen_sparse_wrapper<Scalar_T>::const_iterator::
  col() const
  { return m_inner.col(); }

  /*
   * @brief Dereference
   * @details
   * @tparam Scalar_T
   * @return Product
   */
  template< typename Scalar_T >
  inline Scalar_T
  eigen_sparse_wrapper<Scalar_T>::const_iterator::
  operator* () const
  { return m_inner.value(); }

  /*
   * @brief Is infinite?
   * @details
   * @tparam Scalar_T
   * @return True if successful or condition met
   */
  template< typename Scalar_T >
  inline bool
  eigen_sparse_wrapper<Scalar_T>::
  isinf() const
  {
    for (matrix_index_t k = 0; k < static_cast<matrix_index_t>(m_mat.outerSize()); ++k)
      for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
        if (numeric_traits<Scalar_T>::isinf(it.value()))
          return true;
    return false;
  }

  /*
   * @brief Is NaN?
   * @details
   * @tparam Scalar_T
   * @return True if successful or condition met
   */
  template< typename Scalar_T >
  inline bool
  eigen_sparse_wrapper<Scalar_T>::
  isnan() const
  {
    for (matrix_index_t k = 0; k < static_cast<matrix_index_t>(m_mat.outerSize()); ++k)
      for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
        if (numeric_traits<Scalar_T>::isnan(it.value()))
          return true;
    return false;
  }

  /*
   * @brief Trace
   * @details
   * @tparam Scalar_T
   * @return Trace
   */
  template< typename Scalar_T >
  inline Scalar_T
  eigen_sparse_wrapper<Scalar_T>::
  trace() const
  { return m_mat.diagonal().sum(); }

  /*
   * @brief Mixed Kronecker matrix product: Sparse x Dense -> Dense (wrapper)
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Kronecker product
   */
  template< typename Scalar_T >
  template< typename Other_Scalar_T >
  inline eigen_matrix_wrapper<Other_Scalar_T>
  eigen_sparse_wrapper<Scalar_T>::
  kron(const eigen_matrix_wrapper<Other_Scalar_T>& other) const
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

  /*
   * @brief Kronecker matrix product of sparse wrappers
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Kronecker product
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>
  eigen_sparse_wrapper<Scalar_T>::
  kron(const eigen_sparse_wrapper<Scalar_T>& other) const
  {
    eigen_sparse_wrapper<Scalar_T> result(nbr_rows() * other.nbr_rows(), nbr_cols() * other.nbr_cols());
    std::vector<Eigen::Triplet<Scalar_T>> triplets;
    triplets.reserve(nnz() * other.nnz());
    // Iterate lhs
    for (matrix_index_t k = 0; k < static_cast<matrix_index_t>(m_mat.outerSize()); ++k)
    {
      for (typename MatrixType::InnerIterator itA(m_mat, k); itA; ++itA)
      {
        auto rA = itA.row();
        auto cA = itA.col();
        auto vA = itA.value();

        // Iterate rhs
        for (matrix_index_t l = 0; l < static_cast<matrix_index_t>(other.m_mat.outerSize()); ++l)
          for (typename MatrixType::InnerIterator itB(other.m_mat, l); itB; ++itB)
            triplets.emplace_back(rA * other.nbr_rows() + itB.row(), cA * other.nbr_cols() + itB.col(), vA * itB.value());
      }
    }
    result.m_mat.setFromTriplets(triplets.begin(), triplets.end());

    return result;
  }

  /*
   * @brief Monomial Kronecker matrix product of sparse wrappers
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Kronecker product
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>
  eigen_sparse_wrapper<Scalar_T>::
  mono_kron(const eigen_sparse_wrapper<Scalar_T>& other) const
  {
    eigen_sparse_wrapper<Scalar_T> result(nbr_rows() * other.nbr_rows(), nbr_cols() * other.nbr_cols());
    std::vector<Eigen::Triplet<Scalar_T>> triplets;
    triplets.reserve(nnz() * other.nnz());
    
    const matrix_index_t other_rows = other.nbr_rows();
    const matrix_index_t other_cols = other.nbr_cols();

    for (matrix_index_t k = 0; k < static_cast<matrix_index_t>(m_mat.outerSize()); ++k)
    {
      typename MatrixType::InnerIterator itA(m_mat, k);
      if (itA)
      {
        auto rA = itA.row();
        auto cA = itA.col();
        auto vA = itA.value();

        for (matrix_index_t l = 0; l < static_cast<matrix_index_t>(other.m_mat.outerSize()); ++l)
        {
          typename MatrixType::InnerIterator itB(other.m_mat, l);
          if (itB)
          {
            triplets.emplace_back(rA * other_rows + itB.row(), cA * other_cols + itB.col(), vA * itB.value());
          }
        }
      }
    }
    result.m_mat.setFromTriplets(triplets.begin(), triplets.end());

    return result;
  }

  /*
   * @brief Monomial product of sparse wrappers
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Product
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>
  eigen_sparse_wrapper<Scalar_T>::
  mono_prod(const eigen_sparse_wrapper<Scalar_T>& other) const
  {
    eigen_sparse_wrapper<Scalar_T> result(nbr_rows(), nbr_cols());
    std::vector<Eigen::Triplet<Scalar_T>> triplets;
    triplets.reserve(nnz());

    for (matrix_index_t cB = 0; cB < static_cast<matrix_index_t>(other.m_mat.outerSize()); ++cB)
    {
      typename MatrixType::InnerIterator itB(other.m_mat, cB);
      if (itB)
      {
        auto rB = itB.row();
        auto vB = itB.value();

        typename MatrixType::InnerIterator itA(m_mat, rB);
        if (itA)
        {
          triplets.emplace_back(itA.row(), cB, itA.value() * vB);
        }
      }
    }
    result.m_mat.setFromTriplets(triplets.begin(), triplets.end());

    return result;
  }

  /*
   * @brief Infinity norm
   * @details
   * @tparam Scalar_T
   * @return Inf norm
   */
  template< typename Scalar_T >
  inline Scalar_T
  eigen_sparse_wrapper<Scalar_T>::
  norm_inf() const
  {
    Eigen::Vector<typename Eigen::NumTraits<Scalar_T>::Real, Eigen::Dynamic> row_sums(nbr_rows());
    row_sums.setZero();
    for (matrix_index_t k = 0; k < static_cast<matrix_index_t>(m_mat.outerSize()); ++k)
      for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
        row_sums(it.row()) += numeric_traits<Scalar_T>::abs(it.value());
    return static_cast<Scalar_T>(row_sums.maxCoeff());
  }

  /*
   * @brief Squared Frobenius norm
   * @details
   * @tparam Scalar_T
   * @return Frob2 norm
   */
  template< typename Scalar_T >
  inline Scalar_T
  eigen_sparse_wrapper<Scalar_T>::
  norm_frob2() const
  { return static_cast<Scalar_T>(m_mat.squaredNorm()); }

  /*
   * @brief Number of non-zeros
   * @details
   * @tparam Scalar_T
   * @return Number of non-zeros
   */
  template< typename Scalar_T >
  inline matrix_index_t
  eigen_sparse_wrapper<Scalar_T>::
  nnz() const
  { return m_mat.nonZeros(); }

  template< typename Scalar_T >
  inline bool
  eigen_sparse_wrapper<Scalar_T>::
  is_zero() const
  {
    for (int k=0; k<m_mat.outerSize(); ++k)
      for (typename MatrixType::InnerIterator it(m_mat,k); it; ++it)
        if (it.value() != Scalar_T(0)) return false;
    return true;
  }
 
  template< typename Scalar_T >
  inline bool
  eigen_sparse_wrapper<Scalar_T>::
  operator==(const eigen_sparse_wrapper& other) const
  {
    // Sparse comparison: ensure we don't return false due to stored zeros
    if (m_mat.rows() != other.m_mat.rows() || m_mat.cols() != other.m_mat.cols()) return false;
    Eigen::SparseMatrix<Scalar_T> diff = m_mat - other.m_mat;
    // Iterative check is often safer than nnz() after subtraction
    for (int k=0; k<diff.outerSize(); ++k)
      for (typename Eigen::SparseMatrix<Scalar_T>::InnerIterator it(diff,k); it; ++it)
        if (it.value() != Scalar_T(0)) return false;
    return true;
  }

  /*
   * @brief Inner product
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Inner product
   */
  template< typename Scalar_T >
  template< typename Result_Scalar_T, typename Other >
  inline Result_Scalar_T
  eigen_sparse_wrapper<Scalar_T>::
  inner(const Other& other) const
  {
     Result_Scalar_T sum = Result_Scalar_T(0);
     for (matrix_index_t k=0; k < static_cast<matrix_index_t>(m_mat.outerSize()); ++k)
       for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
       {
          sum += static_cast<Result_Scalar_T>(it.value()) * static_cast<Result_Scalar_T>(other(it.row(), it.col()));
       }
     if (nbr_rows() == 0) return Result_Scalar_T(0);
     return sum / Result_Scalar_T(static_cast<double>(nbr_rows()));
  }

  /*
   * @brief Trace of product
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Trace(A*B) / Dim
   */
  template< typename Scalar_T >
  template< typename Other >
  inline Scalar_T
  eigen_sparse_wrapper<Scalar_T>::
  trace_product(const Other& other) const
  {
    using traits_t = numeric_traits<Scalar_T>;
    const auto dim = matrix_index_t(m_mat.rows());
    if (dim == 0) return Scalar_T(0);
    // Trace(AB) = sum_i sum_j A_ij * B_ji
    Scalar_T tr = Scalar_T(0);
    for (int k=0; k < m_mat.outerSize(); ++k)
      for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
        tr += it.value() * other(it.col(), it.row());
    return tr / traits_t::to_scalar_t(dim);
  }

  template< typename Scalar_T >
  inline void
  eigen_sparse_wrapper<Scalar_T>::
  involute()
  {
    for (int k=0; k < m_mat.outerSize(); ++k)
      for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
        if (std::popcount((unsigned int)(it.row() ^ it.col())) % 2 != 0)
          it.valueRef() *= -1.0;
  }

  template< typename Scalar_T >
  inline void
  eigen_sparse_wrapper<Scalar_T>::
  reverse(index_t p, index_t q)
  {
    if (q == 0)
    {
      m_mat = m_mat.transpose().eval();
      return;
    }
    m_mat = m_mat.transpose().eval();
    for (int k=0; k < m_mat.outerSize(); ++k)
      for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
      {
        int grade = std::popcount((unsigned int)(it.row() ^ it.col()));
        if (((grade * (grade - 1)) / 2) % 2 != 0)
          it.valueRef() *= -1.0;
      }
  }

  template< typename Scalar_T >
  inline void
  eigen_sparse_wrapper<Scalar_T>::
  similarity_transform(const std::vector<matrix_index_t>& perm, const std::vector<Scalar_T>& signs)
  {
    const matrix_index_t n = m_mat.rows();
    std::vector<matrix_index_t> inv_perm(n);
    for (matrix_index_t i = 0; i < n; ++i) inv_perm[perm[i]] = i;

    std::vector<Eigen::Triplet<Scalar_T>> triplets;
    triplets.reserve(m_mat.nonZeros());

    for (int k=0; k < m_mat.outerSize(); ++k)
      for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
      {
        const matrix_index_t new_r = inv_perm[it.row()];
        const matrix_index_t new_c = inv_perm[it.col()];
        triplets.push_back(Eigen::Triplet<Scalar_T>(new_r, new_c, signs[new_r] * it.value() * signs[new_c]));
      }
    m_mat.setFromTriplets(triplets.begin(), triplets.end());
  }

  template< typename Scalar_T >
  inline void
  eigen_sparse_wrapper<Scalar_T>::
  transpose_similarity_transform(const std::vector<matrix_index_t>& perm, const std::vector<Scalar_T>& signs)
  {
    const matrix_index_t n = m_mat.rows();
    std::vector<matrix_index_t> inv_perm(n);
    for (matrix_index_t i = 0; i < n; ++i) inv_perm[perm[i]] = i;

    std::vector<Eigen::Triplet<Scalar_T>> triplets;
    triplets.reserve(m_mat.nonZeros());

    for (int k=0; k < m_mat.outerSize(); ++k)
      for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
      {
        const matrix_index_t new_r = inv_perm[it.col()];
        const matrix_index_t new_c = inv_perm[it.row()];
        triplets.push_back(Eigen::Triplet<Scalar_T>(new_r, new_c, signs[new_r] * it.value() * signs[new_c]));
      }
    m_mat.setFromTriplets(triplets.begin(), triplets.end());
  }

  /*
   * @brief Left Kronecker quotient
   * @details
   * @tparam Scalar_T
   * @param rhs Right hand side
   * @param mono Value
   * @return Left Kronecker quotient
   */
  template< typename Scalar_T >
  template< typename RHS_T >
  inline RHS_T
  eigen_sparse_wrapper<Scalar_T>::
  nork(const RHS_T& rhs, bool mono) const
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
         if constexpr (requires { result.m_mat; rhs.m_mat; })
         {
           result.m_mat += static_cast<typename RHS_T::value_type>(val) * rhs.m_mat.block(start_row, start_col, blk_rows, blk_cols);
         }
         else
         {
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
       if constexpr (requires { result.m_mat /= norm_sq; })
         result.m_mat /= norm_sq;
       else if constexpr (requires { result(0, 0); })
         for (matrix_index_t i = 0; i < result.nbr_rows(); ++i)
           for (matrix_index_t j = 0; j < result.nbr_cols(); ++j)
             result(i, j) /= norm_sq;
     }
     return result;
  }

  /*
   * @brief Eigenvalues
   * @details
   * @tparam Scalar_T
   * @return Eigenvalues
   */
  template< typename Scalar_T >
  inline std::vector<std::complex<double>>
  eigen_sparse_wrapper<Scalar_T>::
  eigenvalues() const
  { throw std::runtime_error("Not implemented for sparse"); } // Usually not computed directly on sparse

  // =========================================================================
  // unit_helper Specializations

  // =========================================================================

  // Specialization for Eigen Wrapper
  template< typename Scalar_T >
  struct unit_helper< eigen_matrix_wrapper<Scalar_T> >
  {
    static inline eigen_matrix_wrapper<Scalar_T> apply(matrix_index_t dim)
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
    static inline eigen_sparse_wrapper<Scalar_T> apply(matrix_index_t dim)
    {
      eigen_sparse_wrapper<Scalar_T> result(dim, dim);
      result.m_mat.setIdentity();
      return result;
    }
  };

  /*
   * @brief Product of sparse wrapper and scalar
   * @details
   * @tparam Scalar_T
   * @param m Matrix
   * @param s Scalar
   * @return Product
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>
  operator* (const eigen_sparse_wrapper<Scalar_T>& m, Scalar_T s)
  {
    eigen_sparse_wrapper<Scalar_T> res(m);
    res *= s;
    return res;
  }

  /*
   * @brief Product of scalar and sparse wrapper
   * @details
   * @tparam Scalar_T
   * @param s Scalar
   * @param m Matrix
   * @return Product
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>
  operator* (Scalar_T s, const eigen_sparse_wrapper<Scalar_T>& m)
  {
    return m * s;
  }

  /*
   * @brief Sum of sparse wrappers
   * @details
   * @tparam Scalar_T
   * @param lhs Left Hand Side
   * @param rhs Right Hand Side
   * @return Sum
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>
  operator+ (const eigen_sparse_wrapper<Scalar_T>& lhs, const eigen_sparse_wrapper<Scalar_T>& rhs)
  {
    eigen_sparse_wrapper<Scalar_T> res(lhs);
    res += rhs;
    return res;
  }

  /*
   * @brief Difference of sparse wrappers
   * @details
   * @tparam Scalar_T
   * @param lhs Left Hand Side
   * @param rhs Right Hand Side
   * @return Difference
   */
  template< typename Scalar_T >
  inline eigen_sparse_wrapper<Scalar_T>
  operator- (const eigen_sparse_wrapper<Scalar_T>& lhs, const eigen_sparse_wrapper<Scalar_T>& rhs)
  {
    eigen_sparse_wrapper<Scalar_T> res(lhs);
    res -= rhs;
    return res;
  }



} }
#ifdef GLUCAT_DOCTEST
#include <doctest/doctest.h>
#include <sstream>
#include <utility>

namespace glucat { namespace matrix {

  template <typename Scalar_T>
  void test_eigen_wrappers_templated()
  {
    using Matrix_T = eigen_matrix_wrapper<Scalar_T>;
    using Sparse_T = eigen_sparse_wrapper<Scalar_T>;

    SUBCASE("Dense Matrix: Constructors and Assignment") {
      Matrix_T m1(2, 2);
      m1.zeros();
      CHECK(m1.nbr_rows() == 2);
      CHECK(m1.nbr_cols() == 2);
      CHECK(m1.nnz() == 0);
      CHECK(m1.is_zero() == true);

      // Test nnz unreliability: (m - m).nnz() might be > 0 but is_zero() should be true
      Sparse_T sm1(2,2);
      sm1(0,0) = Scalar_T(1);
      Sparse_T sm2 = sm1 - sm1;
      // Depending on Eigen version/flags, sm2.nonZeros() might be > 0
      // but is_zero() and operator== must handle it correctly.
      CHECK(sm2.is_zero() == true);
      CHECK(sm1 == sm1);
      // Compare sm2 (sparse zero) with m1 (dense zero)
      // Wait! operator== between sparse and dense?
      // m1 is Dense, sm2 is Sparse. We should check sm2.is_zero()
      CHECK(sm2.is_zero() == true);

      m1(0, 0) = Scalar_T(1);
      m1(1, 1) = Scalar_T(2);

      // Copy constructor
      Matrix_T m2(m1);
      CHECK(m2(0, 0) == Scalar_T(1));
      CHECK(m2(1, 1) == Scalar_T(2));

      // Move constructor
      Matrix_T m3(std::move(m2));
      CHECK(m3(0, 0) == Scalar_T(1));

      // Assignment
      Matrix_T m4;
      m4 = m3;
      CHECK(m4(1, 1) == Scalar_T(2));

      // Move assignment
      Matrix_T m5;
      m5 = std::move(m4);
      CHECK(m5(0, 0) == Scalar_T(1));
    }

    SUBCASE("Dense Matrix: Arithmetic Operators") {
      Matrix_T a(2, 2), b(2, 2);
      a.unit(2, 2);
      b.unit(2, 2);

      auto c = a + b;
      CHECK(c(0, 0) == Scalar_T(2));

      auto d = a - b;
      CHECK(d(0, 0) == Scalar_T(0));

      auto e = a * Scalar_T(3);
      CHECK(e(0, 0) == Scalar_T(3));

      auto f = -a;
      CHECK(f(0, 0) == Scalar_T(-1));

      Matrix_T g(2, 2);
      g(0, 0) = Scalar_T(1); g(0, 1) = Scalar_T(2);
      g(1, 0) = Scalar_T(3); g(1, 1) = Scalar_T(4);
      auto h = g.t();
      CHECK(h(0, 1) == Scalar_T(3));
      CHECK(h(1, 0) == Scalar_T(2));

      // Dense-Dense product
      Matrix_T i = g * h;
      CHECK(i(0, 0) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(5))));

      // Compound assignment
      g += a;
      CHECK(g(0, 0) == Scalar_T(2));
      g -= a;
      CHECK(g(0, 0) == Scalar_T(1));
      g *= Scalar_T(2);
      CHECK(g(0, 0) == Scalar_T(2));

      // unit_helper::apply
      auto u = matrix::unit<Matrix_T>(3);
      CHECK(u.nbr_rows() == 3);
      CHECK(u.trace() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(3))));

      a.zeros(3, 3);
      CHECK(a.nbr_rows() == 3);
      a.clear();
      CHECK(a.is_finite());
    }

    SUBCASE("Dense Matrix: Analysis and Solve") {
      Matrix_T m(2, 2);
      m(0, 0) = Scalar_T(2); m(0, 1) = Scalar_T(1);
      m(1, 0) = Scalar_T(1); m(1, 1) = Scalar_T(2);

      CHECK(m.is_finite());
      CHECK_FALSE(m.isnan());
      CHECK_FALSE(m.isinf());

      Matrix_T rhs(2, 1);
      rhs(0, 0) = Scalar_T(3);
      rhs(1, 0) = Scalar_T(3);

      Matrix_T x(2, 1);
      bool success = solve(x, m, rhs);
      CHECK(success);
      CHECK(x(0, 0) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1))));
      CHECK(x(1, 0) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1))));

      auto ev = m.eigenvalues();
      CHECK(ev.size() == 2);
    }

    SUBCASE("Sparse Matrix: Basic Operations") {
      Sparse_T s(4, 4);
      s.zeros();
      s(0, 0) = Scalar_T(5);
      s(3, 3) = Scalar_T(-2);

      CHECK(s.nnz() == 2);
      CHECK(s.trace() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(3))));

      Sparse_T s2;
      s2.unit(4, 4);
      CHECK(s2.nnz() == 4);
      CHECK(s2.trace() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(4))));

      // Sparse iterator operator!=
      Sparse_T s3(2, 2);
      s3(0, 0) = Scalar_T(1);
      s3(1, 0) = Scalar_T(2); // Two elements in first column (if ColMajor)
      auto it1 = s3.begin();
      auto it2 = s3.begin();
      ++it2;
      CHECK(it1 != it2);
    }

    SUBCASE("Kronecker Product and Nork") {
      Matrix_T a(2, 2), b(2, 2);
      a.unit(2, 2);
      b.unit(2, 2);

      // Dense x Dense
      auto c = a.kron(b);
      CHECK(c.nbr_rows() == 4);
      CHECK(c.trace() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(4))));

      // Nork (Dense)
      auto q = a.nork(c, true);
      CHECK(q.nbr_rows() == 2);
      CHECK(q.trace() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(2))));

      // Mixed: Dense x Sparse
      Sparse_T s(2, 2);
      s.unit(2, 2);
      auto mixed = a.kron(s);
      CHECK(mixed.nbr_rows() == 4);

      // Mixed: Sparse x Dense
      auto mixed2 = s.kron(b);
      CHECK(mixed2.nbr_rows() == 4);

      // Sparse x Sparse
      auto ss = s.kron(s);
      CHECK(ss.nbr_rows() == 4);
      CHECK(ss.nnz() == 4);

      // Nork (Sparse)
      auto qs = s.nork(ss, true);
      CHECK(qs.nbr_rows() == 2);
    }

    SUBCASE("Norms and Inner Product") {
      Matrix_T m(2, 2);
      m(0,0) = Scalar_T(1); m(0,1) = Scalar_T(-2);
      m(1,0) = Scalar_T(3); m(1,1) = Scalar_T(4);

      CHECK(m.norm_inf() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(7)))); // max(1+2, 3+4)
      CHECK(m.norm_frob2() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1+4+9+16))));

      Sparse_T s(2, 2);
      s(0,0) = Scalar_T(1); s(1,1) = Scalar_T(4);
      CHECK(s.norm_inf() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(4))));
      CHECK(s.norm_frob2() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1+16))));

      // Inner product
      auto in = m.template inner<Scalar_T>(m);
      CHECK(in == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1+4+9+16)/Scalar_T(2))));

      auto ins = s.template inner<Scalar_T>(s);
      CHECK(ins == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1+16)/Scalar_T(2))));
    }

    SUBCASE("Interop and Conversion") {
      Sparse_T s(2, 2);
      s.unit(2, 2);

      // Sparse to Dense
      Matrix_T d;
      d = s; // Assignment
      CHECK(d.trace() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(2))));

      Matrix_T d2(s); // Constructor
      CHECK(d2.trace() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(2))));

      // Dense to Sparse? (Need to check if constructor exists)
      // eigen_sparse_wrapper doesn't have a direct constructor from dense wrapper in the code I saw,
      // but it might have via Eigen interop.

      // operator<< (mostly for coverage)
      std::ostringstream oss;
      oss << s;
      CHECK(!oss.str().empty());
    }

    SUBCASE("Iterator and Edge Cases") {
      Sparse_T empty(0, 0);
      CHECK(empty.begin() == empty.end());

      Sparse_T s(2, 2);
      s(1, 1) = Scalar_T(5);
      auto it = s.begin();
      CHECK(it != s.end());
      CHECK(it.row() == 1);
      CHECK(it.col() == 1);
      CHECK(*it == Scalar_T(5));
      ++it;
      CHECK(it == s.end());

      // Cross-matrix comparison
      Sparse_T s2(2, 2);
      s2(1, 1) = Scalar_T(5);
      CHECK(s.begin() != s2.begin());
      CHECK_FALSE(s.begin() == s2.begin());
    }

    SUBCASE("Transforms") {
      // Test similarity_transform and transpose_similarity_transform
      // Permutation: (0 1) -> (1 0), Signs: (1, -1)
      std::vector<matrix_index_t> perm = {1, 0};
      std::vector<Scalar_T> signs = {Scalar_T(1), Scalar_T(-1)};

      // Dense
      Matrix_T d(2, 2);
      d(0, 0) = Scalar_T(1); d(0, 1) = Scalar_T(2);
      d(1, 0) = Scalar_T(3); d(1, 1) = Scalar_T(4);

      Matrix_T d_trans = d;
      d_trans.similarity_transform(perm, signs);
      // Expected: signs[i] * val(perm[i], perm[j]) * signs[j]
      // (0,0): signs[0]*val(1,1)*signs[0] = 1*4*1 = 4
      // (0,1): signs[0]*val(1,0)*signs[1] = 1*3*-1 = -3
      // (1,0): signs[1]*val(0,1)*signs[0] = -1*2*1 = -2
      // (1,1): signs[1]*val(0,0)*signs[1] = -1*1*-1 = 1
      CHECK(d_trans(0, 0) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(4))));
      CHECK(d_trans(0, 1) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(-3))));
      CHECK(d_trans(1, 0) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(-2))));
      CHECK(d_trans(1, 1) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1))));

      // Sparse
      Sparse_T s(2, 2);
      s(0, 0) = Scalar_T(1); s(0, 1) = Scalar_T(2);
      s(1, 0) = Scalar_T(3); s(1, 1) = Scalar_T(4);

      Sparse_T s_trans = s;
      s_trans.similarity_transform(perm, signs);
      CHECK(s_trans(0, 0) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(4))));
      CHECK(s_trans(0, 1) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(-3))));
      CHECK(s_trans(1, 0) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(-2))));
      CHECK(s_trans(1, 1) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1))));

      // Transpose Transform (Dense)
      Matrix_T d_t_trans = d;
      d_t_trans.transpose_similarity_transform(perm, signs);
      // Expected: signs[i] * val(perm[j], perm[i]) * signs[j]
      // (0,0): signs[0]*val(p[0],p[0])*signs[0] = 1*val(1,1)*1 = 4
      // (0,1): signs[0]*val(p[1],p[0])*signs[1] = 1*val(0,1)*-1 = -2
      // (1,0): signs[1]*val(p[0],p[1])*signs[0] = -1*val(1,0)*1 = -3
      // (1,1): signs[1]*val(p[1],p[1])*signs[1] = -1*val(0,0)*-1 = 1
      CHECK(d_t_trans(0, 0) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(4))));
      CHECK(d_t_trans(0, 1) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(-2))));
      CHECK(d_t_trans(1, 0) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(-3))));
      CHECK(d_t_trans(1, 1) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1))));

      // Transpose Transform (Sparse)
      Sparse_T s_t_trans = s;
      s_t_trans.transpose_similarity_transform(perm, signs);
      CHECK(s_t_trans(0, 0) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(4))));
      CHECK(s_t_trans(0, 1) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(-2))));
      CHECK(s_t_trans(1, 0) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(-3))));
      CHECK(s_t_trans(1, 1) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1))));
    }
  }

} } // glucat::matrix

TEST_CASE("matrix::eigen_wrappers") {
  SUBCASE("float")       { glucat::matrix::test_eigen_wrappers_templated<float>(); }
  SUBCASE("double")      { glucat::matrix::test_eigen_wrappers_templated<double>(); }
  SUBCASE("long double") { glucat::matrix::test_eigen_wrappers_templated<long double>(); }
#ifdef _GLUCAT_USE_QD
  SUBCASE("dd_real")     { glucat::matrix::test_eigen_wrappers_templated<dd_real>(); }
  SUBCASE("qd_real")     { glucat::matrix::test_eigen_wrappers_templated<qd_real>(); }
#endif
}
#endif

  // =========================================================================
  // Eigen Internal Cast Specializations for High Precision
  // =========================================================================


#if defined(_GLUCAT_USE_QD) && defined(EIGEN_MAJOR_VERSION)
namespace Eigen { namespace internal {

  // Resolves ambiguity for qd_real(long) and qd_real(unsigned long)
  // by explicitly routing through double for cast<T>(index).
  template<>
  struct cast_impl<long, qd_real> {
    static inline qd_real run(const long& x) {
      return qd_real(static_cast<double>(x));
    }
  };

  template<>
  struct cast_impl<unsigned long, qd_real> {
    static inline qd_real run(const unsigned long& x) {
      return qd_real(static_cast<double>(x));
    }
  };

  // Resolves ambiguity for dd_real(long) and dd_real(unsigned long)
  template<>
  struct cast_impl<long, dd_real> {
    static inline dd_real run(const long& x) {
      return dd_real(static_cast<double>(x));
    }
  };

  template<>
  struct cast_impl<unsigned long, dd_real> {
    static inline dd_real run(const unsigned long& x) {
      return dd_real(static_cast<double>(x));
    }
  };

} } // namespace Eigen::internal
#endif

#endif  // _GLUCAT_MATRIX_EIGEN_IMP_H
