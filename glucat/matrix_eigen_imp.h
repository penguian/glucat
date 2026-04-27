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
  inline auto
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
  eigen_matrix_wrapper(eigen_matrix_wrapper<Scalar_T>&& other) noexcept
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
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  operator= (const eigen_matrix_wrapper<Scalar_T>& other) -> eigen_matrix_wrapper<Scalar_T>&
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
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  operator= (eigen_matrix_wrapper<Scalar_T>&& other) noexcept -> eigen_matrix_wrapper<Scalar_T>&
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
  inline auto
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
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  nbr_rows() const -> matrix_index_t
  { return static_cast<matrix_index_t>(m_mat.rows()); }

  /*
   * @brief Number of columns
   * @details
   * @tparam Scalar_T
   * @return Number of cols
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  nbr_cols() const -> matrix_index_t
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
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  is_finite() const -> bool
  { return m_mat.allFinite(); }

  /*
   * @brief Has NaN?
   * @details
   * @tparam Scalar_T
   * @return True if has nan
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  has_nan() const -> bool
  { return m_mat.hasNaN(); }

  /*
   * @brief Element access
   * @details
   * @tparam Scalar_T
   * @return Element
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&
  { return m_mat(i, j); }

  /*
   * @brief Const element access
   * @details
   * @tparam Scalar_T
   * @return Element
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  operator() (matrix_index_t i, matrix_index_t j) const -> const Scalar_T&
  { return m_mat(i, j); }

  /*
   * @brief Add and assign
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  operator+= (const eigen_matrix_wrapper<Scalar_T>& other) -> eigen_matrix_wrapper<Scalar_T>&
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
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  operator-= (const eigen_matrix_wrapper<Scalar_T>& other) -> eigen_matrix_wrapper<Scalar_T>&
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
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  operator*= (const Scalar_T& val) -> eigen_matrix_wrapper<Scalar_T>&
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
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  operator/= (const Scalar_T& val) -> eigen_matrix_wrapper<Scalar_T>&
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
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  operator+ (const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
  { return eigen_matrix_wrapper<Scalar_T>(m_mat + other.m_mat); }

  /*
   * @brief Subtraction
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Difference
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  operator- (const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
  { return eigen_matrix_wrapper<Scalar_T>(m_mat - other.m_mat); }

  /*
   * @brief Matrix Multiplication
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Product
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  operator* (const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
  { return eigen_matrix_wrapper<Scalar_T>(m_mat * other.m_mat); }

  /*
   * @brief Unary -
   * @details
   * @tparam Scalar_T
   * @return Unary minus
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  operator- () const -> eigen_matrix_wrapper<Scalar_T>
  { return eigen_matrix_wrapper<Scalar_T>(-m_mat); }

  /*
   * @brief Transpose
   * @details
   * @tparam Scalar_T
   * @return Transpose
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  t() const -> eigen_matrix_wrapper<Scalar_T>
  { return eigen_matrix_wrapper<Scalar_T>(m_mat.transpose()); }

  /*
   * @brief Output to stream
   * @details
   * @tparam Scalar_T
   * @param os Output stream
   * @param m Matrix
   * @return Output stream
   */
  template< typename Scalar_T >
  inline auto
  operator<< (std::ostream& os, const eigen_matrix_wrapper<Scalar_T>& m) -> std::ostream&
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
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  trace() const -> Scalar_T
  { return m_mat.trace(); }

  /*
   * @brief Eigenvalues
   * @details
   * @tparam Scalar_T
   * @return Eigenvalues
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  eigenvalues() const -> std::vector<std::complex<double>>
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

  /*
   * @brief Infinity norm
   * @details
   * @tparam Scalar_T
   * @return Inf norm
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  norm_inf() const -> typename Eigen::NumTraits<Scalar_T>::Real
  { return m_mat.cwiseAbs().rowwise().sum().maxCoeff(); }

  /*
   * @brief Squared Frobenius norm
   * @details
   * @tparam Scalar_T
   * @return Frob2 norm
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  norm_frob2() const -> typename Eigen::NumTraits<Scalar_T>::Real
  { return m_mat.squaredNorm(); }

  /*
   * @brief Is NaN?
   * @details
   * @tparam Scalar_T
   * @return True if successful or condition met
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  isnan() const -> bool
  { return m_mat.hasNaN(); }

  /*
   * @brief Is infinite?
   * @details
   * @tparam Scalar_T
   * @return True if successful or condition met
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  isinf() const -> bool
  { return !m_mat.allFinite() && !m_mat.hasNaN(); }

  /*
   * @brief Number of non-zeros
   * @details
   * @tparam Scalar_T
   * @return Number of non-zeros
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  nnz() const -> matrix_index_t
  { return (m_mat.array() != 0).count(); }

  /*
   * @brief Inner product
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Inner product
   */
  template< typename Scalar_T >
  template< typename Result_Scalar_T, typename Other >
  inline auto
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

  /*
   * @brief Kronecker matrix product
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Kronecker product
   */
  template< typename Scalar_T >
  inline auto
  eigen_matrix_wrapper<Scalar_T>::
  kron(const eigen_matrix_wrapper<Scalar_T>& other) const -> eigen_matrix_wrapper<Scalar_T>
  {
    using namespace Eigen;
    return eigen_matrix_wrapper<Scalar_T>(kroneckerProduct(m_mat, other.m_mat).eval());
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
  inline auto
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
  inline auto
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

  /*
   * @brief Number of rows
   * @details
   * @tparam Scalar_T
   * @return Number of rows
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  nbr_rows() const -> matrix_index_t
  { return static_cast<matrix_index_t>(m_mat.rows()); }

  /*
   * @brief Number of columns
   * @details
   * @tparam Scalar_T
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  nbr_cols() const -> matrix_index_t
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
  eigen_sparse_wrapper(eigen_sparse_wrapper<Scalar_T>&& other) noexcept
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
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  operator= (const eigen_sparse_wrapper<Scalar_T>& other) -> eigen_sparse_wrapper<Scalar_T>&
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
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  operator= (eigen_sparse_wrapper<Scalar_T>&& other) noexcept -> eigen_sparse_wrapper<Scalar_T>&
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
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  begin() const -> const_iterator
  { return const_iterator(&m_mat, true); }

  /*
   * @brief End iterator
   * @details
   * @tparam Scalar_T
   * @return Iterator
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  end() const -> const_iterator
  { return const_iterator(&m_mat, false); }

  /*
   * @brief Const element access
   * @details
   * @tparam Scalar_T
   * @return Element
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  operator() (matrix_index_t i, matrix_index_t j) const -> Scalar_T
  { return m_mat.coeff(i, j); }

  /*
   * @brief Element access
   * @details
   * @tparam Scalar_T
   * @return Element
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  operator() (matrix_index_t i, matrix_index_t j) -> Scalar_T&
  { return m_mat.coeffRef(i, j); }

  /*
   * @brief Add and assign
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  operator+= (const eigen_sparse_wrapper<Scalar_T>& other) -> eigen_sparse_wrapper<Scalar_T>&
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
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  operator-= (const eigen_sparse_wrapper<Scalar_T>& other) -> eigen_sparse_wrapper<Scalar_T>&
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
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  operator* (const eigen_sparse_wrapper<Scalar_T>& other) const -> eigen_sparse_wrapper<Scalar_T>
  { return eigen_sparse_wrapper(m_mat * other.m_mat); }

  /*
   * @brief Multiply by scalar and assign
   * @details
   * @tparam Scalar_T
   * @param val Value
   * @return Reference to this
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  operator*= (const Scalar_T& val) -> eigen_sparse_wrapper<Scalar_T>&
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
  inline auto
  operator<< (std::ostream& os, const eigen_sparse_wrapper<Scalar_T>& m) -> std::ostream&
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
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  const_iterator::
  is_end() const -> bool
  { return m_outer >= static_cast<matrix_index_t>(mp_mat->outerSize()); }

  /*
   * @brief Prefix increment
   * @details
   * @tparam Scalar_T
   * @return Sum
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  const_iterator::
  operator++ () -> const_iterator&
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
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  const_iterator::
  operator!= (const const_iterator& other) const -> bool
  {
    if (m_outer != other.m_outer)
      return true;
    if (m_outer >= static_cast<matrix_index_t>(mp_mat->outerSize()))
      return false;
    return m_inner != other.m_inner;
  }

  /*
   * @brief Row index
   * @details
   * @tparam Scalar_T
   * @return Result
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  const_iterator::
  row() const -> matrix_index_t
  { return m_inner.row(); }

  /*
   * @brief Column index
   * @details
   * @tparam Scalar_T
   * @return Result
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  const_iterator::
  col() const -> matrix_index_t
  { return m_inner.col(); }

  /*
   * @brief Dereference
   * @details
   * @tparam Scalar_T
   * @return Product
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  const_iterator::
  operator* () const -> Scalar_T
  { return m_inner.value(); }

  /*
   * @brief Is infinite?
   * @details
   * @tparam Scalar_T
   * @return True if successful or condition met
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  isinf() const -> bool
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
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  isnan() const -> bool
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
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  trace() const -> Scalar_T
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
  inline auto
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

  /*
   * @brief Kronecker matrix product of sparse wrappers
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Kronecker product
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  kron(const eigen_sparse_wrapper<Scalar_T>& other) const -> eigen_sparse_wrapper<Scalar_T>
  {
    eigen_sparse_wrapper<Scalar_T> result(nbr_rows() * other.nbr_rows(), nbr_cols() * other.nbr_cols());
    std::vector<Eigen::Triplet<Scalar_T>> triplets;
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
   * @brief Infinity norm
   * @details
   * @tparam Scalar_T
   * @return Inf norm
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  norm_inf() const -> typename Eigen::NumTraits<Scalar_T>::Real
  {
    Eigen::Vector<typename Eigen::NumTraits<Scalar_T>::Real, Eigen::Dynamic> row_sums(nbr_rows());
    row_sums.setZero();
    for (matrix_index_t k = 0; k < static_cast<matrix_index_t>(m_mat.outerSize()); ++k)
      for (typename MatrixType::InnerIterator it(m_mat, k); it; ++it)
        row_sums(it.row()) += numeric_traits<Scalar_T>::abs(it.value());
    return row_sums.maxCoeff();
  }

  /*
   * @brief Squared Frobenius norm
   * @details
   * @tparam Scalar_T
   * @return Frob2 norm
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  norm_frob2() const -> typename Eigen::NumTraits<Scalar_T>::Real
  { return m_mat.squaredNorm(); }

  /*
   * @brief Number of non-zeros
   * @details
   * @tparam Scalar_T
   * @return Number of non-zeros
   */
  template< typename Scalar_T >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  nnz() const -> matrix_index_t
  { return m_mat.nonZeros(); }

  /*
   * @brief Inner product
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   * @return Inner product
   */
  template< typename Scalar_T >
  template< typename Result_Scalar_T, typename Other >
  inline auto
  eigen_sparse_wrapper<Scalar_T>::
  inner(const Other& other) const -> Result_Scalar_T
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
   * @brief Left Kronecker quotient
   * @details
   * @tparam Scalar_T
   * @param rhs Right hand side
   * @param mono Value
   * @return Left Kronecker quotient
   */
  template< typename Scalar_T >
  template< typename RHS_T >
  inline auto
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

  /*
   * @brief Eigenvalues
   * @details
   * @tparam Scalar_T
   * @return Eigenvalues
   */
  template< typename Scalar_T >
  inline auto
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
    static inline auto apply(matrix_index_t dim) -> eigen_matrix_wrapper<Scalar_T>
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
    static inline auto apply(matrix_index_t dim) -> eigen_sparse_wrapper<Scalar_T>
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
  inline auto
  operator* (const eigen_sparse_wrapper<Scalar_T>& m, Scalar_T s) -> eigen_sparse_wrapper<Scalar_T>
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
  inline auto
  operator* (Scalar_T s, const eigen_sparse_wrapper<Scalar_T>& m) -> eigen_sparse_wrapper<Scalar_T>
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
  inline auto
  operator+ (const eigen_sparse_wrapper<Scalar_T>& lhs, const eigen_sparse_wrapper<Scalar_T>& rhs) -> eigen_sparse_wrapper<Scalar_T>
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
  inline auto
  operator- (const eigen_sparse_wrapper<Scalar_T>& lhs, const eigen_sparse_wrapper<Scalar_T>& rhs) -> eigen_sparse_wrapper<Scalar_T>
  {
    eigen_sparse_wrapper<Scalar_T> res(lhs);
    res -= rhs;
    return res;
  }



} }
#ifdef GLUCAT_DOCTEST
#include <doctest.h>
#include <iostream>

TEST_CASE("matrix::eigen_matrix_wrapper<Scalar_T>") {
  using namespace glucat::matrix;
  using Scalar_T = double;
  using Matrix_T = eigen_matrix_wrapper<Scalar_T>;

  SUBCASE("Dense initialization and basic operations") {
    const matrix_index_t rows = 2, cols = 3;
    Matrix_T mat(rows, cols);
    mat.zeros();
    CHECK(mat.nbr_rows() == rows);
    CHECK(mat.nbr_cols() == cols);
    CHECK(mat.nnz() == 0);
    CHECK(mat.trace() == doctest::Approx(0.0));

    mat(0, 0) = 1.0;
    mat(1, 1) = 2.0;
    mat(0, 2) = -3.0;
    CHECK(mat.nnz() == 3);
    CHECK(mat.trace() == doctest::Approx(3.0));
    CHECK(mat.norm_inf() == doctest::Approx(4.0)); // |1| + |0| + |-3| = 4
    CHECK(mat.norm_frob2() == doctest::Approx(1.0*1.0 + 2.0*2.0 + (-3.0)*(-3.0)));
  }

    SUBCASE("Dimensions") {
      Matrix_T mat_dims(2, 2);
      CHECK(mat_dims.nbr_rows() == 2);
      CHECK(mat_dims.nbr_cols() == 2);
    }

    SUBCASE("Assignment and Operations") {
    Matrix_T mat;
    mat.unit(3, 3);
    CHECK(mat.nbr_rows() == 3);
    CHECK(mat.nbr_cols() == 3);
    CHECK(mat.nnz() == 3);
    CHECK(mat.trace() == doctest::Approx(3.0));
  }

  SUBCASE("Sparse initialization and basic operations") {
    using Sparse_T = eigen_sparse_wrapper<Scalar_T>;
    const matrix_index_t rows = 4, cols = 4;
    Sparse_T mat(rows, cols);
    mat.zeros();
    CHECK(mat.nbr_rows() == rows);
    CHECK(mat.nbr_cols() == cols);
    CHECK(mat.nnz() == 0);

    mat(0, 0) = 5.0;
    mat(3, 3) = -2.0;
    mat(1, 2) = 1.0;
    CHECK(mat.nnz() == 3);
    CHECK(mat.trace() == doctest::Approx(3.0));
    CHECK(mat.norm_inf() == doctest::Approx(5.0));
    CHECK(mat.norm_frob2() == doctest::Approx(25.0 + 4.0 + 1.0));
  }
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
