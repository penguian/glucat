#ifndef _GLUCAT_MATRIX_ARMA_IMP_H
#define _GLUCAT_MATRIX_ARMA_IMP_H
/***************************************************************************
 G luCat : Generic library of universal Clifford algebra templates         *
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

#include "glucat/errors.h"
#include "glucat/scalar.h"
#include "glucat/matrix.h"
#include "glucat/matrix_arma.h"

#include <set>
#include <vector>
#include <complex>
#include <type_traits>
#include <cstdio>

#include <iostream>
#include <algorithm>

namespace glucat { namespace matrix
{
  /// Solve for arma_matrix_wrapper
  template< typename Scalar_T >
  auto
  solve(arma_matrix_wrapper<Scalar_T>& X, const arma_matrix_wrapper<Scalar_T>& lhs, const arma_matrix_wrapper<Scalar_T>& rhs, int opts = 0) -> bool
  {
    if (lhs.nbr_rows() != lhs.nbr_cols())
      return false;

    return arma::solve(X.m_mat, lhs.m_mat, rhs.m_mat, arma::solve_opts::no_approx);
  }

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

  /// Constructor from eigen_matrix_wrapper
  template< typename Scalar_T >
  template< typename Other_Scalar_T >
  arma_matrix_wrapper<Scalar_T>::
  arma_matrix_wrapper(const eigen_matrix_wrapper<Other_Scalar_T>& other)
  {
    set_size(other.nbr_rows(), other.nbr_cols());
    for (matrix_index_t i = 0; i < nbr_rows(); ++i)
      for (matrix_index_t j = 0; j < nbr_cols(); ++j)
        (*this)(i, j) = static_cast<Scalar_T>(other(i, j));
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
  {  }

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
    arma_matrix_wrapper<Scalar_T> result;
    result.m_mat = arma::kron(m_mat, other.m_mat);
    return result;
  }

  /// Mixed Kronecker matrix product: Dense x Sparse -> Dense (wrapper)
  template< typename Scalar_T >
  template< typename Other_Scalar_T >
  auto
  arma_matrix_wrapper<Scalar_T>::
  kron(const arma_sparse_wrapper<Other_Scalar_T>& other) const -> arma_matrix_wrapper<Other_Scalar_T>
  {
    arma_matrix_wrapper<Other_Scalar_T> result(nbr_rows() * other.nbr_rows(), nbr_cols() * other.nbr_cols());
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
           {
             result(r_offset + it.row(), c_offset + it.col()) = val * static_cast<Other_Scalar_T>(*it);
           }
         }
      }
    }
    return result;
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

  /// Left Kronecker quotient
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
       sum += static_cast<Result_Scalar_T>(*it) * static_cast<Result_Scalar_T>(other(it.row(), it.col()));
    if (nbr_rows() == 0)
      return Result_Scalar_T(0);
    return sum / Result_Scalar_T(static_cast<double>(nbr_rows()));
  }
  /// Mixed Kronecker matrix product: Sparse x Dense -> Dense (wrapper)
  template< typename Scalar_T >
  template< typename Other_Scalar_T >
  auto
  arma_sparse_wrapper<Scalar_T>::
  kron(const arma_matrix_wrapper<Other_Scalar_T>& other) const -> arma_matrix_wrapper<Other_Scalar_T>
  {
    arma_matrix_wrapper<Other_Scalar_T> result(nbr_rows() * other.nbr_rows(), nbr_cols() * other.nbr_cols());
    result.zeros();

    for (auto it = begin(); it != end(); ++it)
    {
       auto val = static_cast<Other_Scalar_T>(*it);
       matrix_index_t r_offset = it.row() * other.nbr_rows();
       matrix_index_t c_offset = it.col() * other.nbr_cols();

       for (matrix_index_t i = 0; i < other.nbr_rows(); ++i)
       {
         for (matrix_index_t j = 0; j < other.nbr_cols(); ++j)
         {
           result(r_offset + i, c_offset + j) = val * static_cast<Other_Scalar_T>(other(i, j));
         }
       }
    }
    return result;
  }

  /// Left Kronecker quotient
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
} }
#endif // _GLUCAT_USE_ARMADILLO
#endif // _GLUCAT_MATRIX_IMP_H
