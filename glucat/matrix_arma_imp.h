#ifndef _GLUCAT_MATRIX_ARMA_IMP_H
#define _GLUCAT_MATRIX_ARMA_IMP_H
/**************************************************************************
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
***************************************************************************
***************************************************************************/

#if defined(_GLUCAT_USE_ARMADILLO)

#include <algorithm>
#include <complex>
#include <cstdio>
#include <iostream>
#include <set>
#include <type_traits>
#include <vector>

#include "glucat/errors.h"
#include "glucat/matrix.h"
#include "glucat/matrix_arma.h"
#include "glucat/scalar.h"

namespace glucat
{
  namespace matrix
  {
    /*
     * @brief Solve for arma_matrix_wrapper
     * @details
     * @tparam Scalar_T
     * @param X Value
     * @param lhs Left hand side
     * @param rhs Right hand side
     * @param opts Options
     */
    template <typename Scalar_T>
    inline bool solve(arma_matrix_wrapper<Scalar_T>& X, const arma_matrix_wrapper<Scalar_T>& lhs,
                      const arma_matrix_wrapper<Scalar_T>& rhs, int opts = 0)
    {
      if (lhs.nbr_rows() != lhs.nbr_cols())
        return false;

      return arma::solve(X.m_mat, lhs.m_mat, rhs.m_mat, arma::solve_opts::no_approx);
    }

    // =========================================================================
    // arma_matrix_wrapper Member Definitions
    // =========================================================================

    /*
     * @brief Constructor with size
     * @details
     * @tparam Scalar_T
     * @param rows Number of rows
     * @param cols Number of columns
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>::arma_matrix_wrapper(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.zeros();
    }

    /*
     * @brief Constructor from other matrix type
     * @details
     * @tparam Scalar_T
     * @tparam Other_Matrix_T
     * @param other Other matrix
     */
    template <typename Scalar_T>
    template <typename Other_Matrix_T>
    inline arma_matrix_wrapper<Scalar_T>::arma_matrix_wrapper(const Other_Matrix_T& other)
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

    /*
     * @brief Constructor from eigen_matrix_wrapper
     * @details
     * @tparam Scalar_T
     * @tparam Other_Scalar_T
     * @param other Other matrix
     */
    // Constructor from eigen_matrix_wrapper
    template <typename Scalar_T>
    template <typename Other_Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>::arma_matrix_wrapper(const eigen_matrix_wrapper<Other_Scalar_T>& other)
    {
      set_size(other.nbr_rows(), other.nbr_cols());
      for (matrix_index_t i = 0; i < other.nbr_rows(); ++i)
        for (matrix_index_t j = 0; j < other.nbr_cols(); ++j)
          m_mat(i, j) = static_cast<Scalar_T>(other(i, j));
    }

    // Constructor from eigen_sparse_wrapper
    template <typename Scalar_T>
    template <typename Other_Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>::arma_matrix_wrapper(const eigen_sparse_wrapper<Other_Scalar_T>& other)
    {
      set_size(other.nbr_rows(), other.nbr_cols());
      m_mat.zeros();
      for (auto it = other.begin(); it != other.end(); ++it)
        m_mat(it.row(), it.col()) = static_cast<Scalar_T>(*it);
    }

    /*
     * @brief Constructor from arma_sparse_wrapper
     * @details
     * @tparam Scalar_T
     * @tparam Other_Scalar_T
     * @param other Other matrix
     */
    template <typename Scalar_T>
    template <typename Other_Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>::arma_matrix_wrapper(const arma_sparse_wrapper<Other_Scalar_T>& other)
    {
      // Use efficient Armadillo conversion if possible
      m_mat = arma::conv_to<MatrixType>::from(other.m_mat);
    }

    /*
     * @brief Copy constructor
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>::arma_matrix_wrapper(const arma_matrix_wrapper<Scalar_T>& other)
        : m_mat(other.m_mat)
    {
    }

    /*
     * @brief Move constructor
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>::arma_matrix_wrapper(arma_matrix_wrapper<Scalar_T>&& other) noexcept(
        std::is_nothrow_move_constructible_v<Scalar_T>)
        : m_mat(std::move(other.m_mat))
    {
    }

    /*
     * @brief Copy assignment
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Reference to this
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>& arma_matrix_wrapper<Scalar_T>::operator=(const arma_matrix_wrapper<Scalar_T>& other)
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
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>& arma_matrix_wrapper<Scalar_T>::operator=(arma_matrix_wrapper<Scalar_T>&& other) noexcept(
        std::is_nothrow_move_assignable_v<Scalar_T>)
    {
      if (this != &other)
        m_mat = std::move(other.m_mat);
      return *this;
    }

    /*
     * @brief Assignment from sparse wrapper
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Reference to this
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>& arma_matrix_wrapper<Scalar_T>::operator=(const arma_sparse_wrapper<Scalar_T>& other)
    {
      m_mat = other.m_mat;
      return *this;
    }

    /*
     * @brief Conversion to const MatrixType reference
     * @details
     * @tparam Scalar_T
     * @return Result
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>::operator const MatrixType&() const
    { return m_mat; }

    /*
     * @brief Conversion to MatrixType reference
     * @details
     * @tparam Scalar_T
     * @return Result
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>::operator MatrixType&()
    { return m_mat; }

    /*
     * @brief Set size
     * @details
     * @tparam Scalar_T
     * @param rows Number of rows
     * @param cols Number of columns
     */
    template <typename Scalar_T>
    inline void arma_matrix_wrapper<Scalar_T>::set_size(matrix_index_t rows, matrix_index_t cols)
    {
      m_mat.set_size(rows, cols);
    }

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
    template <typename Scalar_T>
    inline matrix_index_t arma_matrix_wrapper<Scalar_T>::nbr_rows() const
    { return m_mat.n_rows; }

    /*
     * @brief Number of columns
     * @details
     * @tparam Scalar_T
     * @return Number of cols
     */
    template <typename Scalar_T>
    inline matrix_index_t arma_matrix_wrapper<Scalar_T>::nbr_cols() const
    { return m_mat.n_cols; }

    /*
     * @brief Clear
     * @details
     * @tparam Scalar_T
     */
    template <typename Scalar_T>
    inline void arma_matrix_wrapper<Scalar_T>::clear()
    {
      m_mat.zeros();
    }

    /*
     * @brief Set to zero
     * @details
     * @tparam Scalar_T
     * @param rows Number of rows
     * @param cols Number of columns
     */
    template <typename Scalar_T>
    inline void arma_matrix_wrapper<Scalar_T>::zeros(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      zeros();
    }

    /*
     * @brief Set size then set to zero
     * @details
     * @tparam Scalar_T
     */
    template <typename Scalar_T>
    inline void arma_matrix_wrapper<Scalar_T>::zeros()
    {
      m_mat.zeros();
    }

    /*
     * @brief Set to identity
     * @details
     * @tparam Scalar_T
     * @param rows Number of rows
     * @param cols Number of columns
     */
    template <typename Scalar_T>
    inline void arma_matrix_wrapper<Scalar_T>::unit(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.eye();
    }

    /*
     * @brief Element access
     * @details
     * @tparam Scalar_T
     * @return Element
     */
    template <typename Scalar_T>
    inline Scalar_T& arma_matrix_wrapper<Scalar_T>::operator()(matrix_index_t i, matrix_index_t j)
    { return m_mat(i, j); }

    /*
     * @brief Const element access
     * @details
     * @tparam Scalar_T
     * @return Element
     */
    template <typename Scalar_T>
    inline const Scalar_T& arma_matrix_wrapper<Scalar_T>::operator()(matrix_index_t i, matrix_index_t j) const
    { return m_mat(i, j); }

    /*
     * @brief Add and assign
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Reference to this
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>& arma_matrix_wrapper<Scalar_T>::operator+=(const arma_matrix_wrapper<Scalar_T>& other)
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
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>& arma_matrix_wrapper<Scalar_T>::operator-=(const arma_matrix_wrapper<Scalar_T>& other)
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
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>& arma_matrix_wrapper<Scalar_T>::operator*=(const Scalar_T& val)
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
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>& arma_matrix_wrapper<Scalar_T>::operator/=(const Scalar_T& val)
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
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T> arma_matrix_wrapper<Scalar_T>::operator+(const arma_matrix_wrapper<Scalar_T>& other) const
    { return arma_matrix_wrapper(MatrixType(m_mat + other.m_mat)); }

    /*
     * @brief Subtraction
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Difference
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T> arma_matrix_wrapper<Scalar_T>::operator-(const arma_matrix_wrapper<Scalar_T>& other) const
    { return arma_matrix_wrapper(MatrixType(m_mat - other.m_mat)); }

    /*
     * @brief Multiplication
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Product
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T> arma_matrix_wrapper<Scalar_T>::operator*(const arma_matrix_wrapper<Scalar_T>& other) const
    {
      // Force evaluation to MatrixType (arma::Mat) to avoid resolving to generic template constructor with Glue
      MatrixType res_arma = m_mat * other.m_mat;
      return arma_matrix_wrapper(std::move(res_arma));
    }

    /*
     * @brief Unary negation
     * @details
     * @tparam Scalar_T
     * @return Unary minus
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T> arma_matrix_wrapper<Scalar_T>::operator-() const
    { return arma_matrix_wrapper(MatrixType(-m_mat)); }

    /*
     * @brief Transpose
     * @details
     * @tparam Scalar_T
     * @return Transpose
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T> arma_matrix_wrapper<Scalar_T>::t() const
    { return arma_matrix_wrapper(MatrixType(m_mat.t())); }

    /*
     * @brief Begin iterator
     * @details
     * @tparam Scalar_T
     * @return Iterator
     */
    template <typename Scalar_T>
    inline auto arma_matrix_wrapper<Scalar_T>::begin()
    { return m_mat.begin(); }

    /*
     * @brief End iterator
     * @details
     * @tparam Scalar_T
     * @return Iterator
     */
    template <typename Scalar_T>
    inline auto arma_matrix_wrapper<Scalar_T>::end()
    { return m_mat.end(); }

    /*
     * @brief Begin const iterator
     * @details
     * @tparam Scalar_T
     * @return Iterator
     */
    template <typename Scalar_T>
    inline auto arma_matrix_wrapper<Scalar_T>::begin() const
    { return m_mat.begin(); }

    /*
     * @brief End const iterator
     * @details
     * @tparam Scalar_T
     * @return Iterator
     */
    template <typename Scalar_T>
    inline auto arma_matrix_wrapper<Scalar_T>::end() const
    { return m_mat.end(); }

    /*
     * @brief Has infinity?
     * @details
     * @tparam Scalar_T
     * @return True if has inf
     */
    template <typename Scalar_T>
    inline bool arma_matrix_wrapper<Scalar_T>::has_inf() const
    { return m_mat.has_inf(); }

    /*
     * @brief Has NaN?
     * @details
     * @tparam Scalar_T
     * @return True if has nan
     */
    template <typename Scalar_T>
    inline bool arma_matrix_wrapper<Scalar_T>::has_nan() const
    { return m_mat.has_nan(); }

    /*
     * @brief Is finite?
     * @details
     * @tparam Scalar_T
     * @return True if is finite
     */
    template <typename Scalar_T>
    inline bool arma_matrix_wrapper<Scalar_T>::is_finite() const
    { return m_mat.is_finite(); }

    /*
     * @brief Output to stream
     * @details
     * @tparam Scalar_T
     * @param os Output stream
     * @param m Matrix
     * @return Output stream
     */
    template <typename Scalar_T>
    inline std::ostream& operator<<(std::ostream& os, const arma_matrix_wrapper<Scalar_T>& m)
    { return os << m.m_mat; }

    // New Member Implementations (moved from free functions)
    // ====================================================

    /*
     * @brief Kronecker matrix product
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Kronecker product
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T> arma_matrix_wrapper<Scalar_T>::kron(const arma_matrix_wrapper<Scalar_T>& other) const
    {
      arma_matrix_wrapper<Scalar_T> result;
      result.m_mat = arma::kron(m_mat, other.m_mat);
      return result;
    }

    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T> arma_matrix_wrapper<Scalar_T>::mono_kron(const arma_matrix_wrapper<Scalar_T>& other) const
    { return kron(other); }

    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T> arma_matrix_wrapper<Scalar_T>::mono_prod(const arma_matrix_wrapper<Scalar_T>& other) const
    {
      arma_matrix_wrapper<Scalar_T> result;
      result.m_mat = m_mat * other.m_mat;
      return result;
    }

    /*
     * @brief Mixed Kronecker matrix product: Dense x Sparse -> Dense (wrapper)
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Kronecker product
     */
    template <typename Scalar_T>
    template <typename Other_Scalar_T>
    inline arma_matrix_wrapper<Other_Scalar_T> arma_matrix_wrapper<Scalar_T>::kron(
        const arma_sparse_wrapper<Other_Scalar_T>& other) const
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
            result.m_mat.submat(i * other.nbr_rows(), j * other.nbr_cols(), (i + 1) * other.nbr_rows() - 1,
                                (j + 1) * other.nbr_cols() - 1) = val * other.m_mat;
          }
        }
      }
      return result;
    }

    /*
     * @brief Trace
     * @details
     * @tparam Scalar_T
     * @return Trace
     */
    template <typename Scalar_T>
    inline Scalar_T arma_matrix_wrapper<Scalar_T>::trace() const
    { return arma::trace(m_mat); }

    /*
     * @brief Eigenvalues
     * @details
     * @tparam Scalar_T
     * @return Eigenvalues
     */
    template <typename Scalar_T>
    inline std::vector<std::complex<double>> arma_matrix_wrapper<Scalar_T>::eigenvalues() const
    {
      using pod_type = typename arma::get_pod_type<Scalar_T>::result;
      arma::Col<std::complex<pod_type>> eigval;
      arma::eig_gen(eigval, m_mat);
      std::vector<std::complex<double>> result(eigval.n_elem);
      for (matrix_index_t i = 0; i < eigval.n_elem; ++i)
        result[i] = std::complex<double>(eigval[i].real(), eigval[i].imag());

      return result;
    }

    /*
     * @brief Is NaN?
     * @details
     * @tparam Scalar_T
     * @return True if successful or condition met
     */
    template <typename Scalar_T>
    inline bool arma_matrix_wrapper<Scalar_T>::isnan() const
    { return m_mat.has_nan(); }

    /*
     * @brief Is infinite?
     * @details
     * @tparam Scalar_T
     * @return True if successful or condition met
     */
    template <typename Scalar_T>
    inline bool arma_matrix_wrapper<Scalar_T>::isinf() const
    { return m_mat.has_inf(); }

    /*
     * @brief Infinity norm
     * @details
     * @tparam Scalar_T
     * @return Inf norm
     */
    template <typename Scalar_T>
    inline Scalar_T arma_matrix_wrapper<Scalar_T>::norm_inf() const
    { return static_cast<Scalar_T>(arma::norm(m_mat, "inf")); }

    /*
     * @brief Squared Frobenius norm
     * @details
     * @tparam Scalar_T
     * @return Frob2 norm
     */
    template <typename Scalar_T>
    inline Scalar_T arma_matrix_wrapper<Scalar_T>::norm_frob2() const
    { return static_cast<Scalar_T>(arma::accu(arma::square(m_mat))); }

    /*
     * @brief Number of non-zeros
     * @details
     * @tparam Scalar_T
     */
    template <typename Scalar_T>
    inline matrix_index_t arma_matrix_wrapper<Scalar_T>::nnz() const
    { return arma::accu(m_mat != 0); }

    /*
     * @brief Helper to construct from raw arma mat
     * @details
     * @tparam Scalar_T
     * @param m Matrix
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>::arma_matrix_wrapper(const MatrixType& m)
        : m_mat(m)
    {
    }

    /*
     * @brief Helper to construct from raw arma mat
     * @details
     * @tparam Scalar_T
     * @param m Matrix
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T>::arma_matrix_wrapper(MatrixType&& m)
        : m_mat(std::move(m))
    {
    }

    /*
     * @brief Trace of product
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Trace(A*B) / Dim
     */
    template <typename Scalar_T>
    template <typename Other>
    inline Scalar_T arma_matrix_wrapper<Scalar_T>::trace_product(const Other& other) const
    {
      using traits_t = numeric_traits<Scalar_T>;
      const auto dim = matrix_index_t(m_mat.n_rows);
      if (dim == 0)
        return Scalar_T(0);
      // Trace(AB) = accu(A % B.t())
      Scalar_T tr = arma::accu(m_mat % other.m_mat.t());
      return tr / traits_t::to_scalar_t(dim);
    }

    /*
     * @brief Inner product
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Inner product
     */
    template <typename Scalar_T>
    template <typename Result_Scalar_T, typename Other>
    inline Result_Scalar_T arma_matrix_wrapper<Scalar_T>::inner(const Other& other) const
    {
      Result_Scalar_T sum = Result_Scalar_T(0);
      for (matrix_index_t i = 0; i < nbr_rows(); ++i)
        for (matrix_index_t j = 0; j < nbr_cols(); ++j)
          sum += static_cast<Result_Scalar_T>(m_mat(i, j)) * static_cast<Result_Scalar_T>(other(i, j));

      if (nbr_rows() == 0)
        return Result_Scalar_T(0);
      return sum / Result_Scalar_T(static_cast<double>(nbr_rows()));
    }

    template <typename Scalar_T>
    inline void arma_matrix_wrapper<Scalar_T>::similarity_transform(const std::vector<matrix_index_t>& perm,
                                                                    const std::vector<Scalar_T>& signs)
    {
      const matrix_index_t n = nbr_rows();
      std::vector<matrix_index_t> inv_perm(n);
      for (matrix_index_t i = 0; i < n; ++i)
        inv_perm[perm[i]] = i;
      MatrixType tmp = arma::zeros<MatrixType>(n, n);
      for (matrix_index_t k = 0; k < n; ++k)
        for (matrix_index_t l = 0; l < n; ++l)
        {
          const Scalar_T val = m_mat(k, l);
          if (val == Scalar_T(0))
            continue;
          const matrix_index_t i = inv_perm[k];
          const matrix_index_t j = inv_perm[l];
          tmp(i, j) = signs[i] * val * signs[j];
        }
      m_mat = std::move(tmp);
    }

    template <typename Scalar_T>
    inline void arma_matrix_wrapper<Scalar_T>::transpose_similarity_transform(const std::vector<matrix_index_t>& perm,
                                                                              const std::vector<Scalar_T>& signs)
    {
      const matrix_index_t n = nbr_rows();
      std::vector<matrix_index_t> inv_perm(n);
      for (matrix_index_t i = 0; i < n; ++i)
        inv_perm[perm[i]] = i;
      MatrixType tmp = arma::zeros<MatrixType>(n, n);
      for (matrix_index_t k = 0; k < n; ++k)
        for (matrix_index_t l = 0; l < n; ++l)
        {
          const Scalar_T val = m_mat(k, l);
          if (val == Scalar_T(0))
            continue;
          const matrix_index_t i = inv_perm[l];
          const matrix_index_t j = inv_perm[k];
          tmp(i, j) = signs[i] * val * signs[j];
        }
      m_mat = std::move(tmp);
    }

    template <typename Scalar_T>
    inline void arma_matrix_wrapper<Scalar_T>::involute()
    {
      for (matrix_index_t r = 0; r < m_mat.n_rows; ++r)
        for (matrix_index_t c = 0; c < m_mat.n_cols; ++c)
        {
          const Scalar_T val = m_mat(r, c);
          if (val == Scalar_T(0))
            continue;
          if (std::popcount((unsigned int)(r ^ c)) % 2 != 0)
            m_mat(r, c) = -val;
        }
    }

    template <typename Scalar_T>
    inline void arma_matrix_wrapper<Scalar_T>::reverse(index_t p, index_t q)
    {
      m_mat = m_mat.t();
    }

    /*
     * @brief Left Kronecker quotient
     * @details
     * @tparam Scalar_T
     * @param rhs Right hand side
     * @param mono Value
     * @return Left Kronecker quotient
     */
    template <typename Scalar_T>
    template <typename RHS_T>
    inline RHS_T arma_matrix_wrapper<Scalar_T>::nork(const RHS_T& rhs, bool mono) const
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
            result.m_mat += static_cast<typename RHS_T::value_type>(val)
                            * rhs.m_mat.submat(r * blk_rows, c * blk_cols, (r + 1) * blk_rows - 1, (c + 1) * blk_cols - 1);
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

    /*
     * @brief Number of rows
     * @details
     * @tparam Scalar_T
     * @return Number of rows
     */
    template <typename Scalar_T>
    inline matrix_index_t arma_sparse_wrapper<Scalar_T>::nbr_rows() const
    { return m_mat.n_rows; }

    /*
     * @brief Number of columns
     * @details
     * @tparam Scalar_T
     */
    template <typename Scalar_T>
    inline matrix_index_t arma_sparse_wrapper<Scalar_T>::nbr_cols() const
    { return m_mat.n_cols; }

    /*
     * @brief Constructor from Armadillo SpMat
     * @details
     * @tparam Scalar_T
     * @param m Matrix
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T>::arma_sparse_wrapper(const MatrixType& m)
        : m_mat(m)
    {
    }

    /*
     * @brief Constructor with size
     * @details
     * @tparam Scalar_T
     * @param rows Number of rows
     * @param cols Number of columns
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T>::arma_sparse_wrapper(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
    }

    /*
     * @brief Copy constructor
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T>::arma_sparse_wrapper(const arma_sparse_wrapper<Scalar_T>& other)
        : m_mat(other.m_mat)
    {
    }

    /*
     * @brief Move constructor
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T>::arma_sparse_wrapper(arma_sparse_wrapper<Scalar_T>&& other) noexcept(
        std::is_nothrow_move_constructible_v<Scalar_T>)
        : m_mat(std::move(other.m_mat))
    {
    }

    /*
     * @brief Copy assignment
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Reference to this
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T>& arma_sparse_wrapper<Scalar_T>::operator=(const arma_sparse_wrapper<Scalar_T>& other)
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
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T>& arma_sparse_wrapper<Scalar_T>::operator=(arma_sparse_wrapper<Scalar_T>&& other) noexcept(
        std::is_nothrow_move_assignable_v<Scalar_T>)
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
    template <typename Scalar_T>
    inline void arma_sparse_wrapper<Scalar_T>::set_size(matrix_index_t rows, matrix_index_t cols)
    {
      m_mat.set_size(rows, cols);
    }

    /*
     * @brief Resize
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
    template <typename Scalar_T>
    inline void arma_sparse_wrapper<Scalar_T>::clear()
    {
      m_mat.zeros();
    }

    /*
     * @brief Set to zero
     * @details
     * @tparam Scalar_T
     * @param rows Number of rows
     * @param cols Number of columns
     */
    template <typename Scalar_T>
    inline void arma_sparse_wrapper<Scalar_T>::zeros(matrix_index_t rows, matrix_index_t cols)
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
    template <typename Scalar_T>
    inline void arma_sparse_wrapper<Scalar_T>::unit(matrix_index_t rows, matrix_index_t cols)
    {
      set_size(rows, cols);
      m_mat.eye(rows, cols);
    }

    /*
     * @brief Set to zero
     * @details
     * @tparam Scalar_T
     */
    template <typename Scalar_T>
    inline void arma_sparse_wrapper<Scalar_T>::zeros()
    {
      m_mat.zeros();
    }

    /*
     * @brief Begin iterator
     * @details
     * @tparam Scalar_T
     * @return Iterator
     */
    template <typename Scalar_T>
    inline typename arma_sparse_wrapper<Scalar_T>::const_iterator arma_sparse_wrapper<Scalar_T>::begin() const
    { return m_mat.begin(); }

    /*
     * @brief End iterator
     * @details
     * @tparam Scalar_T
     * @return Iterator
     */
    template <typename Scalar_T>
    inline typename arma_sparse_wrapper<Scalar_T>::const_iterator arma_sparse_wrapper<Scalar_T>::end() const
    { return m_mat.end(); }

    /*
     * @brief Const element access
     * @details
     * @tparam Scalar_T
     * @return Element
     */
    template <typename Scalar_T>
    inline Scalar_T arma_sparse_wrapper<Scalar_T>::operator()(matrix_index_t i, matrix_index_t j) const
    { return m_mat(i, j); }

    /*
     * @brief Element access
     * @details
     * @tparam Scalar_T
     * @return Element
     */
    template <typename Scalar_T>
    inline decltype(auto) arma_sparse_wrapper<Scalar_T>::operator()(matrix_index_t i, matrix_index_t j)
    { return m_mat(i, j); }

    /*
     * @brief Add and assign
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Reference to this
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T>& arma_sparse_wrapper<Scalar_T>::operator+=(const arma_sparse_wrapper<Scalar_T>& other)
    {
      m_mat += other.m_mat;
      return *this;
    }

    /*
     * @brief Multiply by sparse wrapper
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Product
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T> arma_sparse_wrapper<Scalar_T>::operator*(const arma_sparse_wrapper<Scalar_T>& other) const
    { return arma_sparse_wrapper(m_mat * other.m_mat); }

    /*
     * @brief Multiply by scalar and assign
     * @details
     * @tparam Scalar_T
     * @param val Value
     * @return Reference to this
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T>& arma_sparse_wrapper<Scalar_T>::operator*=(const Scalar_T& val)
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
    template <typename Scalar_T>
    inline std::ostream& operator<<(std::ostream& os, const arma_sparse_wrapper<Scalar_T>& m)
    { return os << m.m_mat; }

    // Armadillo Sparse Helper Member Functions
    /*
     * @brief Is infinite?
     * @details
     * @tparam Scalar_T
     * @return True if successful or condition met
     */
    template <typename Scalar_T>
    inline bool arma_sparse_wrapper<Scalar_T>::isinf() const
    { return m_mat.has_inf(); }

    /*
     * @brief Is NaN?
     * @details
     * @tparam Scalar_T
     * @return True if successful or condition met
     */
    template <typename Scalar_T>
    inline bool arma_sparse_wrapper<Scalar_T>::isnan() const
    { return m_mat.has_nan(); }

    /*
     * @brief Number of non-zeros
     * @details
     * @tparam Scalar_T
     * @return Number of non-zeros
     */
    template <typename Scalar_T>
    inline matrix_index_t arma_sparse_wrapper<Scalar_T>::nnz() const
    { return m_mat.n_nonzero; }

    /*
     * @brief Infinity norm
     * @details
     * @tparam Scalar_T
     * @return Inf norm
     */
    template <typename Scalar_T>
    inline Scalar_T arma_sparse_wrapper<Scalar_T>::norm_inf() const
    { return static_cast<Scalar_T>(arma::norm(m_mat, "inf")); }

    /*
     * @brief Squared Frobenius norm
     * @details
     * @tparam Scalar_T
     * @return Frob2 norm
     */
    template <typename Scalar_T>
    inline Scalar_T arma_sparse_wrapper<Scalar_T>::norm_frob2() const
    { return static_cast<Scalar_T>(arma::accu(arma::square(m_mat))); }

    /*
     * @brief Trace
     * @details
     * @tparam Scalar_T
     * @return Trace
     */
    template <typename Scalar_T>
    inline Scalar_T arma_sparse_wrapper<Scalar_T>::trace() const
    { return arma::trace(m_mat); }

    /*
     * @brief Eigenvalues
     * @details
     * @tparam Scalar_T
     * @return Eigenvalues
     */
    template <typename Scalar_T>
    inline std::vector<std::complex<double>> arma_sparse_wrapper<Scalar_T>::eigenvalues() const
    {
      throw std::runtime_error("Not implemented for sparse");
    }

    /*
     * @brief Trace of product
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Trace(A*B) / Dim
     */
    template <typename Scalar_T>
    template <typename Other>
    inline Scalar_T arma_sparse_wrapper<Scalar_T>::trace_product(const Other& other) const
    {
      using traits_t = numeric_traits<Scalar_T>;
      const auto dim = matrix_index_t(m_mat.n_rows);
      if (dim == 0)
        return Scalar_T(0);
      // Trace(AB) = sum_i sum_j A_ij * B_ji
      Scalar_T tr = Scalar_T(0);
      for (auto it = m_mat.begin(); it != m_mat.end(); ++it)
        tr += (*it) * other(it.col(), it.row());
      return tr / traits_t::to_scalar_t(dim);
    }

    /*
     * @brief Inner product
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Inner product
     */
    template <typename Scalar_T>
    template <typename Result_Scalar_T, typename Other>
    inline Result_Scalar_T arma_sparse_wrapper<Scalar_T>::inner(const Other& other) const
    {
      Result_Scalar_T sum = Result_Scalar_T(0);
      for (auto it = m_mat.begin(); it != m_mat.end(); ++it)
        sum += static_cast<Result_Scalar_T>(*it) * static_cast<Result_Scalar_T>(other(it.row(), it.col()));
      if (nbr_rows() == 0)
        return Result_Scalar_T(0);
      return sum / Result_Scalar_T(static_cast<double>(nbr_rows()));
    }

    template <typename Scalar_T>
    inline void arma_sparse_wrapper<Scalar_T>::similarity_transform(const std::vector<matrix_index_t>& perm,
                                                                    const std::vector<Scalar_T>& signs)
    {
      const matrix_index_t n = nbr_rows();
      std::vector<matrix_index_t> inv_perm(n);
      for (matrix_index_t i = 0; i < n; ++i)
        inv_perm[perm[i]] = i;

      arma::umat locations(2, m_mat.n_nonzero);
      arma::Col<Scalar_T> values(m_mat.n_nonzero);
      matrix_index_t count = 0;

      for (auto it = m_mat.begin(); it != m_mat.end(); ++it)
      {
        const matrix_index_t i = inv_perm[it.row()];
        const matrix_index_t j = inv_perm[it.col()];
        locations(0, count) = i;
        locations(1, count) = j;
        values[count] = signs[i] * (*it) * signs[j];
        count++;
      }
      m_mat = MatrixType(locations, values, n, n);
    }

    template <typename Scalar_T>
    inline void arma_sparse_wrapper<Scalar_T>::transpose_similarity_transform(const std::vector<matrix_index_t>& perm,
                                                                              const std::vector<Scalar_T>& signs)
    {
      const matrix_index_t n = nbr_rows();
      std::vector<matrix_index_t> inv_perm(n);
      for (matrix_index_t i = 0; i < n; ++i)
        inv_perm[perm[i]] = i;

      arma::umat locations(2, m_mat.n_nonzero);
      arma::Col<Scalar_T> values(m_mat.n_nonzero);
      matrix_index_t count = 0;

      for (auto it = m_mat.begin(); it != m_mat.end(); ++it)
      {
        const matrix_index_t i = inv_perm[it.col()];
        const matrix_index_t j = inv_perm[it.row()];
        locations(0, count) = i;
        locations(1, count) = j;
        values[count] = signs[i] * (*it) * signs[j];
        count++;
      }
      m_mat = MatrixType(locations, values, n, n);
    }

    template <typename Scalar_T>
    inline void arma_sparse_wrapper<Scalar_T>::involute()
    {
      for (auto it = m_mat.begin(); it != m_mat.end(); ++it)
        if (std::popcount((unsigned int)(it.row() ^ it.col())) % 2 != 0)
          (*it) *= -1.0;
    }

    template <typename Scalar_T>
    inline void arma_sparse_wrapper<Scalar_T>::reverse(index_t p, index_t q)
    {
      m_mat = m_mat.t();
    }
    /*
     * @brief Mixed Kronecker matrix product: Sparse x Dense -> Dense (wrapper)
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Kronecker product
     */
    template <typename Scalar_T>
    template <typename Other_Scalar_T>
    inline arma_matrix_wrapper<Other_Scalar_T> arma_sparse_wrapper<Scalar_T>::kron(
        const arma_matrix_wrapper<Other_Scalar_T>& other) const
    {
      arma_matrix_wrapper<Other_Scalar_T> result(nbr_rows() * other.nbr_rows(), nbr_cols() * other.nbr_cols());
      result.zeros();

      for (auto it = this->begin(); it != this->end(); ++it)
      {
        auto val = static_cast<Other_Scalar_T>(*it);
        result.m_mat.submat(it.row() * other.nbr_rows(), it.col() * other.nbr_cols(), (it.row() + 1) * other.nbr_rows() - 1,
                            (it.col() + 1) * other.nbr_cols() - 1) = val * other.m_mat;
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
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T> arma_sparse_wrapper<Scalar_T>::kron(const arma_sparse_wrapper<Scalar_T>& other) const
    {
      arma_sparse_wrapper<Scalar_T> result;
      result.m_mat = arma::kron(m_mat, other.m_mat);
      return result;
    }

    /*
     * @brief Monomial Kronecker matrix product of sparse wrappers
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Kronecker product
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T> arma_sparse_wrapper<Scalar_T>::mono_kron(const arma_sparse_wrapper<Scalar_T>& other) const
    {
      const matrix_index_t other_rows = other.nbr_rows();
      const matrix_index_t other_cols = other.nbr_cols();
      const matrix_index_t n_nz = m_mat.n_nonzero * other.m_mat.n_nonzero;

      arma::umat locations(2, n_nz);
      arma::Col<Scalar_T> values(n_nz);
      matrix_index_t count = 0;

      for (auto itA = m_mat.begin(); itA != m_mat.end(); ++itA)
      {
        auto rA = itA.row();
        auto cA = itA.col();
        auto vA = *itA;
        for (auto itB = other.m_mat.begin(); itB != other.m_mat.end(); ++itB)
        {
          locations(0, count) = rA * other_rows + itB.row();
          locations(1, count) = cA * other_cols + itB.col();
          values[count] = vA * (*itB);
          count++;
        }
      }
      arma_sparse_wrapper<Scalar_T> result;
      result.m_mat = MatrixType(locations, values, nbr_rows() * other_rows, nbr_cols() * other_cols);
      return result;
    }

    /*
     * @brief Monomial product of sparse wrappers
     * @details
     * @tparam Scalar_T
     * @param other Other matrix
     * @return Product
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T> arma_sparse_wrapper<Scalar_T>::mono_prod(const arma_sparse_wrapper<Scalar_T>& other) const
    {
      const matrix_index_t n_nz = other.m_mat.n_nonzero;

      arma::umat locations(2, n_nz);
      arma::Col<Scalar_T> values(n_nz);
      matrix_index_t count = 0;

      for (auto itB = other.m_mat.begin(); itB != other.m_mat.end(); ++itB)
      {
        auto rB = itB.row();
        auto cB = itB.col();
        auto vB = *itB;

        arma::uword col_start = m_mat.col_ptrs[rB];
        arma::uword col_end = m_mat.col_ptrs[rB + 1];
        if (col_start < col_end)
        {
          locations(0, count) = m_mat.row_indices[col_start];
          locations(1, count) = cB;
          values[count] = m_mat.values[col_start] * vB;
          count++;
        }
      }
      if (count < n_nz)
      {
        locations.resize(2, count);
        values.resize(count);
      }

      arma_sparse_wrapper<Scalar_T> result;
      result.m_mat = MatrixType(locations, values, nbr_rows(), nbr_cols());
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
    template <typename Scalar_T>
    template <typename RHS_T>
    inline RHS_T arma_sparse_wrapper<Scalar_T>::nork(const RHS_T& rhs, bool mono) const
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
      for (auto it = this->begin(); it != this->end(); ++it)
      {
        auto val = *it;
        result.m_mat += static_cast<typename RHS_T::value_type>(val)
                        * rhs.m_mat.submat(it.row() * blk_rows, it.col() * blk_cols, (it.row() + 1) * blk_rows - 1,
                                           (it.col() + 1) * blk_cols - 1);
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
    // unit_helper Specializations
    // =========================================================================

    // Specialization for Armadillo Wrapper
    template <typename Scalar_T>
    struct unit_helper<arma_matrix_wrapper<Scalar_T>>
    {
      static inline arma_matrix_wrapper<Scalar_T> apply(matrix_index_t dim)
      {
        arma_matrix_wrapper<Scalar_T> result(dim, dim);
        result.m_mat.eye();
        return result;
      }
    };

    // Specialization for Armadillo Sparse Wrapper
    template <typename Scalar_T>
    struct unit_helper<arma_sparse_wrapper<Scalar_T>>
    {
      static inline arma_sparse_wrapper<Scalar_T> apply(matrix_index_t dim)
      {
        arma_sparse_wrapper<Scalar_T> result(dim, dim);
        result.m_mat.eye(dim, dim);
        return result;
      }
    };

    /*
     * @brief Product of scalar and matrix wrapper
     * @details
     * @tparam Scalar_T
     * @param s Scalar
     * @param m Matrix
     * @return Product
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T> operator*(Scalar_T s, const arma_matrix_wrapper<Scalar_T>& m)
    {
      arma_matrix_wrapper<Scalar_T> result;
      result.m_mat = s * m.m_mat;
      return result;
    }

    /*
     * @brief Product of matrix wrapper and scalar
     * @details
     * @tparam Scalar_T
     * @param m Matrix
     * @param s Scalar
     * @return Product
     */
    template <typename Scalar_T>
    inline arma_matrix_wrapper<Scalar_T> operator*(const arma_matrix_wrapper<Scalar_T>& m, Scalar_T s)
    { return s * m; }

    /*
     * @brief Product of scalar and sparse wrapper
     * @details
     * @tparam Scalar_T
     * @param s Scalar
     * @param m Matrix
     * @return Product
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T> operator*(Scalar_T s, const arma_sparse_wrapper<Scalar_T>& m)
    {
      arma_sparse_wrapper<Scalar_T> result(m);
      result *= s;
      return result;
    }

    /*
     * @brief Product of sparse wrapper and scalar
     * @details
     * @tparam Scalar_T
     * @param m Matrix
     * @param s Scalar
     * @return Product
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T> operator*(const arma_sparse_wrapper<Scalar_T>& m, Scalar_T s)
    { return s * m; }

    /*
     * @brief Sum of sparse wrappers
     * @details
     * @tparam Scalar_T
     * @param lhs Left Hand Side
     * @param rhs Right Hand Side
     * @return Sum
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T> operator+(const arma_sparse_wrapper<Scalar_T>& lhs,
                                                   const arma_sparse_wrapper<Scalar_T>& rhs)
    {
      arma_sparse_wrapper<Scalar_T> result(lhs);
      result += rhs;
      return result;
    }

    /*
     * @brief Difference of sparse wrappers
     * @details
     * @tparam Scalar_T
     * @param lhs Left Hand Side
     * @param rhs Right Hand Side
     * @return Difference
     */
    template <typename Scalar_T>
    inline arma_sparse_wrapper<Scalar_T> operator-(const arma_sparse_wrapper<Scalar_T>& lhs,
                                                   const arma_sparse_wrapper<Scalar_T>& rhs)
    {
      arma_sparse_wrapper<Scalar_T> result(lhs);
      result.m_mat -= rhs.m_mat;  // Armadillo supports -=
      return result;
    }

  }  // namespace matrix
}  // namespace glucat
#ifdef GLUCAT_DOCTEST
#include <doctest/doctest.h>

#include <sstream>
#include <utility>

namespace glucat
{
  namespace matrix
  {

    template <typename Scalar_T>
    void test_arma_wrappers_templated()
    {
      using Matrix_T = arma_matrix_wrapper<Scalar_T>;
      using Sparse_T = arma_sparse_wrapper<Scalar_T>;

      SUBCASE("Dense Matrix: Constructors and Assignment")
      {
        Matrix_T m1(2, 2);
        m1.zeros();
        CHECK(m1.nbr_rows() == 2);
        CHECK(m1.nbr_cols() == 2);
        CHECK(m1.nnz() == 0);

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

      SUBCASE("Dense Matrix: Arithmetic Operators")
      {
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
        g(0, 0) = Scalar_T(1);
        g(0, 1) = Scalar_T(2);
        g(1, 0) = Scalar_T(3);
        g(1, 1) = Scalar_T(4);
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

      SUBCASE("Dense Matrix: Analysis and Solve")
      {
        Matrix_T m(2, 2);
        m(0, 0) = Scalar_T(2);
        m(0, 1) = Scalar_T(1);
        m(1, 0) = Scalar_T(1);
        m(1, 1) = Scalar_T(2);

        CHECK(m.is_finite());
        CHECK_FALSE(m.has_nan());

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

      SUBCASE("Sparse Matrix: Basic Operations")
      {
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
      }

      SUBCASE("Kronecker Product and Nork")
      {
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

        // Nork (Sparse)
        auto ss = s.kron(s);
        CHECK(ss.nbr_rows() == 4);
        CHECK(ss.nnz() == 4);

        auto qs = s.nork(ss, true);
        CHECK(qs.nbr_rows() == 2);
      }

      SUBCASE("Norms and Inner Product")
      {
        Matrix_T m(2, 2);
        m(0, 0) = Scalar_T(1);
        m(0, 1) = Scalar_T(-2);
        m(1, 0) = Scalar_T(3);
        m(1, 1) = Scalar_T(4);

        CHECK(m.norm_inf() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(7))));
        CHECK(m.norm_frob2() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1 + 4 + 9 + 16))));
        CHECK(m.isnan() == false);
        CHECK(m.isinf() == false);

        Sparse_T s(2, 2);
        s(0, 0) = Scalar_T(1);
        s(1, 1) = Scalar_T(4);
        CHECK(s.norm_inf() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(4))));
        CHECK(s.norm_frob2() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1 + 16))));
        CHECK(s.isnan() == false);
        CHECK(s.isinf() == false);

        // Sparse iterator operator!=
        Sparse_T s3(2, 2);
        s3(0, 0) = Scalar_T(1);
        s3(1, 0) = Scalar_T(2);
        auto it1 = s3.begin();
        auto it2 = s3.begin();
        ++it2;
        CHECK(it1 != it2);
        // Inner product
        auto in = m.template inner<Scalar_T>(m);
        CHECK(in == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1 + 4 + 9 + 16) / Scalar_T(2))));

        auto ins = s.template inner<Scalar_T>(s);
        CHECK(ins == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1 + 16) / Scalar_T(2))));
      }

      SUBCASE("Interop and Utility")
      {
        Sparse_T s(2, 2);
        s.unit(2, 2);

        // Sparse to Dense
        Matrix_T d(s);
        CHECK(d.trace() == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(2))));

        // operator<<
        std::ostringstream oss;
        oss << s;
        CHECK(!oss.str().empty());

        // Empty iterator
        Sparse_T empty(0, 0);
        CHECK(empty.begin() == empty.end());
      }

      SUBCASE("Cross-backend construction")
      {
        using Eigen_T = eigen_matrix_wrapper<Scalar_T>;
        Eigen_T e_mat(2, 2);
        e_mat.zeros();
        e_mat(0, 0) = Scalar_T(1);
        e_mat(1, 1) = Scalar_T(2);

        Matrix_T a_mat(e_mat);
        CHECK(a_mat.nbr_rows() == 2);
        CHECK(a_mat.nbr_cols() == 2);
        CHECK(a_mat(0, 0) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(1))));
        CHECK(a_mat(1, 1) == doctest::Approx(numeric_traits<Scalar_T>::to_double(Scalar_T(2))));
        CHECK(a_mat.nnz() == 2);
      }
    }

  }  // namespace matrix
}  // namespace glucat

TEST_CASE("matrix::arma_wrappers")
{
  SUBCASE("float")
  {
    glucat::matrix::test_arma_wrappers_templated<float>();
  }
  SUBCASE("double")
  {
    glucat::matrix::test_arma_wrappers_templated<double>();
  }
}
#endif

#endif  // _GLUCAT_USE_ARMADILLO
#endif  // _GLUCAT_MATRIX_ARMA_IMP_H
