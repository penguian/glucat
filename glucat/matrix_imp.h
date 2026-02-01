#ifndef _GLUCAT_MATRIX_IMP_H
#define _GLUCAT_MATRIX_IMP_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix.h : Declare common matrix classes and functions
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
#include "glucat/matrix_base_imp.h"

#include "glucat/matrix_eigen_imp.h"
#if defined(_GLUCAT_USE_ARMADILLO)
#include "glucat/matrix_arma_imp.h"
#endif

namespace glucat { namespace matrix
{

  /**
   * @brief Constructor from other matrix type
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   */
  template< typename Scalar_T >
  template< typename Other_Matrix_T >
  inline
  dense_matrix<Scalar_T>::
  dense_matrix(const Other_Matrix_T& other)
  : Base(other)
  { }

  /**
   * @brief Constructor from other matrix type
   * @details
   * @tparam Scalar_T
   * @param other Other matrix
   */
  template< typename Scalar_T >
  template< typename Other_Matrix_T >
  inline
  sparse_matrix<Scalar_T>::
  sparse_matrix(const Other_Matrix_T& other)
  : Base(other)
  { }

} }

#endif // _GLUCAT_MATRIX_IMP_H
