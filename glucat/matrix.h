#ifndef _GLUCAT_MATRIX_H
#define _GLUCAT_MATRIX_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix.h : Declare common matrix functions
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001 by Paul C. Leopardi
    email                : leopardi@bigpond.net.au
 ***************************************************************************
 *   This library is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2.1 of the  *
 *   License, or (at your option) any later version.                       *
 *   See http://www.fsf.org/copyleft/lesser.html for details               *
 ***************************************************************************
 This library is based on a prototype written by Arvind Raja and was
 licensed under the LGPL with permission of the author. See Arvind Raja,
 "Object-oriented implementations of Clifford algebras in C++: a prototype",
 in Ablamowicz, Lounesto and Parra (eds.)
 "Clifford algebras with numeric and symbolic computations", Birkhauser, 1996.
 ***************************************************************************
     See also Arvind Raja's original header comments in glucat.h
 ***************************************************************************/

using std::vector;

namespace glucat
{
  /// Kronecker tensor product of matrices - as per Matlab kron
  template< class Matrix_T >
  void
  kron(const Matrix_T& x, const Matrix_T& y, Matrix_T& z);

  /// Unit matrix - as per Matlab eye
  template< class Matrix_T >
  void
  unit(int n, Matrix_T& result);

  /// Does a matrix have only one non-zero per row (or column) ?
  // Note: This should return false for dense matrices,
  // but this is OK for the current use of this function,
  // which is to help determine when a compressed multiply should have a dense result
  template< class Matrix_T >
  bool
  is_singlet(const Matrix_T& x);

  /// Does a matrix have only one non-zero per row and column ?
  // Note: This should return false for dense matrices,
  // but this is OK for the current use of this function,
  // which is to help determine when a compressed multiply should have a dense result
  template< class Matrix_T >
  bool
  is_perm_shaped(const Matrix_T& x);

  /// Generate the next generation of generators, given the previous generation
  template< class Matrix_T >
  void
  gen(const vector<Matrix_T>& old, vector< vector< Matrix_T > >& generators);

  /// Generate a specific generation of generators, given the whole family tree
  template< class Matrix_T >
  const vector<Matrix_T>&
  gengen(const int n, vector< vector<Matrix_T> >& generators);

  /// Inner product: sum(x(i,j)*y(i,j))/x.nrows()
  template< class Matrix_T, class Scalar_T >
  Scalar_T
  inner(const Matrix_T& x, const Matrix_T& y);
}
#endif  // _GLUCAT_MATRIX_H
