#ifndef _GLUCAT_GENERATION_H
#define _GLUCAT_GENERATION_H
/***************************************************************************
	  GluCat : Generic library of universal Clifford algebra templates
    generation.h : Declare functions for generation of the matrix representation
                             -------------------
    begin                : Wed Jan 23 2002
    copyright            : (C) 2002 by Paul C. Leopardi
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

namespace glucat
{
  /// Modulo function which works reliably for lhs < 0
  const int
  pos_mod(int lhs, int rhs);

  /// A signature is a pair of indices, p, q, with p == frame.max(), q == -frame.min()
  typedef std::pair<index_t, index_t> signature_t;

  /// Table of generators for specific signatures
  template< class Matrix_T >
  class generator_table :
  private std::map< signature_t, std::vector<Matrix_T> >
  {
  public:
    /// Pointer to generators for a specific signature
    const Matrix_T* operator() (const index_t p, const index_t q);
    /// Single instance of generator table
    friend generator_table& generator<Matrix_T>();
  private:
    /// Construct a vector of generators for a specific signature
    const std::vector<Matrix_T>& gen_vector(const index_t p, const index_t q);
    /// Construct generators for p,q given generators for p-1,q-1
    void gen_from_pm1_qm1(const std::vector<Matrix_T>& old, const signature_t sig);
    /// Construct generators for p,q given generators for p-4,q+4
    void gen_from_pm4_qp4(const std::vector<Matrix_T>& old, const signature_t sig);
    /// Construct generators for p,q given generators for p+4,q-4
    void gen_from_pp4_qm4(const std::vector<Matrix_T>& old, const signature_t sig);
    /// Construct generators for p,q given generators for q+1,p-1
    void gen_from_qp1_pm1(const std::vector<Matrix_T>& old, const signature_t sig);

    // Enforce singleton
    // Reference: A. Alexandrescu, "Modern C++ Design", Chapter 6
    generator_table() {}
    ~generator_table() {}
    generator_table(const generator_table&);
    generator_table& operator= (const generator_table&);
  };

  /// Single instance of generator table
  template< class Matrix_T >
  generator_table<Matrix_T>&
  generator();

  /// Determine the matrix dimension of the fold of a frame
  template< typename Matrix_Index_T, const index_t LO, const index_t HI >
  const Matrix_Index_T
  folded_dim( const index_set<LO,HI>& frm );
}
#endif  // _GLUCAT_GENERATION_H
