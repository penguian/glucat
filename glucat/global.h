#ifndef _GLUCAT_GLOBAL_H
#define _GLUCAT_GLOBAL_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    global.h : Global declarations
                             -------------------
    begin                : Sun Dec 9 2001
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
 "Clifford algebras with numeric and symbolic computations, Birkhauser, 1996."
 ***************************************************************************
 See also Arvind Raja's original header comments and references in glucat.h
 ***************************************************************************/

namespace glucat
{
  /// Compile time assertion
  // Reference: A. Alexandrescu, "Modern C++ Design", Addison-Wesley, 2001, p25
  template<bool> struct CTAssertion;
  template<> struct CTAssertion<true> { };
  #define _GLUCAT_CTAssert(expr, msg) namespace { CTAssertion<(expr)> ERROR_##msg; }

  // Global types which determine sizes
  /// Size of index_t should be enough to represent LO, HI
  typedef short         index_t;
  /// Size of set_value_t should be enough to contain index_set<LO,HI>
  typedef unsigned long set_value_t;

  // Global constants
  const double MS_PER_S = 1000.0;

  // Constants which determine sizes
  /// If radix of unsigned char is not 2, we can't easily determine number of bits from sizeof
  _GLUCAT_CTAssert(std::numeric_limits<unsigned char>::radix == 2, CannotDetermineBitsPerChar);

  /// Number of bits per char is used to determine number of bits in set_value_t
  const index_t BITS_PER_CHAR = std::numeric_limits<unsigned char>::digits;

  /// Number of bits in set_value_t
  const index_t BITS_PER_SET_VALUE = BITS_PER_CHAR * index_t(sizeof(set_value_t));

  // Constants which are determined by size
  /// Default lowest index in an index set
  const index_t DEFAULT_LO = -index_t(BITS_PER_SET_VALUE)/2 + 1;
  /// Default highest index in an index set
  const index_t DEFAULT_HI =  index_t(BITS_PER_SET_VALUE)/2 - 1;

  /// Default for truncation
  const double DEFAULT_TRUNCATION = std::numeric_limits<float>::epsilon();

  /// Tuning policy
  template
  <
    int Mult_Matrix_Threshold = 9,
    int Div_Max_Steps = 3,
    int Sqrt_Max_Steps = 7,
    int Log_Max_Outer_Steps = 256,
    int Log_Max_Inner_Steps = 8
  >
  struct tuning
  {
  // Tuning for multiplication
    /// Minimum index count needed to invoke matrix multiplication algorithm
    enum { mult_matrix_threshold = Mult_Matrix_Threshold };
  // Tuning for division
    /// Maximum steps of iterative refinement in division algorithm
    enum { div_max_steps = Div_Max_Steps };
  // Tuning for sqrt
    /// Maximum number of steps in square root iteration
    enum { sqrt_max_steps = Sqrt_Max_Steps };
  // Tuning for log
    /// Maximum number of incomplete square roots in cascade log algorithm
    enum { log_max_outer_steps = Log_Max_Outer_Steps };
    /// Maximum number of steps in incomplete square root within cascade log algorithm
    enum { log_max_inner_steps = Log_Max_Inner_Steps };
  };
}
#endif // _GLUCAT_GLOBAL_H
