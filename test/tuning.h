#ifndef GLUCAT_TEST_TUNING_H
#define GLUCAT_TEST_TUNING_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    tuning.h : Use preprocessor macros to control test tuning
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2016 by Paul C. Leopardi
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

// If radix of int is not 2, we can't easily set thresholds
_GLUCAT_CTAssert(std::numeric_limits<unsigned int>::radix == 2, CannotSetThresholds)
const unsigned int Test_Tuning_Int_Digits = std::numeric_limits<int>::digits;
const unsigned int Test_Tuning_Max_Threshold = 1 << Test_Tuning_Int_Digits;

typedef glucat::precision_t precision_t;

#define __TEST_TUNING_DEFAULT_CONSTANT(SUFFIX) \
const unsigned int Test_Tuning_##SUFFIX                = glucat::DEFAULT_##SUFFIX

// Tuning policy constants
#if defined ( _GLUCAT_TEST_TUNING_SLOW )
const unsigned int Test_Tuning_Mult_Matrix_Threshold   = Test_Tuning_Max_Threshold;
__TEST_TUNING_DEFAULT_CONSTANT(Div_Max_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Sqrt_Max_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Log_Max_Outer_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Log_Max_Inner_Steps);
const unsigned int Test_Tuning_Basis_Max_Count         =       0;
const unsigned int Test_Tuning_Fast_Size_Threshold     = Test_Tuning_Max_Threshold;
const unsigned int Test_Tuning_Inv_Fast_Dim_Threshold  = Test_Tuning_Max_Threshold;
const unsigned int Test_Tuning_Products_Size_Threshold = Test_Tuning_Max_Threshold;
const precision_t  Test_Tuning_Function_Precision      = glucat::DEFAULT_Function_Precision;
#elif defined ( _GLUCAT_TEST_TUNING_NAIVE )
const unsigned int Test_Tuning_Mult_Matrix_Threshold   =       0;
__TEST_TUNING_DEFAULT_CONSTANT(Div_Max_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Sqrt_Max_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Log_Max_Outer_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Log_Max_Inner_Steps);
const unsigned int Test_Tuning_Basis_Max_Count         = Test_Tuning_Max_Threshold;
const unsigned int Test_Tuning_Fast_Size_Threshold     = Test_Tuning_Max_Threshold;
const unsigned int Test_Tuning_Inv_Fast_Dim_Threshold  = Test_Tuning_Max_Threshold;
__TEST_TUNING_DEFAULT_CONSTANT(Products_Size_Threshold);
const precision_t  Test_Tuning_Function_Precision      = glucat::DEFAULT_Function_Precision;
#elif defined ( _GLUCAT_TEST_TUNING_FAST )
const unsigned int Test_Tuning_Mult_Matrix_Threshold   =       0;
const unsigned int Test_Tuning_Div_Max_Steps           =       0;
const unsigned int Test_Tuning_Sqrt_Max_Steps          =      16;
const unsigned int Test_Tuning_Log_Max_Outer_Steps     =      16;
const unsigned int Test_Tuning_Log_Max_Inner_Steps     =       8;
const unsigned int Test_Tuning_Basis_Max_Count         =       1;
const unsigned int Test_Tuning_Fast_Size_Threshold     =       0;
const unsigned int Test_Tuning_Inv_Fast_Dim_Threshold  =       0;
const unsigned int Test_Tuning_Products_Size_Threshold =       0;
const precision_t  Test_Tuning_Function_Precision      = glucat::DEFAULT_Function_Precision;
#elif defined ( _GLUCAT_TEST_TUNING_PROMOTED )
__TEST_TUNING_DEFAULT_CONSTANT(Mult_Matrix_Threshold);
__TEST_TUNING_DEFAULT_CONSTANT(Div_Max_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Sqrt_Max_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Log_Max_Outer_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Log_Max_Inner_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Basis_Max_Count);
__TEST_TUNING_DEFAULT_CONSTANT(Fast_Size_Threshold);
__TEST_TUNING_DEFAULT_CONSTANT(Inv_Fast_Dim_Threshold);
__TEST_TUNING_DEFAULT_CONSTANT(Products_Size_Threshold);
const precision_t  Test_Tuning_Function_Precision      = glucat::precision_promoted;
#elif defined ( _GLUCAT_TEST_TUNING_DEMOTED )
__TEST_TUNING_DEFAULT_CONSTANT(Mult_Matrix_Threshold);
__TEST_TUNING_DEFAULT_CONSTANT(Div_Max_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Sqrt_Max_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Log_Max_Outer_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Log_Max_Inner_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Basis_Max_Count);
__TEST_TUNING_DEFAULT_CONSTANT(Fast_Size_Threshold);
__TEST_TUNING_DEFAULT_CONSTANT(Inv_Fast_Dim_Threshold);
__TEST_TUNING_DEFAULT_CONSTANT(Products_Size_Threshold);
const precision_t  Test_Tuning_Function_Precision      = glucat::precision_demoted;
#else
__TEST_TUNING_DEFAULT_CONSTANT(Mult_Matrix_Threshold);
__TEST_TUNING_DEFAULT_CONSTANT(Div_Max_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Sqrt_Max_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Log_Max_Outer_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Log_Max_Inner_Steps);
__TEST_TUNING_DEFAULT_CONSTANT(Basis_Max_Count);
__TEST_TUNING_DEFAULT_CONSTANT(Fast_Size_Threshold);
__TEST_TUNING_DEFAULT_CONSTANT(Inv_Fast_Dim_Threshold);
__TEST_TUNING_DEFAULT_CONSTANT(Products_Size_Threshold);
const precision_t  Test_Tuning_Function_Precision      = glucat::DEFAULT_Function_Precision;
#endif

/// Tuning policy
typedef glucat::tuning
  <
    Test_Tuning_Mult_Matrix_Threshold,
    Test_Tuning_Div_Max_Steps,
    Test_Tuning_Sqrt_Max_Steps,
    Test_Tuning_Log_Max_Outer_Steps,
    Test_Tuning_Log_Max_Inner_Steps,
    Test_Tuning_Basis_Max_Count,
    Test_Tuning_Fast_Size_Threshold,
    Test_Tuning_Inv_Fast_Dim_Threshold,
    Test_Tuning_Products_Size_Threshold,
    Test_Tuning_Function_Precision
  > Tune_P;

#undef __TEST_TUNING_DEFAULT_CONSTANT

#endif // GLUCAT_TEST_TUNING_H
