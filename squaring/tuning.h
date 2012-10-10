#ifndef GLUCAT_TEST_SQUARING_TUNING_H
#define GLUCAT_TEST_SQUARING_TUNING_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    tuning.h : Tuning for squaring test
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2012 by Paul C. Leopardi
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

// Tuning policy constants
#if defined ( _GLUCAT_SQUARING_SLOW )
const unsigned int Driver_Mult_Matrix_Threshold  =      64;
const unsigned int Driver_Basis_Max_Count        =       1;
const unsigned int Driver_Fast_Size_Threshold    = 1 << 30;
const unsigned int Driver_Inv_Fast_Dim_Threshold = 1 << 30;
#elif defined ( _GLUCAT_SQUARING_MEDIUM )
const unsigned int Driver_Mult_Matrix_Threshold  =       0;
const unsigned int Driver_Basis_Max_Count        =      64;
const unsigned int Driver_Fast_Size_Threshold    = 1 << 30;
const unsigned int Driver_Inv_Fast_Dim_Threshold = 1 << 30;
#elif defined ( _GLUCAT_SQUARING_FAST )
const unsigned int Driver_Mult_Matrix_Threshold  =       0;
const unsigned int Driver_Basis_Max_Count        =       1;
const unsigned int Driver_Fast_Size_Threshold    = 1 <<  0;
const unsigned int Driver_Inv_Fast_Dim_Threshold = 1 <<  0;
#elif defined ( _GLUCAT_SQUARING_OPT )
const unsigned int Driver_Mult_Matrix_Threshold  =       8;
const unsigned int Driver_Basis_Max_Count        =      10;
const unsigned int Driver_Fast_Size_Threshold    = 1 << 10;
const unsigned int Driver_Inv_Fast_Dim_Threshold = 1 <<  7;
#else
const unsigned int Driver_Mult_Matrix_Threshold  = glucat::DEFAULT_Mult_Matrix_Threshold;
const unsigned int Driver_Basis_Max_Count        = glucat::DEFAULT_Basis_Max_Count;
const unsigned int Driver_Fast_Size_Threshold    = glucat::DEFAULT_Fast_Size_Threshold;
const unsigned int Driver_Inv_Fast_Dim_Threshold = glucat::DEFAULT_Inv_Fast_Dim_Threshold;
#endif

/// Tuning policy
typedef glucat::tuning
  <
    Driver_Mult_Matrix_Threshold,
    glucat::DEFAULT_Div_Max_Steps,
    glucat::DEFAULT_Sqrt_Max_Steps,
    glucat::DEFAULT_Log_Max_Outer_Steps,
    glucat::DEFAULT_Log_Max_Inner_Steps,
    Driver_Basis_Max_Count,
    Driver_Fast_Size_Threshold,
    Driver_Inv_Fast_Dim_Threshold
  > Tune_P;

#endif // GLUCAT_TEST_SQUARING_TUNING_H
