#ifndef GLUCAT_TEST_DRIVER_H
#define GLUCAT_TEST_DRIVER_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    driver.h : Header for example and timing test driver
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

#include "glucat/glucat.h"
const int DRIVER_BASIS_MAX_COUNT = 9;
typedef glucat::tuning
  <
    glucat::DEFAULT_Mult_Matrix_Threshold,
    glucat::DEFAULT_Div_Max_Steps,
    glucat::DEFAULT_Sqrt_Max_Steps,
    glucat::DEFAULT_Log_Max_Outer_Steps,
    glucat::DEFAULT_Log_Max_Inner_Steps,
    DRIVER_BASIS_MAX_COUNT,
    glucat::DEFAULT_Fast_Size_Threshold,
    glucat::DEFAULT_Inv_Fast_Dim_Threshold
  >
  Tune_P;
#include "glucat/glucat_imp.h"
#include <stdio.h>
#include "test/peg01.h"
#include "test/peg02.h"
#include "test/peg03.h"
#include "test/peg04.h"
#include "test/peg05.h"
#include "test/peg06.h"
#include "test/peg07.h"
#include "test/peg08.h"
#include "test/peg09.h"
#include "test/peg10.h"
#include "test/peg11.h"
#include "test/peg12.h"
#include "test/peg13.h"
#include "test/peg14.h"
#include "test/peg15.h"
#include "test/peg16.h"
#include "test/squaring.h"
#endif // GLUCAT_TEST_DRIVER_H
