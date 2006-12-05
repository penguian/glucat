/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg06.cpp : Driver for test 06
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2003 by Paul C. Leopardi
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

#include "test/driver.h"
#include "test06/peg06.h"

int test06()
{
  using namespace peg06;
  cout << "Programming example 6 : outer exponential" << endl;
  cout << endl;
  cout << "Compute the outer exponential of the bivector B = q1*q2 + q3*q4 + q5*q6" << endl;
  cout << endl;
  cout << "framed_multi<double>" << endl;
  do_test6< framed_multi<double> >();
  cout << "matrix_multi<double>" << endl;
  do_test6< matrix_multi<double> >();
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return try_catch(test06); }
