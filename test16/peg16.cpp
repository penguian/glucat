/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg16.cpp : Driver for test 16
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
#include <boost/numeric/ublas/io.hpp>
#include "test16/peg16.h"

int test16()
{
  using namespace peg16;
  cout << "Programming example 16 : Matrices of multivectors" << endl;
  cout << endl;
  cout << "framed_multi<double,-3,3>" << endl;
  do_test16< framed_multi<double,-3,3> >();
  cout << "matrix_multi<double,-3,3>" << endl;
  do_test16< matrix_multi<double,-3,3> >();
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return try_catch(test16); }
