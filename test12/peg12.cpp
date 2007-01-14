/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg12.cpp : Driver for test 12
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2007 by Paul C. Leopardi
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
#include "test12/peg12.h"

int test12()
{
  using namespace peg12;
  cout << "Programming example 12 : Frame group multiplication tables" << endl;
  cout << endl;
  cout << "framed_multi<float,-2,2>" << endl;
  do_test12< framed_multi<float,-2,2> >();
  cout << "matrix_multi<float,-2,2>" << endl;
  do_test12< matrix_multi<float,-2,2> >();
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return try_catch(test12); }
