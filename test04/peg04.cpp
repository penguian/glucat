/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg04.cpp : Driver for test 04
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
#include "test04/peg04.h"

int test04()
{
  using namespace peg04;
  cout << "Programming example 4 : Rotations and electromagnetic fields" << endl;
  cout << endl;
  cout << "Demonstration calculations from P. Lounesto," << endl;
  cout << "''Clifford algebra calculations with a microcomputer'' in" << endl;
  cout << "A. Micali, R. Boudet, J. Helmstetter, (eds)," << endl;
  cout << "Clifford algebras and their applications in mathematical physics" << endl;
  cout << "Kluwer, 1992, pp44-46." << endl;
  cout << "See also: P. Lounesto, et. al., Clical demo" << endl;
  cout << "http://www.teli.stadia.fi/~lounesto/CLICAL.htm" << endl;
  cout << endl;
  cout << "framed_multi<double>" << endl;
  do_test4< framed_multi<double> >();
  cout << "matrix_multi<double>" << endl;
  do_test4< matrix_multi<double> >();
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return try_catch(test04); }
