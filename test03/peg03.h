#ifndef GLUCAT_TEST_PEG03_H
#define GLUCAT_TEST_PEG03_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg03.cpp : programming example 3 : Powers, output to file
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

namespace peg03
{
  using namespace glucat;
  using namespace std;

  template< class Multivector_T >
  static
  void 
  do_test3()
  {
    typedef Multivector_T number;
    const number e1("{1}"), e2("{2}"), e3("{-1}");

    const number x = 2.0*e1 + 6.2*e2 - 4.1*e2*e3;
    const number y = 7.0*e1 - 2.1*e2 + 9.6*e1*e2 + 6.0*e2*e3;

    x.write("x = ");  y.write("y =");
    cout << endl;
    (2.0*x + 3.0*y).write("2x + 3y =");
    (x.pow(5)).write("x to the power 5 =");
    (y.outer_pow(3)).write("y to the outer power 3 =");

    const number ans( (x*y).even() );
    ans.write("even part of xy = ");
    cout << endl;

    ofstream out_file;
    out_file.open("eg3.res");
    ans.write(out_file, "even part of xy = ");
    out_file.close();
  }
}

int test03();

#endif
