#ifndef GLUCAT_TEST_PEG02_H
#define GLUCAT_TEST_PEG02_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg02.cpp : programming example 2 : inner, outer and geometric products, contraction
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

namespace peg02
{
  using namespace glucat;
  using namespace std;
  template< class Multivector_T >
  void do_test2()
  {
    typedef Multivector_T number;
    const number  a("2{1}+3{2}+3{3}"),
                  b("3{1}+2{2}+4{3}"),
                  c("5+{1}+3{2,3}");
    a.write("a =");
    b.write("b =");
    c.write("c =");

    (a & b).write("a&b =");
    (a % b).write("a%b =");
    (a ^ b).write("a^b =");
    (a * b).write("a*b =");

    ((a & b) + (a ^ b)).write("a&b + a^b =");
    ((a % b) + (a ^ b)).write("a%b + a^b =");
    cout << endl;
  }
}

int test02()
{
  using namespace peg02;
  cout << "Programming example 2 : inner, outer and geometric products, contraction" << endl;
  cout << endl;
  cout << "Let a, b be elements in Cl(3,0) where" << endl;
  cout << "a = 2{1} + 3{2} + 3{3} and b = 3{1} + {2} + 4{3}." << endl;
  cout << "Determine ab, a.b, a^b; Verify that ab == a.b + a^b." << endl;
  cout << endl;
  cout << "framed_multi<double>" << endl;
  do_test2< framed_multi<double> >();
  cout << "matrix_multi<double>" << endl;
  do_test2< matrix_multi<double> >();
  return 0;
}
#endif
