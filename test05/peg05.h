#ifndef GLUCAT_TEST_PEG05_H
#define GLUCAT_TEST_PEG05_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg05.cpp : programming example 5 : Octonions
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

namespace peg05
{
  using namespace glucat;
  using namespace std;
  template< class Multivector_T >
  Multivector_T operator| (const Multivector_T& lhs, const Multivector_T& rhs)
  {
    typedef Multivector_T number;
    const number o = number("1+{-4,-2,-1}") * number("1+{-5,-3,-2}") *
                     number("1+{-6,-4,-3}") * number("1+{-7,-5,-4}");
    const number result = lhs * rhs * o;
    return result(0) + result(1);
  }

  template< class Multivector_T >
  void do_test5()
  {
    typedef Multivector_T number;
    const number a("3+{-1}+4{-2}"),
                 b("2+3{-2}+5{-3}");
    a.write("a = ");
    b.write("b = ");
    (a | b).write("octonion product of a and b = ");
    cout << endl;
  }
}

int test05()
{
  using namespace peg05;
  cout << "Programming example 5 : Octonions" << endl;
  cout << endl;
  cout << "In the Cayley algebra (=Cl(0,7)) compute the octonion product ab, where" << endl;
  cout << "a = 3 + e1 + 4e2 and b = 2 + 3e2 + 5e3." << endl;
  cout << endl;
  cout << "framed_multi<double>" << endl;
  do_test5< framed_multi<double> >();
  cout << "matrix_multi<double>" << endl;
  do_test5< matrix_multi<double> >();
  return 0;
}
#endif
