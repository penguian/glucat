#ifndef GLUCAT_TEST_PEG14_H
#define GLUCAT_TEST_PEG14_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg14.cpp : programming example 14 : framed_multi <-> matrix_multi
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

int test14()
{
  cout << "Programming example 14 : framed_multi <-> matrix_multi" << endl;
  using namespace glucat;
  typedef matrix_multi<long double> cm;
  typedef framed_multi<long double> cf;
  cf a("{-3}+{-2}+{-1}");
  cf b("1.e4{-1}+1.e-4{-2}");
  cf c;
  cm::index_set_t sub = a.frame() | b.frame();
  cm A( a, sub );
  cm B( b, sub );
  cm C;
  a.write("a =");
  A.write("A =");
  b.write("b =");
  B.write("B =");
  B.truncated().write("B.truncated()");
  (B - b).write("B - b =");
  (b - B).write("b - B =");
  c = (a * b) / b;
  c.write("c = (a * b) / b =");
  C = c;
  C.write("C = c =");
  (C - c).write("C - c =");
  (c - C).write("c - C =");

  cout << "star(a, b)    = " << star(a, b)    << endl;
  cout << "scalar(a * b) = " << scalar(a * b) << endl;
  cout << "star(b, a)    = " << star(b, a)    << endl;
  cout << "scalar(b * a) = " << scalar(b * a) << endl;
  cout << "star(A, B)    = " << star(A, B)    << endl;;
  cout << "scalar(A * B) = " << scalar(A * B) << endl;
  cout << "star(B, A)    = " << star(B, A)    << endl;;
  cout << "scalar(B * A) = " << scalar(B * A) << endl;
  return 0;
}
#endif
