#ifndef GLUCAT_TEST_PEG04_H
#define GLUCAT_TEST_PEG04_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg04.cpp : programming example 4 : Rotations and electromagnetic fields
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

namespace peg04
{
  using namespace glucat;
  using namespace std;

  template< class Multivector_T >
  static
  void 
  do_test4()
  {
    typedef Multivector_T number;
    const number p1("{1}"), p2("{2}"), p3("{3}");
    const number q1("{-1}");
    const number j("{1,2,3}");

    cout << "Rotations in three dimensions" << endl;
    {
      const number x = 2.0*p1 + 3.0*p2 + 5.0*p3;
      x.write("x =");
      const number a =     p1 - 3.0*p2 + 4.0*p3;
      a.write("a =");
      const number A = j*a;
      A.write("A =");
      const number s = exp(A / 2.0);
      s.write("s = exp(A/2) =");
      cout << "Rotation about the axis a" << endl;
      const number y = s*x/s;
      y.write("y = s*x/s =");
      (y*y).write("y*y =");
      (x*x).write("x*x =");
      cout << endl;
    }
    cout << "Lorentz transformations of space-time" << endl;
    {
      const number x = 2.0*p1 + 3.0*p2 + 5.0*p3 + 2.0*q1;
      x.write("x =");
      const number A = 4.0*p1*p2 + 3.0*p1*p3 + 3.0*p1*q1 + p2*p3 + p2*q1 + 2.0*p3*q1;
      A.write("A =");
      const number s = exp(A/2.0);
      s.write("s = exp(A/2) =");
      const number y = s*x/involute(s);
      y.write("y = s*x/involute(s) =");
      (y*y).write("y*y =");
      (x*x).write("x*x =");

      const number F = 5.0*p1*p2 + p1*p3 + 2.0*p1*q1 + 7.0*p2*p3 + p2*q1 + 2.0*p3*q1;
      F.write("F =");
      const number G = s*F/s;
      G.write("G = s*F/s =");
      (G*G).write("G*G =");
      (F*F).write("F*F =");
      (-F*q1*F / 2.0).write("-F*q1*F / 2 =");
      cout << endl;
    }
    cout << "Electromagnetism in Cl_3" << endl;
    cout << "Let F = E - jB be the electromagnetic field, where the electric component" << endl;
    cout << "E = {1}+2{2}+4{3}, the magnetic component B = 3{1}+5{2}+7{3}." << endl;
    cout << "j ={1,2,3} is the unit director of R_3." << endl;
    cout << "Compute the Lorentz invariants, energy density and Poynting vector of F." << endl;
    number E =     p1 + 2.0*p2 + 4.0*p3;
    E.write("E =");
    number B = 3.0*p1 + 5.0*p2 + 7.0*p3;
    B.write("B =");
    {
      number F = E - j*B;
      F.write("Electromagnetic field: F = E - j*B =");
      (F*F / 2.0).write
        ("Lorentz invariants: F*F / 2.0 =");
      // Correction below comes from Clical demo http://www.teli.stadia.fi/~lounesto/CLICAL.zip
      (-involute(F)*F / 2.0).write
        ("Energy density and Poynting vector: -involute(F)*F / 2.0 =");
      cout << "Boost at half the velocity of light" << endl;
      const number v = 0.5*p1;
      v.write("v =");
      const number A = atanh(v);
      A.write("A =");

      const number s = exp(A / 2.0);
      s.write("s = exp(A / 2.0) =");

      const number G = s * F / s;
      G.write("Electric field: G = s * F / s =");

      (j * G(2)).write("Magnetic induction:");

      const number x = 10.0 + p1 + p2;
      x.write("x =");
      const number y = s * x / involute(s);
      y.write("y = s * x / involute(s) =");
      (y * conj(y)).write("y * conj(y) =");
      (x * conj(x)).write("x * conj(x) = ");
      cout << endl;
    }
    cout << "Electromagnetism in Cl_{3,1}" << endl;
    {
      const number i = j * q1;
      i.write("i =");
      E *= q1;
      E.write("E =");
      B *= q1;
      B.write("B =");
      const number F = E - i*B;
      F.write("F = E - i*B =");
      (F*F / 2.0).write
        ("Lorentz invariants: F*F / 2.0 =");
      (-F * q1 * F / 2.0).write
        ("Energy density and Poynting vector: -F * q1 * F / 2.0 =");
      cout << endl;
    }
  }
}

int test04();

#endif
