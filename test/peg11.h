#ifndef GLUCAT_TEST_PEG11_H
#define GLUCAT_TEST_PEG11_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg11.cpp : programming example 11 : Square root and transcendental functions
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

namespace peg11
{
  using namespace glucat;
  template< class Multivector_T >
  void transcendtest(const Multivector_T& A)
  {
    typedef Multivector_T m_;
    typedef typename m_::scalar_t s_;
    const streamsize prec = cout.precision(numeric_limits<s_>::digits10);
    A.write("A");
    sqrt(A).write("sqrt(A)");
    exp(A).write("exp(A)");
    cos(A).write("cos(A)");
    cosh(A).write("cosh(A)");
    sin(A).write("sin(A)");
    sinh(A).write("sinh(A)");
    tan(A).write("tan(A)");
    tanh(A).write("tanh(A)");
    cout << endl;

    (sqrt(A)*sqrt(A)).write("sqrt(A)*sqrt(A)");
    exp(log(A)).write("exp(log(A))");
    cos(acos(A)).write("cos(acos(A))");
    cosh(acosh(A)).write("cosh(acosh(A))");
    sin(asin(A)).write("sin(asin(A))");
    sinh(asinh(A)).write("sinh(asinh(A))");
    tan(atan(A)).write("tan(atan(A))");
    tanh(atanh(A)).write("tanh(atanh(A))");
    cout << endl;
    cout.precision(prec);
  }

  template< class Multivector_T >
  void do_test11()
  {
    typedef Multivector_T m_;
    typedef typename m_::scalar_t s_;
    const s_ pi = s_(l_pi);
    const m_ i("{-1}");
    transcendtest(m_(0));
    transcendtest(m_(1));
    transcendtest(m_(2));
    transcendtest(m_(s_(0.5)));
    transcendtest(m_(pi));
    // exp(i*pi) should be -1
    transcendtest(pi*m_("{1,2}"));
    transcendtest(pi*i/s_(2.0));
    transcendtest(-pi*i*s_(100.0));
    transcendtest(m_("{1}"));
    transcendtest(m_("{-1}+{1}"));
    // the following should produce NaN values
    transcendtest(m_("{-1}+{1}")*s_(1.0e500L));
  }
}

int test11()
{
  cout << "Programming example 11 : Square root and transcendental functions" << endl;
  using namespace peg11;
  cout << endl;
  cout << "framed_multi<double>" << endl;
  do_test11< framed_multi<double> >();
  cout << "matrix_multi<double>" << endl;
  do_test11< matrix_multi<double> >();
  return 0;
}
#endif
