#ifndef GLUCAT_TEST_PEG06_H
#define GLUCAT_TEST_PEG06_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg06.cpp : programming example 6 : outer exponential
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2007 by Paul C. Leopardi
 ***************************************************************************

    This library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library.  If not, see <http://www.gnu.org/licenses/>.

 ***************************************************************************
 This library is based on a prototype written by Arvind Raja and was
 licensed under the LGPL with permission of the author. See Arvind Raja,
 "Object-oriented implementations of Clifford algebras in C++: a prototype",
 in Ablamowicz, Lounesto and Parra (eds.)
 "Clifford algebras with numeric and symbolic computations", Birkhauser, 1996.
 ***************************************************************************
     See also Arvind Raja's original header comments in glucat.h
 ***************************************************************************/

namespace peg06
{
  using namespace glucat;
  using namespace std;

  template< class Multivector_T >
  static
  Multivector_T 
  expo(const Multivector_T& f)
  {
    typedef Multivector_T number;
    const number ff = f^f;
    return 1.0 + f + ff/2.0 + (ff^(f/6.0 + ff/24.0));
  }

  template< class Multivector_T >
  static
  void 
  do_test6()
  {
    typedef Multivector_T number;

    const number q1("{-1}"), q2("{-2}"), q3("{-3}"), q4("{-4}"), q5("{-5}"), q6("{-6}");
    const number B = q1*q2 + q3*q4 + q5*q6;
    B.write("B = ");
    expo(B).write("outer exponential of B = ");
    cout << endl;
  }
}

int test06();

#endif
