/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg14.cpp : Driver for test 14
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

#include "test/driver.h"
#include "test14/peg14.h"

int test14()
{
  using namespace glucat;
  using namespace std;
  cout << "Programming example 14 : framed_multi <-> matrix_multi" << endl;
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

  framed_multi<double> d( 2.0, sub );
  cout << "d " << ((d == 2.0) ? "==" : "!=") << " 2.0; ";
  cout << "d " << ((d != 2.0) ? "!=" : "==") << " 2.0; ";
  cout << "2.0 " << ((2.0 != d) ? "!=" : "==") << " d" << endl;
  cout << "d " << ((d == 4.0) ? "==" : "!=") << " 4.0; ";
  cout << "d " << ((d != 4.0) ? "!=" : "==") << " 4.0; ";
  cout << "4.0 " << ((4.0 != d) ? "!=" : "==") << " d" << endl;

  matrix_multi<double> D( 2.0, sub );
  cout << "D " << ((D == 2.0) ? "==" : "!=") << " 2.0; ";
  cout << "D " << ((D != 2.0) ? "!=" : "==") << " 2.0; ";
  cout << "2.0 " << ((2.0 != D) ? "!=" : "==") << " D" << endl;
  cout << "D " << ((D == 4.0) ? "==" : "!=") << " 4.0; ";
  cout << "D " << ((D != 4.0) ? "!=" : "==") << " 4.0; ";
  cout << "4.0 " << ((4.0 != D) ? "!=" : "==") << " D" << endl;

  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return control_t::control(argc, argv).call(test14); }
