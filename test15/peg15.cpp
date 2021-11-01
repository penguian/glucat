/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg15.cpp : Driver for test 15
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
#include "test15/peg15.h"

int test15()
{
  using namespace glucat;
  using namespace std;
  cout << "Programming example 15 : Subscripts and parts" << endl;
  using cm = matrix_multi<double>;
  using cf = framed_multi<double>;
  using e = cm::index_set_t;
  cf a("{-3}+{-2}+{-1}+2.0{1,2,3}");
  cf b("1.0+1.e8{-1}+1.e-8{1}");
  cf c = a * b;
  e sub = a.frame() | b.frame();
  cm A( a, sub );
  cm B( b, sub );
  cm C = A * B;
  cout << "a: " << a << endl;
  cout << "A: " << A << endl;
  cout << "b: " << b << endl;
  cout << "B: " << B << endl;
  cout << "c = a*b: " << c << endl;
  cout << "C = A*B: " << C << endl;
  for (index_t k = -3; k != 4; ++k)
  {
    cout << "c[{" << k << "}]: " << c[e(k)] << endl;
    cout << "C[{" << k << "}]: " << C[e(k)] << endl;
  }
  for (index_t k = 0; k != 6; ++k)
  {
    cout << "c(" << k << "): " << c(k) << endl;
    cout << "C(" << k << "): " << C(k) << endl;
  }
  cout << "even(c): "   << even(c) << endl;
  cout << "scalar(c): " << scalar(c) << endl;
  cout << "pure(c): "   << pure(c) << endl;
  cout << "real(c): "   << real(c) << endl;
  cout << "imag(c): "   << imag(c) << endl;
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return control_t::control(argc, argv).call(test15); }
