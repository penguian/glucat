/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg11.cpp : Driver for test 11
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2010 by Paul C. Leopardi
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

#if !defined(NDEBUG)
#define NDEBUG
#endif
#include "test/driver.h"
#include "test11/peg11.h"

int test11()
{
  using namespace peg11;
  cout << "Programming example 11 : Square root and transcendental functions" << endl;
  cout << endl;
  cout << "framed_multi<float>" << endl;
  do_test11< framed_multi<float> >();
  cout << "matrix_multi<float>" << endl;
  do_test11< matrix_multi<float> >();
  cout << "framed_multi<double>" << endl;
  do_test11< framed_multi<double> >();
  cout << "matrix_multi<double>" << endl;
  do_test11< matrix_multi<double> >();
  cout << "framed_multi<long double>" << endl;
  do_test11< framed_multi<long double> >();
  cout << "matrix_multi<long double>" << endl;
  do_test11< matrix_multi<long double> >();
#ifdef _GLUCAT_USE_QD
  cout << "framed_multi<dd_real>" << endl;
  do_test11< framed_multi<dd_real> >();
  cout << "matrix_multi<dd_real>" << endl;
  do_test11< matrix_multi<dd_real> >();
  cout << "framed_multi<qd_real>" << endl;
  do_test11< framed_multi<qd_real> >();
  cout << "matrix_multi<qd_real>" << endl;
  do_test11< matrix_multi<qd_real> >();
#endif
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return try_catch(test11); }
