/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg10.cpp : Driver for test 10
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
#include "test10/peg10.h"

int test10()
{
  using namespace peg10;
  cout << "Programming example 10 : Norms and involutions" << endl;
  string a_str = "1{-3,-2}+10{-3,-1}+100{-2,-1}+1000{-3,-2,-1}";
  string b_str = "41{-3,-2}+43{-3,-1}+47{-2,-1}+53{-3,-2,-1}";
  cout << endl;
  cout << "framed_multi<double>" << endl;
  const framed_multi<double> a(a_str);
  const framed_multi<double> b(b_str);
  test_norms_and_involutions< framed_multi<double> >(a, b);
  cout << "matrix_multi<double>" << endl;
  const framed_multi<double> A(a_str);
  const framed_multi<double> B(b_str);
  test_norms_and_involutions< matrix_multi<double> >(A, B);
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return try_catch(test10); }
