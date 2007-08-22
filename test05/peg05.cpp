/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg05.cpp : Driver for test 05
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
#include "test05/peg05.h"

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

using namespace glucat;

int main(int argc, char ** argv)
{ return try_catch(test05); }
