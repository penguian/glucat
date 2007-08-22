/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg03.cpp : Driver for test 03
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
#include "test03/peg03.h"

int test03()
{
  using namespace peg03;
  cout << "Programming example 3 : Powers, output to file" << endl;
  cout << endl;
  cout << "In Cl_{2,1} let x = 2e1 + 6.2e2 - 4.1e23 and" << endl;
  cout << "                y = 7e1 - 2.1e2 + 9.6e12 + 6e23." << endl;
  cout << "Compute 2x + 3y, x to the power 5, y to the outer power 3." << endl;
  cout << "Find the even part of x*y and write the result to a file called eg3.res" << endl;
  cout << endl;
  cout << "framed_multi<double>" << endl;
  do_test3< framed_multi<double> >();
  cout << "matrix_multi<double>" << endl;
  do_test3< matrix_multi<double> >();
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return try_catch(test03); }
