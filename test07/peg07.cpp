/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg07.cpp : Driver for test 07
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
#include "test07/peg07.h"

int test07()
{
  using namespace peg07;
  cout << "Programming example 7 : Triality" << endl;
  cout << endl;
  cout << "framed_multi<double>" << endl;
  do_test7< framed_multi<double> >();
  cout << "matrix_multi<double>" << endl;
  do_test7< matrix_multi<double> >();
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return try_catch(test07); }
