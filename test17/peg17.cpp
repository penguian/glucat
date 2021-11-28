/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg17.cpp : Driver for test 17
                             -------------------
    begin                : Sun 2021-11-24
    copyright            : (C) 2001-2025 by Paul C. Leopardi
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
#include "test17/peg17.h"

using namespace glucat;
using namespace std;

int test17()
{
  cout << "Programming example 17 : Truncation and printing" << endl;

  cout << "float:" << endl;
  printing_test<float>(7);
  cout << "double:" << endl;
  printing_test<double>(14);
  cout << "long double:" << endl;
  printing_test<long double>(17);
  return 0;
}

int main(int argc, char ** argv)
{ return control_t::control(argc, argv).call(test17); }
