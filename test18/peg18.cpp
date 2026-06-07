/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg18.cpp : Driver for test 18
                             -------------------
    begin                : Wed May 06 2026
    copyright            : (C) 2001-2026 by Paul C. Leopardi
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

 ***************************************************************************/

#include "test/driver.h"
#include "test18/peg18.h"

using namespace glucat;
using namespace std;

int test18()
{
  cout << "Programming example 18 : Consistency across signatures" << endl;

  cout << "double:" << endl;
  consistency_test<double>();

  return 0;
}

int main(int argc, char ** argv)
{ return control_t::control(argc, argv).call(test18); }
