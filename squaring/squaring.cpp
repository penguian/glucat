/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    squaring.cpp : Timing test driver
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

#include "squaring/squaring.h"

int squaring(const int n)
{
  using namespace glucat_mult_test;
  if (n > max_n)
  {
    cout << "Value " << n << " is too big." << endl;
    cout << "Maximum value allowed is " << max_n << "." << endl;
    return 1;
  }
  cout << endl << "framed_multi<double>" << endl;
  mult_test< framed_multi<double> >(n, max_n);
  cout << endl << "matrix_multi<double>" << endl;
  mult_test< matrix_multi<double> >(n, max_n);

  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{
  using namespace std;
  for (argc--, argv++; argc != 0; argc--, argv++)
  {
    int n = 0;
    sscanf(*argv, "%d", &n);
    try_catch(squaring, n);
  }
  return 0;
}
