/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    products.cpp : Products timing test driver
                             -------------------
    begin                : Sat 2007-09-01
    copyright            : (C) 2007-2021 by Paul C. Leopardi
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
#include "products/products.h"

int products(const int n)
{
  using namespace glucat_prod_test;
  if (n > max_n)
  {
    cout << "Value " << n << " is too big." << endl;
    cout << "Maximum value allowed is " << max_n << "." << endl;
    return 1;
  }
  const int small_n = std::min(n,4);
  cout << endl << "framed_multi<float,-4,4>" << endl;
  prod_test< framed_multi<float,-4,4> >(small_n, max_n);
  cout << endl << "matrix_multi<float,-4,4>" << endl;
  prod_test< matrix_multi<float,-4,4> >(small_n, max_n);
  cout << endl << "framed_multi<double>" << endl;
  prod_test< framed_multi<double> >(n, max_n);
  cout << endl << "matrix_multi<double>" << endl;
  prod_test< matrix_multi<double> >(n, max_n);

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
    try_catch(products, n);
  }
  return 0;
}
