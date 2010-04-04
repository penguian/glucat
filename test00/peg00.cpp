/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg00.cpp : Driver for test 00
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

#include "test/driver.h"
#include <boost/numeric/ublas/io.hpp>
#include "test00/peg00.h"

int test00()
{
  using namespace peg00;
  std::cout << "Programming example 0 : Geometric algebra identities" << std::endl;

  std::cout << std::endl << "framed_multi<float,-3,4>" << std::endl;
  do_test00< framed_multi<float,-3,4> >(3);

  std::cout << std::endl << "matrix_multi<float,-3,4>" << std::endl;
  do_test00< matrix_multi<float,-3,4> >(3);

  std::cout << std::endl << "framed_multi<float,-4,4>" << std::endl;
  do_test00< framed_multi<float,-4,4> >(4);

  std::cout << std::endl << "matrix_multi<float,-4,4>" << std::endl;
  do_test00< matrix_multi<float,-4,4> >(4);

  std::cout << std::endl << "framed_multi<float,-5,6>" << std::endl;
  do_test00< framed_multi<float,-5,6> >(5);

  std::cout << std::endl << "matrix_multi<float,-5,6>" << std::endl;
  do_test00< matrix_multi<float,-5,6> >(5);

  std::cout << std::endl << "framed_multi<double>" << std::endl;
  do_test00< framed_multi<double> >(5);

  std::cout << std::endl << "matrix_multi<double>" << std::endl;
  do_test00< matrix_multi<double> >(5);

  std::cout << std::endl << "framed_multi<long double>" << std::endl;
  do_test00< framed_multi<long double> >(5);

  std::cout << std::endl << "matrix_multi<long double>" << std::endl;
  do_test00< matrix_multi<long double> >(5);
#ifdef _GLUCAT_USE_QD
  unsigned int old_fpu_control;
  fpu_fix_start(&old_fpu_control);
  std::cout << std::endl << "framed_multi<dd_real>" << std::endl;
  do_test00< framed_multi<dd_real> >(5);

  std::cout << std::endl << "matrix_multi<dd_real>" << std::endl;
  do_test00< matrix_multi<dd_real> >(5);

  std::cout << std::endl << "framed_multi<qd_real>" << std::endl;
  do_test00< framed_multi<qd_real> >(5);

  std::cout << std::endl << "matrix_multi<qd_real>" << std::endl;
  do_test00< matrix_multi<qd_real> >(5);
  fpu_fix_end(&old_fpu_control);
#endif
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return try_catch(test00); }
