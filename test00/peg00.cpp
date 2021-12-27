/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg00.cpp : Driver for test 00
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2021 by Paul C. Leopardi
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
#include "test00/peg00.h"

int test00()
{
  using namespace peg00;
  std::cout << "Programming example 0 : Geometric algebra identities" << std::endl;

  do_test00_all_tunes<float,-3,4>("float", 3);
  do_test00_all_tunes<float,-4,4>("float", 4);
  do_test00_all_tunes<float,-5,6>("float", 5);
  do_test00_all_tunes<double>("double", 5);
  do_test00_all_tunes<long double>("long double", 5);
#ifdef _GLUCAT_USE_QD
  unsigned int old_fpu_control;
  fpu_fix_start(&old_fpu_control);
  do_test00_all_tunes<dd_real>("dd_real", 5);
  do_test00_all_tunes<qd_real>("qd_real", 5);
  fpu_fix_end(&old_fpu_control);
#endif
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return control_t::control(argc, argv).call(test00); }
