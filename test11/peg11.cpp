/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg11.cpp : Driver for test 11
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

#if !defined(NDEBUG)
#define NDEBUG
#endif
#include "test/driver.h"
#include "test11/peg11.h"

int test11()
{
  using namespace peg11;
  cout <<   "Programming example 11 : Square root and transcendental functions" << endl;
  cout << endl;
  do_test11_tuned<float,DEFAULT_LO,DEFAULT_HI,default_tuning_promoted_p>("Promoted:", "float");
  do_test11_tuned<double,DEFAULT_LO,DEFAULT_HI,default_tuning_promoted_p>("Promoted:", "double");
  do_test11_tuned<long double,DEFAULT_LO,DEFAULT_HI,default_tuning_promoted_p>("Promoted:", "long double");
#  if defined(_GLUCAT_USE_QD)
  do_test11_tuned<dd_real,DEFAULT_LO,DEFAULT_HI>("Default:", "dd_real");
  do_test11_tuned<dd_real,DEFAULT_LO,DEFAULT_HI,default_tuning_promoted_p>("Promoted:", "dd_real");
  do_test11_tuned<qd_real,DEFAULT_LO,DEFAULT_HI,default_tuning_promoted_p>("Promoted:", "qd_real");
#  endif
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{ return control_t::control(argc, argv).call(test11); }
