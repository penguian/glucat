/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    test19.cpp : Driver for test 19
                             -------------------
    begin                : Thu Jun 11 2026
    copyright            : (C) 2001-2026 by Paul C. Leopardi
 ***************************************************************************/

#include "test/driver.h"
#include "test19/peg19.h"

using namespace glucat;
using namespace std;

int test19()
{
  cout << "Programming example 19 : Sandwich products equivalence" << endl;

  cout << "double:" << endl;
  versor_equivalence_test<double>();

  return 0;
}

int main(int argc, char ** argv)
{ return control_t::control(argc, argv).call(test19); }
