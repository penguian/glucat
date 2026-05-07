/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    peg18.cpp : Driver for test 18
                             -------------------
    begin                : Wed May 06 2026
    copyright            : (C) 2001-2026 by Paul C. Leopardi
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
