/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    versor.cpp : Timing benchmark for versor products
                             -------------------
    begin                : Thu Jun 11 2026
    copyright            : (C) 2001-2026 by Paul C. Leopardi
 ***************************************************************************/

#include "test/driver.h"
#include "versor/versor.h"

int versor_test(const int n)
{
  using namespace glucat_versor_test;
  if (n > max_n)
  {
    cout << "Value " << n << " is too big." << endl;
    cout << "Maximum value allowed is " << max_n << "." << endl;
    return 1;
  }
  cout << "Benchmarking versor product performance up to Cl(" << n << ",0)" << endl;
  run_benchmark<double>(n);
  return 0;
}

using namespace glucat;

int main(int argc, char ** argv)
{
  using namespace std;
  for (argc--, argv++; argc != 0; argc--, argv++)
  {
    int max_index = 0;
    sscanf(*argv, "%d", &max_index);
    if (max_index > 0)
      try_catch( versor_test, max_index );
  }
  return 0;
}
