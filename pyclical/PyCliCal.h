#include "glucat/glucat.h"
typedef glucat::tuning<> Tune_P;
#include "glucat/glucat_imp.h"
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace glucat;
using namespace std;

typedef index_set<DEFAULT_LO,DEFAULT_HI> IndexSet;
typedef double scalar_t;
typedef framed_multi<scalar_t> Clifford;

string index_set_to_repr(const IndexSet& ist)
{
  ostringstream os;
  os << ist;
  return os.str();
}

string index_set_to_str(const IndexSet& ist)
{
  ostringstream os;
  os << ist;
  return os.str();
}

string clifford_to_repr(const Clifford& mv)
{
  ostringstream os;
  os << setprecision(17) << mv;
  return os.str();
}

string clifford_to_str(const Clifford& mv)
{
  ostringstream os;
  if (abs(mv) < 1.0e-8)
    os << 0.0;
  else
    os << setprecision(4) << mv.truncated(1.0e-4);
  return os.str();
}
