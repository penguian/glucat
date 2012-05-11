/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    PyClical.h : C++ definitions needed by PyClical
                             -------------------
    copyright            : (C) 2008-2012 by Paul C. Leopardi
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
     See also Arvind Raja's original header comments in glucat/glucat.h
 ***************************************************************************/
// References for algorithms:
// [DL]:
// C. Doran and A. Lasenby, "Geometric algebra for physicists", Cambridge, 2003.

#include "glucat/glucat.h"
typedef glucat::tuning<> Tune_P;
#include "glucat/glucat_imp.h"
#include <iostream>
#include <sstream>
#include <iomanip>

/// Create a PyFloatObject object from Scalar_T v.
/// Needed because Scalar_T might not be the same as double.
template <typename Scalar_T>
inline PyObject* PyFloat_FromDouble(Scalar_T v)
{ return ::PyFloat_FromDouble(glucat::numeric_traits<Scalar_T>::to_double(v)); }


// String representations for use by PyClical Python classes.

typedef std::string String;

/// The “official” string representation of Index_Set_T ist.
template<typename Index_Set_T>
inline String index_set_to_repr(const Index_Set_T& ist)
{
  std::ostringstream os;
  os << "index_set(\"" << ist << "\")";
  return os.str();
}

/// The "informal" string representation of Index_Set_T ist.
template<typename Index_Set_T>
inline String index_set_to_str(const Index_Set_T& ist)
{
  std::ostringstream os;
  os << ist;
  return os.str();
}

/// The “official” string representation of Multivector_T mv.
template<typename Multivector_T>
inline String clifford_to_repr(const Multivector_T& mv)
{
  typedef typename Multivector_T::scalar_t scalar_t;
  std::ostringstream os;
  os << std::setprecision(std::numeric_limits<scalar_t>::digits10 + 1);
  os << "clifford(\"" << mv << "\")";
  return os.str();
}

/// The "informal" string representation of Multivector_T mv.
template<typename Multivector_T>
inline String clifford_to_str(const Multivector_T& mv)
{
  typedef typename Multivector_T::scalar_t scalar_t;
  std::ostringstream os;
  if (abs(mv) < std::numeric_limits<scalar_t>::epsilon())
    os << 0.0;
  else
    os << std::setprecision(4) << mv.truncated(scalar_t(1.0e-4));
  return os.str();
}


/// Definitions for 3D Conformal Geometric Algebra [DL].
namespace cga3
{
  /// Convert Euclidean 3D vector to Conformal Geometric Algebra null vector [DL (10.50)].
  template<typename Multivector_T>
  inline Multivector_T cga3(const Multivector_T& x)
  {
    typedef Multivector_T cl;
    typedef typename cl::index_set_t ist;
    static const cl ninf3 = cl(ist(4)) + cl(ist(-1));

    return (cl(ist(4)) - x) * ninf3 * (x - cl(ist(4)));
  }

  /// Convert CGA3 null vector to standard Conformal Geometric Algebra  null vector [DL (10.52)].
  template<typename Multivector_T>
  inline Multivector_T cga3std(const Multivector_T& X)
  {
    typedef Multivector_T cl;
    typedef typename cl::index_set_t ist;
    typedef typename cl::scalar_t scalar_t;
    static const cl ninf3 = cl(ist(4)) + cl(ist(-1));

    return scalar_t(-2.0) * X / (X & ninf3);
  }

  /// Convert CGA3 null vector to Euclidean 3D vector [DL (10.50)].
  template<typename Multivector_T>
  inline Multivector_T agc3(const Multivector_T& X)
  {
    typedef Multivector_T cl;
    typedef typename cl::index_set_t ist;
    typedef typename cl::scalar_t scalar_t;

    const cl& cga3stdX = cga3std(X)(1);
    return (cl(ist(1))*cga3stdX[ist(1)] + cl(ist(2))*cga3stdX[ist(2)] + cl(ist(3))*cga3stdX[ist(3)]) / scalar_t(2.0);
  }
}


// Specifications of the IndexSet and Clifford C++ classes for use with PyClical.

using namespace glucat;
const index_t lo_ndx = DEFAULT_LO;
const index_t hi_ndx = DEFAULT_HI;
typedef index_set<lo_ndx,hi_ndx> IndexSet;

typedef double scalar_t;
typedef matrix_multi<scalar_t> Clifford;
