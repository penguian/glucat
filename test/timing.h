#ifndef GLUCAT_TEST_TIMING_H
#define GLUCAT_TEST_TIMING_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    timing.h : Common definitions for timing tests
                             -------------------
    begin                : Tue 2012-03-27
    copyright            : (C) 2012 by Paul C. Leopardi
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

namespace glucat
{
  namespace timing
  {
    /// Timing constant: milliseconds per second
    const double MS_PER_SEC = 1000.0;

    /// Timing constant: milliseconds per clock
    const double MS_PER_CLOCK = MS_PER_SEC / double(CLOCKS_PER_SEC);

    /// Timing constant: trial expansion factor
    const int EXTRA_TRIALS = 2;

    /// Elapsed time in milliseconds
    inline
    static
    double
    elapsed(clock_t cpu_time)
    { return double(clock() - cpu_time) * MS_PER_CLOCK; }

  }
}
#endif // GLUCAT_TEST_TIMING_H
