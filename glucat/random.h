#ifndef _GLUCAT_RANDOM_H
#define _GLUCAT_RANDOM_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    random.h : Random number generator with single instance per Scalar_T
                             -------------------
    begin                : 2010-03-28
    copyright            : (C) 2001-2012 by Paul C. Leopardi
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
 "Clifford algebras with numeric and symbolic computations, Birkhauser, 1996."
 ***************************************************************************
 See also Arvind Raja's original header comments and references in glucat.h
 ***************************************************************************/

#if   defined(_GLUCAT_USE_GSL_RANDOM)
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
#else
# include <random>
#endif

namespace glucat
{
  /// Random number generator with single instance per Scalar_T
  // Enforce singleton
  // Reference: A. Alexandrescu, "Modern C++ Design", Chapter 6
  template< typename Scalar_T >
  class random_generator
  {
  private:
    /// Friend declaration to avoid compiler warning:
    /// "... only defines a private destructor and has no friends"
    /// Ref: Carlos O'Ryan, ACE http://doc.ece.uci.edu
    friend class friend_for_private_destructor;
  public:
    /// Single instance of Random number generator
    static random_generator& generator() { static random_generator g; return g;}
    random_generator(const random_generator&) = delete;
    random_generator& operator= (const random_generator&) = delete;
  private:
    static const unsigned long seed = 19590921UL;
#if defined(_GLUCAT_USE_GSL_RANDOM)

    gsl_rng* gen;

    random_generator() :
    gen(gsl_rng_alloc(gsl_rng_mt19937))
    { gsl_rng_set(this->gen, seed); }

    ~random_generator()
    { gsl_rng_free(this->gen); }

  public:
    Scalar_T uniform()
    { return Scalar_T(gsl_ran_flat(this->gen, 0.0, 1.0)); }
    Scalar_T normal()
    { return Scalar_T(gsl_ran_gaussian(this->gen, 1.0)); }

#else

    std::mt19937 uint_gen;
    std::uniform_real_distribution<double> uniform_dist;
    std::normal_distribution<double> normal_dist;

    random_generator() :
    uint_gen(), uniform_dist(0.0, 1.0), normal_dist(0.0, 1.0)
    { this->uint_gen.seed(seed); }

    ~random_generator() = default;

  public:
    Scalar_T uniform()
    { return Scalar_T(this->uniform_dist(this->uint_gen)); }
    Scalar_T normal()
    { return Scalar_T(this->normal_dist(this->uint_gen)); }

#endif
  };
}

#endif // _GLUCAT_RANDOM_H
