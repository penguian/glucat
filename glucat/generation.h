#ifndef _GLUCAT_GENERATION_H
#define _GLUCAT_GENERATION_H
/***************************************************************************
	  GluCat : Generic library of universal Clifford algebra templates
    generation.h : Declare functions for generation of the matrix representation
                             -------------------
    begin                : Wed Jan 23 2002
    copyright            : (C) 2002-2012 by Paul C. Leopardi
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

#include <utility>
#include <array>
#include <map>
#include <vector>

namespace glucat { namespace gen
{
  namespace ublas = boost::numeric::ublas;
  
  /// A signature is a pair of indices, p, q, with p == frame.max(), q == -frame.min()
  using signature_t = std::pair<index_t, index_t>;

  /// Table of generators for specific signatures
  template< class Matrix_T >
  class generator_table :
  private std::map< signature_t, std::vector<Matrix_T> >
  {
  public:
    /// Pointer to generators for a specific signature
    const Matrix_T* operator() (const index_t p, const index_t q);
    /// Single instance of generator table
    static generator_table<Matrix_T>& generator();
  private:
    /// Construct a vector of generators for a specific signature
    const std::vector<Matrix_T>& gen_vector(const index_t p, const index_t q);
    /// Construct generators for p,q given generators for p-1,q-1
    void gen_from_pm1_qm1(const std::vector<Matrix_T>& old, const signature_t sig);
    /// Construct generators for p,q given generators for p-4,q+4
    void gen_from_pm4_qp4(const std::vector<Matrix_T>& old, const signature_t sig);
    /// Construct generators for p,q given generators for p+4,q-4
    void gen_from_pp4_qm4(const std::vector<Matrix_T>& old, const signature_t sig);
    /// Construct generators for p,q given generators for q+1,p-1
    void gen_from_qp1_pm1(const std::vector<Matrix_T>& old, const signature_t sig);

    /// Friend declaration to avoid compiler warning:
    /// "... only defines a private destructor and has no friends"
    /// Ref: Carlos O'Ryan, ACE http://doc.ece.uci.edu
    friend class friend_for_private_destructor;
    // Enforce singleton
    // Reference: A. Alexandrescu, "Modern C++ Design", Chapter 6
    generator_table() = default;
    ~generator_table() = default;
  public:
    generator_table(const generator_table&) = delete;
    generator_table& operator= (const generator_table&) = delete;
  };

  /// Offsets between the current signature and that of the real superalgebra
  static const std::array<index_t, 8> offset_to_super = {0,-1, 0,-1,-2, 3, 2, 1};

} }
#endif  // _GLUCAT_GENERATION_H
