#ifndef _GLUCAT_GENERATION_IMP_H
#define _GLUCAT_GENERATION_IMP_H
/***************************************************************************
	  GluCat : Generic library of universal Clifford algebra templates
    generation_imp.h : Implement functions for generation of the matrix representation
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

#include "glucat/global.h"
#include "glucat/generation.h"
#include "glucat/matrix.h"

namespace glucat { namespace gen
{
  // References for algorithms:
  // [M]: Scott Meyers, "Effective C++" Second Edition, Addison-Wesley, 1998.
  // [P]: Ian R. Porteous, "Clifford algebras and the classical groups", Cambridge UP, 1995.
  // [L]: Pertti Lounesto, "Clifford algebras and spinors", Cambridge UP, 1997.

  /// Single instance of generator table
  // Reference: [M] Item 47
  template< class Matrix_T >
  auto
  generator_table<Matrix_T>::
  generator() -> generator_table<Matrix_T>&
  { static generator_table<Matrix_T> g; return g;}

  /// Pointer to generators for a specific signature
  // Reference: [P] Table 15.27, p 133
  template< class Matrix_T >
  inline
  auto
  generator_table<Matrix_T>::
  operator() (const index_t p, const index_t q) -> const Matrix_T*
  {
    const auto bott = pos_mod(p-q, 8);
    switch(bott)
    {
    case 0:
    case 2:
      // Construct generators
      return &(gen_vector(p, q)[q]);
    default:
      // Select generators from the vector for a larger frame
      const auto super_p = p + std::max(offset_to_super[bott],index_t(0));
      const auto super_q = q - std::min(offset_to_super[bott],index_t(0));
      return &(gen_vector(super_p, super_q)[super_q]);
    }
  }

  /// Construct a vector of generators for a specific signature
  template< class Matrix_T >
  auto
  generator_table<Matrix_T>::
  gen_vector(const index_t p, const index_t q) -> const std::vector<Matrix_T>&
  {
    using result_t = std::vector<Matrix_T>;
    const auto card = p + q;
    const auto bias = p - q;
    const auto bott = pos_mod(bias, 8);
    const auto sig = signature_t(p, q);
    if (this->find(sig) == this->end())
      switch(bott)
      {
      case 0:
        if (bias < 0)
          // Construct generators for p,q given generators for p+4,q-4
          gen_from_pp4_qm4(gen_vector(p+4, q-4), sig);
        else if (bias > 0)
          // Construct generators for p,q given generators for p-4,q+4
          gen_from_pm4_qp4(gen_vector(p-4, q+4), sig);
        else if (card == 0)
        { // Base case. Save a generator vector containing one matrix, size 1.
          auto result = result_t(1, matrix::unit<Matrix_T>(1));
          this->insert(make_pair(sig, result));
        }
        else
          // Construct generators for p,q given generators for p-1,q-1
          gen_from_pm1_qm1(gen_vector(p-1, q-1), sig);
        break;
      case 2:
        if (bias < 2)
          // Construct generators for p,q given generators for p+4,q-4
          gen_from_pp4_qm4(gen_vector(p+4, q-4), sig);
        else if (bias > 2)
          // Construct generators for p,q given generators for p-4,q+4
          gen_from_pm4_qp4(gen_vector(p-4, q+4), sig);
        else
          // Construct generators for p,q given generators for q+1,p-1
          gen_from_qp1_pm1(gen_vector(q+1, p-1), sig);
        break;
      default:
        break;
      }
    return (*this)[sig];
  }

  /// Construct generators for p,q given generators for p-1,q-1
  // Reference: [P] Proposition 15.17, p 131
  template< class Matrix_T >
  void
  generator_table<Matrix_T>::
  gen_from_pm1_qm1(const std::vector<Matrix_T>& old, const signature_t sig)
  {
    const auto new_size = old.size() + 2;
    using size_t = decltype(new_size);
    using result_t = std::vector<Matrix_T>;
    auto result = result_t(new_size);

    const auto old_dim = old[0].size1();
    const auto& eye = matrix::unit<Matrix_T>(old_dim);

    auto neg = Matrix_T(2,2,2);
    neg(0,1) =    -1;
    neg(1,0) = 1;

    auto pos = neg;
    pos(0,1) =     1;

    auto dup = Matrix_T(2,2,2);
    dup(0,0) = 1;
    dup(1,1) =    -1;

    result[0] = matrix::mono_kron(neg, eye);
    for (auto
        k = size_t(1);
        k != new_size-1;
        ++k)
      result[k] = matrix::mono_kron(dup, old[k-1]);
    result[new_size-1] = matrix::mono_kron(pos, eye);

    // Save the resulting generator array.
    this->insert(make_pair(sig, result));
  }

  /// Construct generators for p,q given generators for p-4,q+4
  // Reference: [L] 16.4 Periodicity of 8, p216
  template< class Matrix_T >
  void
  generator_table<Matrix_T>::
  gen_from_pm4_qp4(const std::vector<Matrix_T>& old, const signature_t sig)
  {
    const auto old_size = old.size();
    using size_t = decltype(old_size);
    using result_t = std::vector<Matrix_T>;
    auto result = result_t(old_size);

    auto h = old[0];
    for (auto
        k = size_t(1);
        k != size_t(4);
        ++k)
      h = matrix::mono_prod(old[k], h);

    for (auto
        k = size_t(0);
        k != old_size-4;
        ++k)
      result[k] = old[k+4];
    for (auto
        k = old_size-4;
        k != old_size;
        ++k)
      result[k] = matrix::mono_prod(old[k+4-old_size], h);
    // Save the resulting generator array.
    this->insert(make_pair(sig, result));
  }

  /// Construct generators for p,q given generators for p+4,q-4
  // Reference: [L] 16.4 Periodicity of 8, p216
  template< class Matrix_T >
  void
  generator_table<Matrix_T>::
  gen_from_pp4_qm4(const std::vector<Matrix_T>& old, const signature_t sig)
  {
    const auto old_size = old.size();
    using size_t = decltype(old_size);
    using result_t = std::vector<Matrix_T>;
    auto result = result_t(old_size);

    auto h = old[old_size-1];
    for (auto
        k = size_t(1);
        k != size_t(4);
        ++k)
      h = matrix::mono_prod(old[old_size-1-k], h);

    for (auto
        k = size_t(0);
        k != size_t(4);
        ++k)
      result[k] = matrix::mono_prod(old[k+old_size-4], h);
    for (auto
        k = size_t(4);
        k != old_size;
        ++k)
      result[k] = old[k-4];
    // Save the resulting generator array.
    this->insert(make_pair(sig, result));
  }

  /// Construct generators for p,q given generators for q+1,p-1
  // Reference: [P] Proposition 15.20, p 131
  template< class Matrix_T >
  void
  generator_table<Matrix_T>::
  gen_from_qp1_pm1(const std::vector<Matrix_T>& old, const signature_t sig)
  {
    const auto old_size = old.size();
    using size_t = decltype(old_size);
    using result_t = std::vector<Matrix_T>;
    auto result = result_t(old_size);

    const auto& h = old[old_size-1];
    for (auto
        k = size_t(0);
        k != old_size-1;
        ++k)
      result[k] = matrix::mono_prod(old[old_size-2-k], h);
    result[old_size-1] = h;

    // Save the resulting generator array.
    this->insert(make_pair(sig, result));
  }

} }
#endif  // _GLUCAT_GENERATION_IMP_H
