#ifndef _GLUCAT_GENERATION_IMP_H
#define _GLUCAT_GENERATION_IMP_H
/***************************************************************************
	  GluCat : Generic library of universal Clifford algebra templates
    generation_imp.h : Implement functions for generation of the matrix representation
                             -------------------
    begin                : Wed Jan 23 2002
    copyright            : (C) 2002 by Paul C. Leopardi
    email                : leopardi@bigpond.net.au
 ***************************************************************************
 *   This library is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU Lesser General Public License as        *
 *   published by the Free Software Foundation; either version 2.1 of the  *
 *   License, or (at your option) any later version.                       *
 *   See http://www.fsf.org/copyleft/lesser.html for details               *
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
  // References for algorithms:
  // [P]:
  // Ian R. Porteous, "Clifford algebras and the classical groups", Cambridge UP, 1995.
  // [L}:
  // Pertti Lounesto, "Clifford algebras and spinors", Cambridge UP, 1997.

  namespace generation
  {
    /// Generate the next generation of generators, given the previous generation
    // Reference: [P] Proposition 15.17, p 131
    template< class Matrix_T >
    void
    gen_from_pm1_qm1(const vector<Matrix_T>& old, const signature_t sig,
                    map< signature_t, vector<Matrix_T> >& generators)
    {
      Matrix_T pos(2,2);
      mtl::set(pos, 0);
      pos(0,1) =     1;
      pos(1,0) = 1;

      Matrix_T dup(2,2);
      mtl::set(dup, 0);
      dup(0,0) = 1;
      dup(1,1) =    -1;

      Matrix_T neg(2,2);
      mtl::set(neg, 0);
      neg(0,1) =    -1;
      neg(1,0) = 1;

      const int new_size = old.size() + 2;
      const int old_dim = old[0].nrows();
  		Matrix_T eye(old_dim, old_dim);
      unit(old_dim, eye);

      vector<Matrix_T> result(new_size);
      kron(neg , eye, result[0]);
      for (int k = 1; k != new_size-1; ++k)
        kron(dup , old[k-1], result[k]);
      kron(pos , eye, result[new_size-1]);

      // Save the resulting generator array.
      generators.insert(make_pair(sig, result));
    }

    /// Reference: [L] 16.4 Periodicity of 8, p216
    template< class Matrix_T >
    void
    gen_from_pm4_qp4(const vector<Matrix_T>& old, const signature_t sig,
                    map< signature_t, vector<Matrix_T> >& generators)
    {
      const int dim = old[0].nrows();
      const int old_size = old.size();
      Matrix_T h(dim, dim);
      mtl::copy(old[0], h);
      for (int k = 1; k != 4; ++k)
      {
        Matrix_T temp(dim, dim);
        mtl::mult(old[k], h, temp);
        mtl::copy(temp, h);
      }
      vector<Matrix_T> result(old_size);
      for (int k = 0; k != 4; ++k)
      {
        result[k+old_size-4] = Matrix_T(dim, dim);
        mtl::mult(old[k], h, result[k+old_size-4]);
      }
      for (int k = 4; k != old_size; ++k)
      {
        result[k-4] = Matrix_T(dim, dim);
        mtl::copy(old[k], result[k-4]);
      }

      // Save the resulting generator array.
      generators.insert(make_pair(sig, result));
    }

    /// Reference: [L] 16.4 Periodicity of 8, p216
    template< class Matrix_T >
    void
    gen_from_pp4_qm4(const vector<Matrix_T>& old, const signature_t sig,
                    map< signature_t, vector<Matrix_T> >& generators)
    {
      const int dim = old[0].nrows();
      const int old_size = old.size();
      Matrix_T h(dim, dim);
      mtl::copy(old[old_size-1], h);
      for (int k = 1; k != 4; ++k)
      {
        Matrix_T temp(dim, dim);
        mtl::mult(old[old_size-1-k], h, temp);
        mtl::copy(temp, h);
      }
      vector<Matrix_T> result(old_size);
      for (int k = 0; k != 4; ++k)
      {
        result[k] = Matrix_T(dim, dim);
        mtl::mult(old[k+old_size-4], h, result[k]);
      }
      for (int k = 4; k != old_size; ++k)
      {
        result[k] = Matrix_T(dim, dim);
        mtl::copy(old[k-4], result[k]);
      }

      // Save the resulting generator array.
      generators.insert(make_pair(sig, result));
    }

    /// Reference: [P] Proposition 15.20, p 131
    template< class Matrix_T >
    void
    gen_from_qp1_pm1(const vector<Matrix_T>& old, const signature_t sig,
                    map< signature_t, vector<Matrix_T> >& generators)
    {
      const int dim = old[0].nrows();
      const int old_size = old.size();
      const Matrix_T& a = old[old_size-1];
      vector<Matrix_T> result(old_size);
      int m = 0;
      for (int k = old_size-1; k != 0; --k, ++m)
      {
        result[m] = Matrix_T(dim, dim);
        mtl::mult(old[k-1], a, result[m]);
      }
      result[old_size-1] = Matrix_T(dim, dim);
      mtl::copy(a, result[old_size-1]);

      // Save the resulting generator array.
      generators.insert(make_pair(sig, result));
    }

    /// The vector of generators for a specific signature
    // Generate a specific generation of generators, given the whole family tree
    template< class Matrix_T >
    const vector<Matrix_T>&
    gen_vector(const index_t p, const index_t q, map< signature_t, vector<Matrix_T> >& generators)
    {
      const index_t card = p + q;
      const index_t bias = p - q;
      const index_t bott = pos_mod(bias, 8);
      const signature_t sig(p, q);
      if (generators.find(sig) == generators.end())
        switch(bott)
        {
        case 0:
          if (bias < 0)
            // Reference: [L] 16.4 Periodicity of 8, p216
            gen_from_pp4_qm4(gen_vector(p+4, q-4, generators), sig, generators);
          else if (bias > 0)
            // Reference: [L] 16.4 Periodicity of 8, p216
            gen_from_pm4_qp4(gen_vector(p-4, q+4, generators), sig, generators);
          else if (card == 0)
          { // Base case. Generator array will contain one matrix, size 1.
            // Save this generator array.
            vector<Matrix_T> result(1);
            result[0] = Matrix_T(1,1);
            mtl::set(result[0], 0);
            generators.insert(make_pair(sig, result));
          }
          else
            // Use recursion to define generator array in terms of smaller generator array
            // Reference: [P] Proposition 15.17, p 131
            gen_from_pm1_qm1(gen_vector(p-1, q-1, generators), sig, generators);
          break;
        case 2:
          if (bias < 2)
            // Reference: [L] 16.4 Periodicity of 8, p216
            gen_from_pp4_qm4(gen_vector(p+4, q-4, generators), sig, generators);
          else if (bias > 2)
            // Reference: [L] 16.4 Periodicity of 8, p216
            gen_from_pm4_qp4(gen_vector(p-4, q+4, generators), sig, generators);
          else
            // Reference: [P] Proposition 15.20, p 131
            gen_from_qp1_pm1(gen_vector(q+1, p-1, generators), sig, generators);
          break;
        default:
          break;
        }
      return generators[sig];
    }
  }

  /// Pointer to generators for a specific signature
  template< class Matrix_T >
  const Matrix_T*
  gen(const index_t p, const index_t q, map< signature_t, vector<Matrix_T> >& generators)
  {
    // Reference: [P] Table 15.27, p 133
    static const int offset_to_super[] = {0, -1, 0, -1, -2, 3, 2, 1};
    const index_t bott = pos_mod(p-q, 8);
    switch(bott)
    {
    case 0:
      // Reference: [P] Proposition 15.17, p 131
    case 2:
      // Reference: [P] Proposition 15.20, p 131
      return &(generation::gen_vector(p, q, generators)[q]);
      break;
    default:
      // Select generators from the vector for a larger frame
      // Reference: [P] Table 15.27, p 133
      const index_t super_p = p + max(offset_to_super[bott],0);
      const index_t super_q = q - min(offset_to_super[bott],0);
      return &(generation::gen_vector(super_p, super_q, generators)[super_q]);
      break;
    }
  }

  /// Determine the matrix dimension of the fold of a subalegbra
  template< typename Matrix_Index_T, const index_t LO, const index_t HI >
  const Matrix_Index_T
  folded_dim( const index_set<LO,HI>& sub )
  {
    // Reference: [P] Table 15.27, p 133
    static const int offset_log2_dim[] = {0, 1 , 0, 1, 1, 2, 1, 1};
    const index_set<LO,HI> folded = sub.fold();
    const int p = max( int(folded.max()), 0);
    const int q = max(-int(folded.min()), 0);
    const int bott = pos_mod(p-q, 8);
    return 1 << ((p+q)/2 + offset_log2_dim[bott]);
  }
}
#endif  // _GLUCAT_GENERATION_IMP_H
