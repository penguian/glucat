#ifndef _GLUCAT_UBLAS_IMP_H
#define _GLUCAT_UBLAS_IMP_H
#ifndef _GLUCAT_HAVE_UBLAS_LU_H

/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    ublas_imp.h : Implement uBLAS interface
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001 by Paul C. Leopardi
                         : uBLAS interface contributed by Joerg Walter
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
namespace boost { namespace numeric {
  namespace ublas
  {
    /// Swap rows of a vector
    template<class PM, class MV>
    BOOST_UBLAS_INLINE
    void
    swap_rows (const PM &pm, MV &mv, ublas::vector_tag)
    {
      typedef typename PM::size_type size_type;
      typedef typename MV::value_type value_type;

      size_type size = pm.size ();
      for (size_type 
          i = 0; 
          i < size; 
          ++ i)
      {
        value_type t (mv (i));
        mv (i) = mv (pm (i));
        mv (pm (i)) = t;
      }
    }

    /// Swap rows of a matrix
    template<class PM, class MV>
    BOOST_UBLAS_INLINE
    void
    swap_rows (const PM &pm, MV &mv, ublas::matrix_tag)
    {
      typedef typename PM::size_type size_type;
      typedef typename MV::value_type value_type;

      size_type size = pm.size ();
      for (size_type 
          i = 0; 
          i < size; 
          ++ i)
        for (size_type 
            j = 0; 
            j < mv.size2 (); 
            ++ j)
        {
          value_type t (mv (i, j));
          mv (i, j) = mv (pm (i), j);
          mv (pm (i), j) = t;
        }
    }

    /// Swap rows of a matrix or a vector: dispatcher
    template<class PM, class MV>
    BOOST_UBLAS_INLINE
    void
    swap_rows (const PM &pm, MV &mv)
    { swap_rows (pm, mv, BOOST_UBLAS_TYPENAME MV::type_category ()); }

    /// LU factorize: return LU in place of matrix
    template<class M, class PM>
    typename M::size_type lu_factorize (M &m, PM &pm)
    {
      typedef M matrix_type;
      typedef BOOST_UBLAS_TYPENAME M::size_type size_type;
      typedef BOOST_UBLAS_TYPENAME M::value_type value_type;

#ifdef BOOST_UBLAS_TYPE_CHECK
      matrix_type cm (m);
#endif
      int singular = 0;
      size_type size1 = m.size1 ();
      size_type size2 = m.size2 ();
      size_type size = std::min (size1, size2);
      for (size_type 
          i = 0; 
          i < size; 
          ++ i) 
      {
        matrix_column<M> mci (column (m, i));
        matrix_row<M> mri (row (m, i));
        size_type i_norm_inf = i + index_norm_inf (project (mci, range (i, size1)));
        BOOST_UBLAS_CHECK (i_norm_inf < m.size1 (), external_logic ());
        if (m (i_norm_inf, i) != value_type ())
        {
          pm (i) = i_norm_inf;
          if (i_norm_inf != i)
              row (m, i_norm_inf).swap (mri);
          project (mci, range (i + 1, size1)) *= value_type (1) / m (i, i);
        } else if (singular == 0)
          singular = i + 1;
        project (m, range (i + 1, size1), range (i + 1, size2)).minus_assign (
          outer_prod (project (mci, range (i + 1, size1)),
                      project (mri, range (i + 1, size2))));
      }
#ifdef BOOST_UBLAS_TYPE_CHECK
      swap_rows (pm, cm);
      BOOST_UBLAS_CHECK (
        singular != 0 ||
        equals (prod (triangular_adaptor<matrix_type, unit_lower> (m),
                      triangular_adaptor<matrix_type, upper> (m)), cm), internal_logic ());
#endif
      return singular;
    }

    /// LU solve: linear solve using LU and permutation
    template<class M, class PM, class MV>
    void
    lu_substitute (M &m, const PM &pm, MV &mv)
    {
      typedef M matrix_type;
      typedef MV matrix_vector_type;

      swap_rows (pm, mv);
#ifdef BOOST_UBLAS_TYPE_CHECK
      matrix_vector_type cmv1 (mv);
#endif
      inplace_solve (triangular_adaptor<matrix_type, unit_lower> (m),
                     mv, unit_lower_tag (),
                     BOOST_UBLAS_TYPENAME MV::type_category ());
#ifdef BOOST_UBLAS_TYPE_CHECK
      BOOST_UBLAS_CHECK (equals
        (prod(triangular_adaptor<matrix_type, unit_lower> (m), mv), cmv1), internal_logic ());
      matrix_vector_type cmv2 (mv);
#endif
      inplace_solve (triangular_adaptor<matrix_type, upper> (m),
                     mv, upper_tag (),
                     BOOST_UBLAS_TYPENAME MV::type_category ());
#ifdef BOOST_UBLAS_TYPE_CHECK
      BOOST_UBLAS_CHECK (equals
        (prod (triangular_adaptor<matrix_type, upper> (m), mv), cmv2), internal_logic ());
#endif
    }
  }

} }
#endif  // _GLUCAT_HAVE_UBLAS_LU_H
#endif  // _GLUCAT_UBLAS_IMP_H
