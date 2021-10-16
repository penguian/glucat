#ifndef _GLUCAT_MATRIX_MULTI_H
#define _GLUCAT_MATRIX_MULTI_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    matrix_multi.h : Declare a class for the matrix representation of a multivector
                             -------------------
    begin                : Sun 2001-12-09
    copyright            : (C) 2001-2016 by Paul C. Leopardi
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
#include "glucat/errors.h"
#include "glucat/index_set.h"
#include "glucat/clifford_algebra.h"
#include "glucat/framed_multi.h"

#include <boost/numeric/ublas/fwd.hpp>

#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace glucat
{
  namespace ublas = boost::numeric::ublas;

  // Forward declarations for friends

  template< typename Scalar_T, const index_t LO, const index_t HI >
  class framed_multi;  // forward

  template< typename Scalar_T, const index_t LO, const index_t HI >
  class matrix_multi;  // forward

  /// Geometric product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  operator* (const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs);

  /// Outer product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  operator^ (const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs);

  /// Inner product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  operator& (const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs);

  /// Left contraction
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  operator% (const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs);

  /// Hestenes scalar product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  Scalar_T
  star(const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs);

  /// Geometric quotient
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  operator/ (const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs);

  /// Transformation via twisted adjoint action
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  operator| (const matrix_multi<Scalar_T,LO,HI>& lhs, const matrix_multi<Scalar_T,LO,HI>& rhs);

  /// Read multivector from input
  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::istream&
  operator>> (std::istream& s, matrix_multi<Scalar_T,LO,HI>& val);

  /// Write multivector to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::ostream&
  operator<< (std::ostream& os, const matrix_multi<Scalar_T,LO,HI>& val);

  /// Find a common frame for operands of a binary operator 
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const index_set<LO,HI>
  reframe (const matrix_multi<Scalar_T,LO,HI>& lhs,    const matrix_multi<Scalar_T,LO,HI>& rhs,
                 matrix_multi<Scalar_T,LO,HI>& lhs_reframed, matrix_multi<Scalar_T,LO,HI>& rhs_reframed);

  /// Square root of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  sqrt(const matrix_multi<Scalar_T,LO,HI>& val, const matrix_multi<Scalar_T,LO,HI>& i, bool prechecked);

  /// Square root of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  matrix_sqrt(const matrix_multi<Scalar_T,LO,HI>& val, const matrix_multi<Scalar_T,LO,HI>& i);

  /// Natural logarithm of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  log(const matrix_multi<Scalar_T,LO,HI>& val, const matrix_multi<Scalar_T,LO,HI>& i, bool prechecked);

  /// Natural logarithm of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  matrix_log(const matrix_multi<Scalar_T,LO,HI>& val, const matrix_multi<Scalar_T,LO,HI>& i);

  /// A matrix_multi<Scalar_T,LO,HI> is a matrix approximation to a multivector
  template< typename Scalar_T = double, const index_t LO = DEFAULT_LO, const index_t HI = DEFAULT_HI >
  class matrix_multi :
  public clifford_algebra< Scalar_T, index_set<LO,HI>, matrix_multi<Scalar_T,LO,HI> >
  {
  public:
    using multivector_t = matrix_multi;
    using matrix_multi_t = multivector_t;
    using scalar_t = Scalar_T;
    using index_set_t = index_set<LO, HI>;
    using term_t = std::pair<const index_set_t, Scalar_T>;
    using vector_t = std::vector<Scalar_T>;
    using error_t = error<multivector_t>;
    using framed_multi_t = framed_multi<Scalar_T, LO, HI>;
    template< typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI >
    friend class framed_multi;
    template< typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI >
    friend class matrix_multi;

  private:
    using orientation_t = ublas::row_major;
    using basis_matrix_t = ublas::compressed_matrix<int, orientation_t>;
#if defined(_GLUCAT_USE_DENSE_MATRICES)
    using matrix_t = ublas::matrix<Scalar_T, orientation_t>;
#else
    typedef ublas::compressed_matrix< Scalar_T, orientation_t >
                                                       matrix_t;
#endif
    using matrix_index_t = typename matrix_t::size_type;

  public:
    /// Class name used in messages
    static const std::string classname();
    /// Destructor
    ~matrix_multi() override = default;
    /// Default constructor
    matrix_multi();
    /// Construct a multivector from a multivector with a different scalar type
    template< typename Other_Scalar_T >
    matrix_multi(const matrix_multi<Other_Scalar_T,LO,HI>& val);
    /// Construct a multivector, within a given frame, from a given multivector
    template< typename Other_Scalar_T >
    matrix_multi(const matrix_multi<Other_Scalar_T,LO,HI>& val,
                 const index_set_t frm, const bool prechecked = false);
    /// Construct a multivector, within a given frame, from a given multivector
    matrix_multi(const multivector_t& val,
                 const index_set_t frm, const bool prechecked = false);
    /// Construct a multivector from an index set and a scalar coordinate
    matrix_multi(const index_set_t ist, const Scalar_T& crd = Scalar_T(1));
    /// Construct a multivector, within a given frame, from an index set and a scalar coordinate
    matrix_multi(const index_set_t ist, const Scalar_T& crd,
                 const index_set_t frm, const bool prechecked = false);
    /// Construct a multivector from a scalar (within a frame, if given)
    matrix_multi(const Scalar_T& scr, const index_set_t frm = index_set_t());
    /// Construct a multivector from an int (within a frame, if given)
    matrix_multi(const int scr, const index_set_t frm = index_set_t());
    /// Construct a multivector, within a given frame, from a given vector
    matrix_multi(const vector_t& vec,
                 const index_set_t frm, const bool prechecked = false);
    /// Construct a multivector from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
    matrix_multi(const std::string& str);
    /// Construct a multivector, within a given frame, from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
    matrix_multi(const std::string& str,
                 const index_set_t frm, const bool prechecked = false);
    /// Construct a multivector from a char*: eg: "3+2{1,2}-6.1e-2{2,3}"
    matrix_multi(const char* str)
    { *this = matrix_multi(std::string(str)); };
    /// Construct a multivector, within a given frame, from a char*: eg: "3+2{1,2}-6.1e-2{2,3}"
    matrix_multi(const char* str,
                 const index_set_t frm, const bool prechecked = false)
    { *this = matrix_multi(std::string(str), frm, prechecked); };
    /// Construct a multivector from a framed_multi_t
    template< typename Other_Scalar_T >
    matrix_multi(const framed_multi<Other_Scalar_T,LO,HI>& val);
    /// Construct a multivector, within a given frame, from a framed_multi_t
    template< typename Other_Scalar_T >
    matrix_multi(const framed_multi<Other_Scalar_T,LO,HI>& val,
                 const index_set_t frm, const bool prechecked = false);
    /// Use generalized FFT to construct a matrix_multi_t
    const matrix_multi_t fast_matrix_multi(const index_set_t frm) const;
    /// Use inverse generalized FFT to construct a framed_multi_t
    template< typename Other_Scalar_T >
    const framed_multi<Other_Scalar_T,LO,HI> fast_framed_multi() const;

  private:
    /// Construct a multivector within a given frame from a given matrix
    template< typename Matrix_T >
    matrix_multi(const Matrix_T& mtx, const index_set_t frm);
    /// Construct a multivector within a given frame from a given matrix
    matrix_multi(const matrix_t& mtx, const index_set_t frm);
    /// Create a basis element matrix within the current frame
    const basis_matrix_t basis_element(const index_set<LO,HI>& ist) const;

  public:
    _GLUCAT_CLIFFORD_ALGEBRA_OPERATIONS

    /// Assignment operator
    multivector_t&     operator= (const multivector_t& rhs);

    /// Random multivector within a frame
    static const matrix_multi_t random(const index_set_t frm, Scalar_T fill = Scalar_T(1));

    // Friend declarations

    friend const matrix_multi_t
      operator* <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);
    friend const matrix_multi_t
      operator^ <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);
    friend const matrix_multi_t
      operator& <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);
    friend const matrix_multi_t
      operator% <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);
    friend Scalar_T
      star      <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);
    friend const matrix_multi_t
      operator/ <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);
    friend const matrix_multi_t
      operator| <>(const matrix_multi_t& lhs, const matrix_multi_t& rhs);

    friend std::istream&
      operator>> <>(std::istream& s, multivector_t& val);
    friend std::ostream&
      operator<< <>(std::ostream& os, const multivector_t& val);
    friend std::ostream&
      operator<< <>(std::ostream& os, const term_t& term);

    template< typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI >
    friend const index_set<Other_LO,Other_HI>
    reframe (const matrix_multi<Other_Scalar_T,Other_LO,Other_HI>& lhs,    const matrix_multi<Other_Scalar_T,Other_LO,Other_HI>& rhs,
                   matrix_multi<Other_Scalar_T,Other_LO,Other_HI>& lhs_reframed, matrix_multi<Other_Scalar_T,Other_LO,Other_HI>& rhs_reframed);
    template< typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI >
    friend const matrix_multi<Other_Scalar_T,Other_LO,Other_HI>
      matrix_sqrt(const matrix_multi<Other_Scalar_T,Other_LO,Other_HI>& val, const matrix_multi<Other_Scalar_T,Other_LO,Other_HI>& i);
    template< typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI >
    friend const matrix_multi<Other_Scalar_T,Other_LO,Other_HI>
      matrix_log(const matrix_multi<Other_Scalar_T,Other_LO,Other_HI>& val, const matrix_multi<Other_Scalar_T,Other_LO,Other_HI>& i);

    /// Add a term, if non-zero
    multivector_t&     operator+= (const term_t& rhs);

   private:
    // Data members

    /// Index set representing the frame for the subalgebra which contains the multivector
    index_set_t        m_frame;
    /// Matrix value representing the multivector within the folded frame
    matrix_t           m_matrix;
  };

  // Non-members

  /// Exponential of multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const matrix_multi<Scalar_T,LO,HI>
  exp(const matrix_multi<Scalar_T,LO,HI>& val);

}

namespace std
{
  /// Numeric limits for matrix_multi inherit limits for the corresponding scalar type
  template <typename Scalar_T, const glucat::index_t LO, const glucat::index_t HI>
  struct numeric_limits< glucat::matrix_multi<Scalar_T,LO,HI> > :
  public numeric_limits<Scalar_T>
  { };
}
#endif  // _GLUCAT_MATRIX_MULTI_H
