#ifndef _GLUCAT_FRAMED_MULTI_H
#define _GLUCAT_FRAMED_MULTI_H
/***************************************************************************
    GluCat : Generic library of universal Clifford algebra templates
    framed_multi.h : Declare a class for the framed representation of a multivector
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
#include "glucat/clifford_algebra.h"

#if defined(_GLUCAT_USE_BOOST_POOL_ALLOC)
// Use the Boost pool allocator
#include <boost/pool/poolfwd.hpp>
#endif

#include <string>
#include <utility>
#include <map>
#include <vector>

// Use the appropriate type of map

#if defined(_GLUCAT_USE_STD_UNORDERED_MAP)
# include <unordered_map>
#endif

#if defined(_GLUCAT_USE_STD_UNORDERED_MAP)
# define _GLUCAT_MAP_IS_HASH
#else
# define _GLUCAT_MAP_IS_ORDERED
#endif

namespace glucat
{
  // Forward declarations for friends

  template< typename Scalar_T, const index_t LO, const index_t HI >
  class framed_multi; // forward

  template< typename Scalar_T, const index_t LO, const index_t HI >
  class matrix_multi; // forward

  /// Geometric product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  operator* (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs);

  /// Outer product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  operator^ (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs);

  /// Inner product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  operator& (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs);

  /// Left contraction
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  operator% (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs);

  /// Hestenes scalar product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  Scalar_T
  star(const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs);

  /// Geometric quotient
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  operator/ (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs);

  /// Transformation via twisted adjoint action
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  operator| (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs);

  /// Read multivector from input
  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::istream&
  operator>> (std::istream& s, framed_multi<Scalar_T,LO,HI>& val);

  /// Write multivector to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::ostream&
  operator<< (std::ostream& os, const framed_multi<Scalar_T,LO,HI>& val);

  /// Write term to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::ostream&
  operator<< (std::ostream& os, const std::pair< const index_set<LO,HI>, Scalar_T >& term);

  /// Exponential of multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  exp(const framed_multi<Scalar_T,LO,HI>& val);

  template< const index_t LO, const index_t HI>
  class index_set_hash
  {
  public:
    typedef index_set<LO,HI> index_set_t;
    inline size_t operator()(index_set_t val) const { return val.hash_fn(); }
  };

  /// A framed_multi<Scalar_T,LO,HI> is a framed approximation to a multivector
  template< typename Scalar_T = double,  const index_t LO = DEFAULT_LO, const index_t HI = DEFAULT_HI >
  class framed_multi :
  public clifford_algebra< Scalar_T, index_set<LO,HI>, framed_multi<Scalar_T,LO,HI> >,
#if defined(_GLUCAT_USE_STD_UNORDERED_MAP)
  private std::unordered_map< index_set<LO,HI>, Scalar_T, index_set_hash<LO,HI> >
#else
  private std::map< index_set<LO,HI>, Scalar_T,
                    std::less< const index_set<LO,HI> >
#if defined(_GLUCAT_USE_BOOST_POOL_ALLOC)
                  , boost::fast_pool_allocator< std::pair<const index_set<LO,HI>, Scalar_T> >
#endif
                  >
#endif
  {
  public:
    typedef framed_multi                               multivector_t;
    typedef multivector_t                              framed_multi_t;
    typedef Scalar_T                                   scalar_t;
    typedef index_set<LO,HI>                           index_set_t;
    typedef std::pair<const index_set_t, Scalar_T>     term_t;
    typedef std::vector<Scalar_T>                      vector_t;
    typedef error<multivector_t>                       error_t;
    typedef      matrix_multi<Scalar_T,LO,HI>          matrix_multi_t;
    template< typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI >
    friend class matrix_multi;
    template< typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI >
    friend class framed_multi;

  private:
    class                                              var_term; // forward
    typedef class var_term                             var_term_t;
    typedef typename matrix_multi_t::matrix_t          matrix_t;
    typedef std::map< index_set_t, Scalar_T,
                      std::less<const index_set_t>
#if defined(_GLUCAT_USE_BOOST_POOL_ALLOC)
                    , boost::fast_pool_allocator<term_t>
#endif
                    >
                                                       sorted_map_t;
#if defined(_GLUCAT_USE_STD_UNORDERED_MAP)
    typedef std::unordered_map< index_set_t, Scalar_T, index_set_hash<LO,HI> >
                                                       map_t;
#else
    typedef sorted_map_t                               map_t;
#endif

    class hash_size_t
    {
    public:
      hash_size_t(size_t hash_size)
      : n(hash_size)
      { };
      size_t operator()() const
      { return n; }
    private:
      size_t n;
    };

    typedef std::pair< const multivector_t, const multivector_t >
                                                       framed_pair_t;
    typedef typename map_t::size_type                  size_type;
    typedef typename map_t::iterator                   iterator;
    typedef typename map_t::const_iterator             const_iterator;

  public:
    /// Class name used in messages
    static const std::string classname();
    /// Destructor
    ~framed_multi() {};
    /// Default constructor
    framed_multi();

  private:
    /// Private constructor using hash_size
    framed_multi(const hash_size_t& hash_size);
  public:
    /// Construct a multivector from a multivector with a different scalar type
    template< typename Other_Scalar_T >
    framed_multi(const framed_multi<Other_Scalar_T,LO,HI>& val);
    /// Construct a multivector, within a given frame, from a given multivector
    template< typename Other_Scalar_T >
    framed_multi(const framed_multi<Other_Scalar_T,LO,HI>& val,
                 const index_set_t frm, const bool prechecked = false);
    /// Construct a multivector, within a given frame, from a given multivector
    framed_multi(const framed_multi_t& val,
                 const index_set_t frm, const bool prechecked = false);
    /// Construct a multivector from an index set and a scalar coordinate
    framed_multi(const index_set_t ist, const Scalar_T& crd = Scalar_T(1));
    /// Construct a multivector, within a given frame, from an index set and a scalar coordinate
    framed_multi(const index_set_t ist, const Scalar_T& crd,
                 const index_set_t frm, const bool prechecked = false);
    /// Construct a multivector from a scalar (within a frame, if given)
    framed_multi(const Scalar_T& scr, const index_set_t frm = index_set_t());
    /// Construct a multivector from an int (within a frame, if given)
    framed_multi(const int scr, const index_set_t frm = index_set_t());
    /// Construct a multivector, within a given frame, from a given vector
    framed_multi(const vector_t& vec,
                 const index_set_t frm, const bool prechecked = false);
    /// Construct a multivector from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
    framed_multi(const std::string& str);
    /// Construct a multivector, within a given frame, from a string: eg: "3+2{1,2}-6.1e-2{2,3}"
    framed_multi(const std::string& str,
                 const index_set_t frm, const bool prechecked = false);
    /// Construct a multivector from a char*: eg: "3+2{1,2}-6.1e-2{2,3}"
    framed_multi(const char* str)
    { *this = framed_multi(std::string(str)); };
    /// Construct a multivector, within a given frame, from a char*: eg: "3+2{1,2}-6.1e-2{2,3}"
    framed_multi(const char* str,
                 const index_set_t frm, const bool prechecked = false)
    { *this = framed_multi(std::string(str), frm, prechecked); };
    /// Construct a multivector from a matrix_multi_t
    template< typename Other_Scalar_T >
    framed_multi(const matrix_multi<Other_Scalar_T,LO,HI>& val);
    /// Use generalized FFT to construct a matrix_multi_t
    template< typename Other_Scalar_T >
    const matrix_multi<Other_Scalar_T,LO,HI> fast_matrix_multi(const index_set_t frm) const;
    /// Use inverse generalized FFT to construct a framed_multi_t
    const framed_multi_t fast_framed_multi() const;

    _GLUCAT_CLIFFORD_ALGEBRA_OPERATIONS

    /// Number of terms
    unsigned long nbr_terms() const;

    /// Random multivector within a frame
    static const framed_multi_t random(const index_set_t frm, Scalar_T fill = Scalar_T(1));

    // Friend declarations

    friend const framed_multi_t
      operator* <>(const framed_multi_t& lhs, const framed_multi_t& rhs);
    friend const framed_multi_t
      operator^ <>(const framed_multi_t& lhs, const framed_multi_t& rhs);
    friend const framed_multi_t
      operator& <>(const framed_multi_t& lhs, const framed_multi_t& rhs);
    friend const framed_multi_t
      operator% <>(const framed_multi_t& lhs, const framed_multi_t& rhs);
    friend Scalar_T
      star      <>(const framed_multi_t& lhs, const framed_multi_t& rhs);
    friend const framed_multi_t
      operator/ <>(const framed_multi_t& lhs, const framed_multi_t& rhs);
    friend const framed_multi_t
      operator| <>(const framed_multi_t& lhs, const framed_multi_t& rhs);

    friend std::istream&
      operator>> <>(std::istream& s, multivector_t& val);
    friend std::ostream&
      operator<< <>(std::ostream& os, const multivector_t& val);
    friend std::ostream&
      operator<< <>(std::ostream& os, const term_t& term);

    friend const framed_multi_t
      exp <>(const framed_multi_t& val);

    /// Add a term, if non-zero
    multivector_t&      operator+= (const term_t& term);

  private:
    /// Subalgebra isomorphism: fold each term within the given frame
    multivector_t       fold(const index_set_t frm) const;
    /// Subalgebra isomorphism: unfold each term within the given frame
    multivector_t       unfold(const index_set_t frm) const;
    /// Subalgebra isomorphism: R_{p,q} to R_{p-4,q+4}
    multivector_t&      centre_pm4_qp4(index_t& p, index_t& q);
    /// Subalgebra isomorphism: R_{p,q} to R_{p+4,q-4}
    multivector_t&      centre_pp4_qm4(index_t& p, index_t& q);
    /// Subalgebra isomorphism: R_{p,q} to R_{q+1,p-1}
    multivector_t&      centre_qp1_pm1(index_t& p, index_t& q);
    /// Divide multivector into part divisible by index_set and remainder
    const framed_pair_t divide(const index_set_t ist) const;
    /// Generalized FFT from framed_multi_t to matrix_t
    const matrix_t      fast(const index_t level, const bool odd) const;

    /// Variable term
    class var_term :
    public std::pair<index_set<LO,HI>, Scalar_T>
    {
    public:
      typedef std::pair<index_set<LO,HI>, Scalar_T>      var_pair_t;

      /// Class name used in messages
      static const std::string classname()
      { return "var_term"; };
      /// Destructor
      ~var_term() {};
      /// Default constructor
      var_term()
      : var_pair_t(index_set_t(), Scalar_T(1))
      { };
      /// Construct a variable term from an index set and a scalar coordinate
      var_term(const index_set_t ist, const Scalar_T& crd = Scalar_T(1))
      : var_pair_t(ist, crd)
      { };
      /// Product of variable term and term
      var_term_t& operator*= (const term_t& rhs)
      {
        this->second *= rhs.second * this->first.sign_of_mult(rhs.first);
        this->first  ^= rhs.first;
        return *this;
      }
    };
  };

  // Non-members

  /// Coordinate of product of terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  inline
  static
  Scalar_T
  crd_of_mult(const std::pair<const index_set<LO,HI>, Scalar_T>& lhs,
              const std::pair<const index_set<LO,HI>, Scalar_T>& rhs);

  /// Product of terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const std::pair<const index_set<LO,HI>, Scalar_T>
  operator*
   (const std::pair<const index_set<LO,HI>, Scalar_T>& lhs,
    const std::pair<const index_set<LO,HI>, Scalar_T>& rhs);

  /// Square root of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  sqrt(const framed_multi<Scalar_T,LO,HI>& val, const framed_multi<Scalar_T,LO,HI>& i, bool prechecked);

  /// Exponential of multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  exp(const framed_multi<Scalar_T,LO,HI>& val);

  /// Natural logarithm of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  const framed_multi<Scalar_T,LO,HI>
  log(const framed_multi<Scalar_T,LO,HI>& val, const framed_multi<Scalar_T,LO,HI>& i, bool prechecked);
}

namespace std
{
  /// Numeric limits for framed_multi inherit limits for the corresponding scalar type
  template <typename Scalar_T, const glucat::index_t LO, const glucat::index_t HI>
  struct numeric_limits< glucat::framed_multi<Scalar_T,LO,HI> > :
  public numeric_limits<Scalar_T>
  { };
}
#endif  // _GLUCAT_FRAMED_MULTI_H
