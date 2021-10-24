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
  auto
  operator* (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> const framed_multi<Scalar_T,LO,HI>;

  /// Outer product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator^ (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> const framed_multi<Scalar_T,LO,HI>;

  /// Inner product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator& (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> const framed_multi<Scalar_T,LO,HI>;

  /// Left contraction
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator% (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> const framed_multi<Scalar_T,LO,HI>;

  /// Hestenes scalar product
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  star(const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> Scalar_T;

  /// Geometric quotient
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator/ (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> const framed_multi<Scalar_T,LO,HI>;

  /// Transformation via twisted adjoint action
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator| (const framed_multi<Scalar_T,LO,HI>& lhs, const framed_multi<Scalar_T,LO,HI>& rhs) -> const framed_multi<Scalar_T,LO,HI>;

  /// Read multivector from input
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator>> (std::istream& s, framed_multi<Scalar_T,LO,HI>& val) -> std::istream&;

  /// Write multivector to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator<< (std::ostream& os, const framed_multi<Scalar_T,LO,HI>& val) -> std::ostream&;

  /// Write term to output
  template< typename Scalar_T, const index_t LO, const index_t HI >
  std::ostream&
  operator<< (std::ostream& os, const std::pair< const index_set<LO,HI>, Scalar_T >& term);

  /// Exponential of multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  exp(const framed_multi<Scalar_T,LO,HI>& val) -> const framed_multi<Scalar_T,LO,HI>;

  template< const index_t LO, const index_t HI>
  class index_set_hash
  {
  public:
    using index_set_t = index_set<LO, HI>;
    inline auto operator()(index_set_t val) const -> size_t { return val.hash_fn(); }
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
    using multivector_t = framed_multi;
    using framed_multi_t = multivector_t;
    using scalar_t = Scalar_T;
    using index_set_t = index_set<LO, HI>;
    using term_t = std::pair<const index_set_t, Scalar_T>;
    using vector_t = std::vector<Scalar_T>;
    using error_t = error<multivector_t>;
    using matrix_multi_t = matrix_multi<Scalar_T, LO, HI>;
    template< typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI >
    friend class matrix_multi;
    template< typename Other_Scalar_T, const index_t Other_LO, const index_t Other_HI >
    friend class framed_multi;

  private:
    class                                              var_term; // forward
    using var_term_t = class var_term;
    using matrix_t = typename matrix_multi_t::matrix_t;
    using sorted_map_t = std::map<index_set_t, Scalar_T, std::less<const index_set_t>>;
#if defined(_GLUCAT_USE_STD_UNORDERED_MAP)
    using map_t = std::unordered_map<index_set_t, Scalar_T, index_set_hash<LO, HI>>;
#else
    typedef sorted_map_t                               map_t;
#endif

    class hash_size_t
    {
    public:
      hash_size_t(size_t hash_size)
      : n(hash_size)
      { };
      auto operator()() const -> size_t
      { return n; }
    private:
      size_t n;
    };

    using framed_pair_t = std::pair<const multivector_t, const multivector_t>;
    using size_type = typename map_t::size_type;
    using iterator = typename map_t::iterator;
    using const_iterator = typename map_t::const_iterator;

  public:
    /// Class name used in messages
    static auto classname() -> const std::string;
    /// Destructor
    ~framed_multi() override = default;
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
    auto fast_matrix_multi(const index_set_t frm) const -> const matrix_multi<Other_Scalar_T,LO,HI>;
    /// Use inverse generalized FFT to construct a framed_multi_t
    auto fast_framed_multi() const -> const framed_multi_t;

    _GLUCAT_CLIFFORD_ALGEBRA_OPERATIONS

    /// Number of terms
    auto nbr_terms() const -> unsigned long;

    /// Random multivector within a frame
    static auto random(const index_set_t frm, Scalar_T fill = Scalar_T(1)) -> const framed_multi_t;

    // Friend declarations

    friend auto
      operator* <>(const framed_multi_t& lhs, const framed_multi_t& rhs) -> const framed_multi_t;
    friend auto
      operator^ <>(const framed_multi_t& lhs, const framed_multi_t& rhs) -> const framed_multi_t;
    friend auto
      operator& <>(const framed_multi_t& lhs, const framed_multi_t& rhs) -> const framed_multi_t;
    friend auto
      operator% <>(const framed_multi_t& lhs, const framed_multi_t& rhs) -> const framed_multi_t;
    friend auto
      star      <>(const framed_multi_t& lhs, const framed_multi_t& rhs) -> Scalar_T;
    friend auto
      operator/ <>(const framed_multi_t& lhs, const framed_multi_t& rhs) -> const framed_multi_t;
    friend auto
      operator| <>(const framed_multi_t& lhs, const framed_multi_t& rhs) -> const framed_multi_t;

    friend auto
      operator>> <>(std::istream& s, multivector_t& val) -> std::istream&;
    friend auto
      operator<< <>(std::ostream& os, const multivector_t& val) -> std::ostream&;
    friend auto
      operator<< <>(std::ostream& os, const term_t& term) -> std::ostream&;

    friend auto
      exp <>(const framed_multi_t& val) -> const framed_multi_t;

    /// Add a term, if non-zero
    auto      operator+= (const term_t& term) -> multivector_t&;

  private:
    /// Subalgebra isomorphism: fold each term within the given frame
    auto       fold(const index_set_t frm) const -> multivector_t;
    /// Subalgebra isomorphism: unfold each term within the given frame
    auto       unfold(const index_set_t frm) const -> multivector_t;
    /// Subalgebra isomorphism: R_{p,q} to R_{p-4,q+4}
    auto      centre_pm4_qp4(index_t& p, index_t& q) -> multivector_t&;
    /// Subalgebra isomorphism: R_{p,q} to R_{p+4,q-4}
    auto      centre_pp4_qm4(index_t& p, index_t& q) -> multivector_t&;
    /// Subalgebra isomorphism: R_{p,q} to R_{q+1,p-1}
    auto      centre_qp1_pm1(index_t& p, index_t& q) -> multivector_t&;
    /// Divide multivector into part divisible by index_set and remainder
    auto      divide(const index_set_t ist) const -> const framed_pair_t;
    /// Generalized FFT from framed_multi_t to matrix_t
    auto      fast(const index_t level, const bool odd) const -> const matrix_t;

    /// Variable term
    class var_term :
    public std::pair<index_set<LO,HI>, Scalar_T>
    {
    public:
      using var_pair_t = std::pair<index_set<LO, HI>, Scalar_T>;

      /// Class name used in messages
      static auto classname() -> const std::string
      { return "var_term"; };
      /// Destructor
      ~var_term() = default;
      /// Default constructor
      var_term()
      : var_pair_t(index_set_t(), Scalar_T(1))
      { };
      /// Construct a variable term from an index set and a scalar coordinate
      var_term(const index_set_t ist, const Scalar_T& crd = Scalar_T(1))
      : var_pair_t(ist, crd)
      { };
      /// Product of variable term and term
      auto operator*= (const term_t& rhs) -> var_term_t&
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
  auto
  crd_of_mult(const std::pair<const index_set<LO,HI>, Scalar_T>& lhs,
              const std::pair<const index_set<LO,HI>, Scalar_T>& rhs) -> Scalar_T;

  /// Product of terms
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  operator*
   (const std::pair<const index_set<LO,HI>, Scalar_T>& lhs,
    const std::pair<const index_set<LO,HI>, Scalar_T>& rhs) -> const std::pair<const index_set<LO,HI>, Scalar_T>;

  /// Square root of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  sqrt(const framed_multi<Scalar_T,LO,HI>& val, const framed_multi<Scalar_T,LO,HI>& i, bool prechecked) -> const framed_multi<Scalar_T,LO,HI>;

  /// Exponential of multivector
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  exp(const framed_multi<Scalar_T,LO,HI>& val) -> const framed_multi<Scalar_T,LO,HI>;

  /// Natural logarithm of multivector with specified complexifier
  template< typename Scalar_T, const index_t LO, const index_t HI >
  auto
  log(const framed_multi<Scalar_T,LO,HI>& val, const framed_multi<Scalar_T,LO,HI>& i, bool prechecked) -> const framed_multi<Scalar_T,LO,HI>;
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
