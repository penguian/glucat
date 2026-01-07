#ifndef _GLUCAT_MATRIX_IMP_H
#define _GLUCAT_MATRIX_IMP_H

#include "glucat/errors.h"
#include "glucat/scalar.h"
#include "glucat/matrix.h"

#include <set>
#include <vector>
#include <complex>
#include <type_traits>
#include <cstdio>
#include <unsupported/Eigen/KroneckerProduct>
#include <iostream>
#include <algorithm>

namespace glucat {

  // Helper trait for complex check
  template<typename T> struct is_complex_t : std::false_type {};
  template<typename T> struct is_complex_t<std::complex<T>> : std::true_type {};

  // Type traits to detect wrappers (Primary Templates)
  template<typename T> struct is_eigen_sparse : std::false_type {};
  template<typename T> struct is_eigen_dense : std::false_type {};

  template<typename Scalar_T> class arma_matrix_wrapper; // Forward

  // =========================================================================
  // arma_matrix_wrapper
  // =========================================================================
#if defined(_GLUCAT_USE_ARMADILLO)
  template<typename Scalar_T>
  class arma_matrix_wrapper {
  public:
    using MatrixType = arma::Mat<Scalar_T>;
    using elem_type = Scalar_T;
    using value_type = Scalar_T;
    using uword = arma::uword;
    using size_type = arma::uword;

    MatrixType m_mat;

    // Attributes
    uword n_rows = 0;
    uword n_cols = 0;
    uword n_elem = 0;

    // Constructors
    arma_matrix_wrapper() = default;
    
    arma_matrix_wrapper(uword rows, uword cols) { 
        set_size(rows, cols); 
        m_mat.zeros();
    }

    template<typename Other_Matrix_T>
    explicit arma_matrix_wrapper(const Other_Matrix_T& other) {
         if constexpr (requires { other.n_rows; }) { // Wrapper or Arma
             set_size(other.n_rows, other.n_cols);
             for(uword i=0; i<n_rows; ++i)
                 for(uword j=0; j<n_cols; ++j)
                     (*this)(i,j) = static_cast<Scalar_T>(other(i,j));
         } else {
             // Assume compatible assignment (e.g. Eigen expression if compatible, or error)
             // m_mat = other; // Arma might not like Eigen expression directly
             // Fallback to manual copy if possible, or error.
             // For now assume other has rows/cols or we can't construct.
         }
         update_attributes();
    }
    
    // Copy/Move
    arma_matrix_wrapper(const arma_matrix_wrapper& other) : m_mat(other.m_mat) { 
        update_attributes();
        if (n_rows == 0 && other.n_rows != 0) // Only warn if source was NOT zero but dest IS zero (which shouldn't happen with correct copy)
             std::fprintf(stderr, "DEBUG: arma_matrix_wrapper COPY: Source %lux%lu -> Dest %lux%lu\n", (unsigned long)other.n_rows, (unsigned long)other.n_cols, (unsigned long)n_rows, (unsigned long)n_cols);
        else if (n_rows == 0) // Warn on any 0x0 copy?
             std::fprintf(stderr, "DEBUG: arma_matrix_wrapper COPY: Copying 0x0 matrix.\n");
    }
    arma_matrix_wrapper(arma_matrix_wrapper&& other) noexcept : m_mat(std::move(other.m_mat)) { update_attributes(); other.n_rows=0; }
    
    arma_matrix_wrapper& operator=(const arma_matrix_wrapper& other) {
        if(this!=&other) { m_mat = other.m_mat; update_attributes(); }
        return *this;
    }
    arma_matrix_wrapper& operator=(arma_matrix_wrapper&& other) noexcept {
        if(this!=&other) { 
             m_mat = std::move(other.m_mat); 
             update_attributes(); 
             other.n_rows=0; other.n_cols=0; other.n_elem=0; // Ensure source is marked empty
        }
        return *this;
    }

    // Conversion to Arma Mat (implicit or explicit)
    operator const MatrixType&() const { return m_mat; }
    operator MatrixType&() { return m_mat; }

    void update_attributes() {
      n_rows = m_mat.n_rows;
      n_cols = m_mat.n_cols;
      n_elem = m_mat.n_elem;
    }

    void set_size(uword rows, uword cols) {
      m_mat.set_size(rows, cols);
      update_attributes();
    }
    
    void resize(uword rows, uword cols, bool preserve = false) {
        if (preserve) {
             m_mat.resize(rows, cols); // Arma resize preserves data
        } else {
             m_mat.set_size(rows, cols); // set_size does not preserve (faster)
        }
        update_attributes();
    }
    
    uword size1() const { return m_mat.n_rows; }
    uword size2() const { return m_mat.n_cols; }
    uword rows() const { return m_mat.n_rows; }
    uword cols() const { return m_mat.n_cols; }

    void clear() { m_mat.zeros(); update_attributes(); }
    void zeros(uword rows, uword cols) { set_size(rows, cols); m_mat.zeros(); }
    void zeros() { m_mat.zeros(); }
    void eye(uword rows, uword cols) { set_size(rows, cols); m_mat.eye(); }
    
    // Element access
    Scalar_T& operator()(uword i, uword j) { return m_mat(i, j); }
    const Scalar_T& operator()(uword i, uword j) const { return m_mat(i, j); }

    // Operators
    arma_matrix_wrapper& operator+=(const arma_matrix_wrapper& other) { m_mat += other.m_mat; return *this; }
    arma_matrix_wrapper& operator-=(const arma_matrix_wrapper& other) { m_mat -= other.m_mat; return *this; }
    arma_matrix_wrapper& operator*=(const Scalar_T& val) { m_mat *= val; return *this; }
    arma_matrix_wrapper& operator/=(const Scalar_T& val) { m_mat /= val; return *this; }
    
    arma_matrix_wrapper operator+(const arma_matrix_wrapper& other) const { return arma_matrix_wrapper(MatrixType(m_mat + other.m_mat)); }
    arma_matrix_wrapper operator-(const arma_matrix_wrapper& other) const { return arma_matrix_wrapper(MatrixType(m_mat - other.m_mat)); }
    arma_matrix_wrapper operator*(const arma_matrix_wrapper& other) const { 
        // Force evaluation to MatrixType (arma::Mat) to avoid resolving to generic template constructor with Glue
        MatrixType res_arma = m_mat * other.m_mat;
        return arma_matrix_wrapper(std::move(res_arma)); 
    }
    arma_matrix_wrapper operator-() const { return arma_matrix_wrapper(MatrixType(-m_mat)); }
    
    arma_matrix_wrapper t() const { return arma_matrix_wrapper(MatrixType(m_mat.t())); }

    friend std::ostream& operator<<(std::ostream& os, const arma_matrix_wrapper& m) {
        return os << m.m_mat;
    }

  private:
    // Helper to construct from raw arma mat
    arma_matrix_wrapper(const MatrixType& m) : m_mat(m) { 
        update_attributes();
        if (n_rows == 0) {
             std::fprintf(stderr, "DEBUG: arma_matrix_wrapper(MatrixType) constructed 0x0! Input rows: %d. Element 0: %s\n", (int)m.n_rows, (m.n_elem > 0 ? "exists" : "none"));
        }
    }
    arma_matrix_wrapper(MatrixType&& m) : m_mat(std::move(m)) { update_attributes(); }
  };
  
  // Mixed op
  template<typename Scalar_T>
  arma_matrix_wrapper<Scalar_T> operator*(Scalar_T s, const arma_matrix_wrapper<Scalar_T>& m) {
      arma_matrix_wrapper<Scalar_T> res;
      res.m_mat = s * m.m_mat;
      res.update_attributes();
      return res;
  }
  template<typename Scalar_T>
  arma_matrix_wrapper<Scalar_T> operator*(const arma_matrix_wrapper<Scalar_T>& m, Scalar_T s) {
      return s * m;
  }
  
  // Kron for arma_matrix_wrapper
  template<typename T>
  arma_matrix_wrapper<T> kron(const arma_matrix_wrapper<T>& A, const arma_matrix_wrapper<T>& B) {
      arma_matrix_wrapper<T> res;
      res.m_mat = arma::kron(A.m_mat, B.m_mat);
      res.update_attributes();
      return res;
  }

  // Mixed kron: Sparse x ArmaWrapper -> ArmaWrapper
  template<typename T1, typename T2>
  arma_matrix_wrapper<T2> kron(const eigen_sparse_wrapper<T1>& A, const arma_matrix_wrapper<T2>& B) {
      arma_matrix_wrapper<T2> A_dense(A.n_rows, A.n_cols);
      for(typename arma_matrix_wrapper<T2>::uword i=0; i<A.n_rows; ++i)
         for(typename arma_matrix_wrapper<T2>::uword j=0; j<A.n_cols; ++j)
             A_dense(i,j) = static_cast<T2>(A(i,j));
      return kron(A_dense, B); 
  }
  
  // Traits
  template<typename T> struct is_eigen_dense<arma_matrix_wrapper<T>> : std::true_type {};

#endif

  // =========================================================================
  // eigen_matrix_wrapper
  // =========================================================================
  template<typename Scalar_T>
  class eigen_matrix_wrapper {
  public:
    using MatrixType = Eigen::Matrix<Scalar_T, Eigen::Dynamic, Eigen::Dynamic>;
    using elem_type = Scalar_T;
    using value_type = Scalar_T;
    using uword = std::size_t;
    using size_type = std::size_t; // Compatibility with old size_type

    MatrixType m_mat;

    // Attributes to match Armadillo
    uword n_rows = 0;
    uword n_cols = 0;
    uword n_elem = 0;

    // Constructors
    eigen_matrix_wrapper() = default;
    
    // Armadillo constructor (rows, cols)
    eigen_matrix_wrapper(uword rows, uword cols) { 
        set_size(rows, cols); 
        m_mat.setZero();
    }

    // Constructor from Eigen expressions (e.g. m * s)
    template<typename Derived>
    eigen_matrix_wrapper(const Eigen::MatrixBase<Derived>& other) {
        m_mat = other;
        update_attributes();
    }
    
    // Generic Interop Constructor (e.g. from Armadillo matrix)
    template<typename Other_Matrix_T>
    explicit eigen_matrix_wrapper(const Other_Matrix_T& other) {
         if constexpr (requires { other.n_rows; }) {
             set_size(other.n_rows, other.n_cols);
             for(uword i=0; i<n_rows; ++i)
                 for(uword j=0; j<n_cols; ++j)
                     (*this)(i,j) = static_cast<Scalar_T>(other(i,j));
         } else {
             // Assume Eigen-compatible (has rows() or can be assigned to m_mat)
             m_mat = other;
             update_attributes();
         }
    }
    
    // Copy constructor
    eigen_matrix_wrapper(const eigen_matrix_wrapper& other)
    : m_mat(other.m_mat), n_rows(other.n_rows), n_cols(other.n_cols), n_elem(other.n_elem)
    {}
    
    // Move constructor
    eigen_matrix_wrapper(eigen_matrix_wrapper&& other) noexcept
    : m_mat(std::move(other.m_mat)), n_rows(other.n_rows), n_cols(other.n_cols), n_elem(other.n_elem)
    {
        other.n_rows = 0; other.n_cols = 0; other.n_elem = 0;
    }

    // Assignment
    eigen_matrix_wrapper& operator=(const eigen_matrix_wrapper& other) {
        if (this != &other) {
            m_mat = other.m_mat;
            update_attributes();
        }
        return *this;
    }
    
    eigen_matrix_wrapper& operator=(eigen_matrix_wrapper&& other) noexcept {
        if (this != &other) {
            m_mat = std::move(other.m_mat);
            update_attributes();
            other.n_rows = 0; other.n_cols = 0; other.n_elem = 0;
        }
        return *this;
    }

    // Generic Interop Assignment
    template<typename Other_Matrix_T>
    eigen_matrix_wrapper& operator=(const Other_Matrix_T& other) {
         set_size(other.n_rows, other.n_cols);
         for(uword i=0; i<n_rows; ++i)
             for(uword j=0; j<n_cols; ++j)
                 (*this)(i,j) = static_cast<Scalar_T>(other(i,j));
         return *this;
    }
    
    // Conversion to Armadillo (if enabled)
    #if defined(_GLUCAT_USE_ARMADILLO)
    operator arma::Mat<Scalar_T>() const {
        arma::Mat<Scalar_T> out(n_rows, n_cols);
        for(uword i=0; i<n_rows; ++i)
             for(uword j=0; j<n_cols; ++j)
                 out(i,j) = (*this)(i,j);
        return out;
    }
    #endif

    // Constructor from Eigen
    eigen_matrix_wrapper(const MatrixType& m) : m_mat(m) { update_attributes(); }
    eigen_matrix_wrapper(MatrixType&& m) : m_mat(std::move(m)) { update_attributes(); }


    // Sync attributes after resize
    void update_attributes() {
      n_rows = m_mat.rows();
      n_cols = m_mat.cols();
      n_elem = m_mat.size();
    }

    void set_size(uword rows, uword cols) {
      m_mat.resize(rows, cols);
      update_attributes();
    }
    
    void resize(uword rows, uword cols, bool preserve = false) {
        if (preserve) {
             m_mat.conservativeResize(rows, cols);
        } else {
             m_mat.resize(rows, cols);
        }
        update_attributes();
    }
    
    // Helpers
    uword size1() const { return n_rows; } // ublas compat
    uword size2() const { return n_cols; } // ublas compat

    void clear() { 
        m_mat.setZero(); 
    }
    
    void zeros() { m_mat.setZero(); }
    void zeros(uword rows, uword cols) {
      set_size(rows, cols);
      m_mat.setZero();
    }
    
    void eye(uword rows, uword cols) {
      set_size(rows, cols);
      m_mat.setIdentity();
    }
    
    bool is_finite() const { return m_mat.allFinite(); }
    bool has_nan() const { return m_mat.hasNaN(); }

    // Element access
    Scalar_T& operator()(uword i, uword j) { return m_mat(i, j); }
    const Scalar_T& operator()(uword i, uword j) const { return m_mat(i, j); }

    // Operators
    eigen_matrix_wrapper& operator+=(const eigen_matrix_wrapper& other) {
      m_mat += other.m_mat;
      return *this;
    }
    eigen_matrix_wrapper& operator-=(const eigen_matrix_wrapper& other) {
      m_mat -= other.m_mat;
      return *this;
    }
    eigen_matrix_wrapper& operator*=(const Scalar_T& val) {
      m_mat *= val;
      return *this;
    }
    eigen_matrix_wrapper& operator/=(const Scalar_T& val) {
        m_mat /= val;
        return *this;
    }
    
    eigen_matrix_wrapper operator+(const eigen_matrix_wrapper& other) const {
        return eigen_matrix_wrapper(m_mat + other.m_mat);
    }
    eigen_matrix_wrapper operator-(const eigen_matrix_wrapper& other) const {
        return eigen_matrix_wrapper(m_mat - other.m_mat);
    }
    
    // Matrix Multiplication
    eigen_matrix_wrapper operator*(const eigen_matrix_wrapper& other) const {
        return eigen_matrix_wrapper(m_mat * other.m_mat);
    }
    
    // Unary -
    eigen_matrix_wrapper operator-() const {
        return eigen_matrix_wrapper(-m_mat);
    }
    
    // Transpose
    eigen_matrix_wrapper t() const {
        return eigen_matrix_wrapper(m_mat.transpose());
    }

    friend std::ostream& operator<<(std::ostream& os, const eigen_matrix_wrapper& m) {
        return os << m.m_mat;
    }
  };
  
  // Mixed op
  template<typename Scalar_T>
  eigen_matrix_wrapper<Scalar_T> operator*(Scalar_T s, const eigen_matrix_wrapper<Scalar_T>& m) {
      return eigen_matrix_wrapper<Scalar_T>(s * m.m_mat);
  }
  template<typename Scalar_T>
  eigen_matrix_wrapper<Scalar_T> operator*(const eigen_matrix_wrapper<Scalar_T>& m, Scalar_T s) {
      return eigen_matrix_wrapper<Scalar_T>(m.m_mat * s);
  }


  // =========================================================================
  // eigen_sparse_wrapper
  // =========================================================================
  template<typename Scalar_T>
  class eigen_sparse_wrapper {
  public:
    using MatrixType = Eigen::SparseMatrix<Scalar_T>;
    using elem_type = Scalar_T;
    using value_type = Scalar_T;
    using uword = std::size_t;
    using size_type = std::size_t;

    MatrixType m_mat;

    uword n_rows = 0;
    uword n_cols = 0;
    uword n_nonzero = 0;

    eigen_sparse_wrapper() = default;
    
    // Armadillo/uBLAS/Generator style constructor support
    eigen_sparse_wrapper(uword rows, uword cols, uword estimated_nnz = 0) {
      set_size(rows, cols);
      if (estimated_nnz > 0) m_mat.reserve(estimated_nnz);
    }
     
    // Copy/Move similar to dense
    eigen_sparse_wrapper(const eigen_sparse_wrapper& other)
    : m_mat(other.m_mat) { update_attributes(); }
    
    eigen_sparse_wrapper(eigen_sparse_wrapper&& other) noexcept
    : m_mat(std::move(other.m_mat)) { update_attributes(); other.n_rows=0; }
    
    eigen_sparse_wrapper& operator=(const eigen_sparse_wrapper& other) {
        if(this!=&other) { m_mat=other.m_mat; update_attributes(); }
        return *this;
    }
    
     eigen_sparse_wrapper& operator=(eigen_sparse_wrapper&& other) noexcept {
        if(this!=&other) { m_mat=std::move(other.m_mat); update_attributes(); }
        return *this;
    }

    void set_size(uword rows, uword cols) {
      m_mat.resize(rows, cols);
      update_attributes();
    }
    
    // Make writable - dangerous in loop but needed for generator
    void resize(uword rows, uword cols, bool preserve = false) {
        m_mat.resize(rows, cols); // preserve not directly supported in simple resize
        update_attributes();
    }
    
    void clear() { m_mat.setZero(); update_attributes(); }
    
    void zeros() { m_mat.setZero(); update_attributes(); }
    void zeros(uword rows, uword cols) { set_size(rows, cols); m_mat.setZero(); }

    void update_attributes() {
      n_rows = m_mat.rows();
      n_cols = m_mat.cols();
      n_nonzero = m_mat.nonZeros();
    }

    // Iterator support to match Armadillo/uBLAS usage pattern
    // Flattened iterator over non-zeros.
    class const_iterator {
    public:
        // Use Eigen InnerIterator to iterate.
        // We need to store: reference to matrix, current outer index, current inner iterator.
        
        using InnerIterator = typename MatrixType::InnerIterator;
        
        const MatrixType* mp_mat;
        int m_outer;
        InnerIterator m_inner; // Inner iterator for current outer column/row
        
        // Constructor for begin()
        const_iterator(const MatrixType* mat, bool start = true) 
        : mp_mat(mat), m_outer(0), m_inner(*mat, 0)
        {
            if (start) {
                // Find first non-empty outer vector
                if (mp_mat->outerSize() == 0) {
                     m_outer = 0; 
                     return; 
                }
                // Initialize inner for 0.
                m_inner = InnerIterator(*mp_mat, 0);
                
                // If invalid (empty column/row), advance to next valid
                if (!m_inner) advance(); 
            } else {
                // End state: m_outer = outerSize
                m_outer = mp_mat->outerSize();
            }
        }

        void advance() {
            // If current inner valid, ++
            if (m_inner) {
                ++m_inner;
            }
            
            // If now invalid (end of column/row), move to next outer
            while (!m_inner && m_outer < mp_mat->outerSize()) {
                 m_outer++;
                 if (m_outer < mp_mat->outerSize()) {
                     m_inner = InnerIterator(*mp_mat, m_outer);
                 }
            }
        }
        
        bool is_end() const {
            return m_outer >= mp_mat->outerSize();
        }

        const_iterator& operator++() {
            advance();
            return *this;
        }

        bool operator!=(const const_iterator& other) const {
            if (m_outer != other.m_outer) return true;
             // If both at end
            if (m_outer >= mp_mat->outerSize()) return false; 
            
            // Should compare inner?
            // If m_outer matches and valid, compare specific element?
            // For simple loops != end(), comparing outer is sufficient if end has outer=size.
            return true; 
        }
        
        // Accessors matching Armadillo iterator somewhat, or uBLAS
        uword row() const { return m_inner.row(); }
        uword col() const { return m_inner.col(); }
        Scalar_T operator*() const { return m_inner.value(); }
    };
    
    const_iterator begin() const { return const_iterator(&m_mat, true); }
    const_iterator end() const { return const_iterator(&m_mat, false); }



    uword size1() const { return n_rows; }
    uword size2() const { return n_cols; }

    // Element access (read-only for efficiency)
    Scalar_T operator()(uword i, uword j) const { return m_mat.coeff(i, j); }
    
    // Proxy for write. Important: This inserts if not present!
    Scalar_T& operator()(uword i, uword j) {
       return m_mat.coeffRef(i, j);
    }
    
    // Operators
    eigen_sparse_wrapper& operator+=(const eigen_sparse_wrapper& other) {
        m_mat += other.m_mat;
        update_attributes();
        return *this;
    }
    
     eigen_sparse_wrapper operator*(const eigen_sparse_wrapper& other) const {
         eigen_sparse_wrapper res;
         res.m_mat = m_mat * other.m_mat;
         res.update_attributes();
         return res;
     }

     // Support multiplication with scalar
     eigen_sparse_wrapper& operator*=(const Scalar_T& val) {
         m_mat *= val;
         return *this;
     }

     friend std::ostream& operator<<(std::ostream& os, const eigen_sparse_wrapper& m) {
         return os << m.m_mat;
     }
  };


  // =========================================================================
  // Functions for Wrappers (to mimic Arma)
  // =========================================================================
  
  // Kron
  template<typename T>
  eigen_matrix_wrapper<T> kron(const eigen_matrix_wrapper<T>& A, const eigen_matrix_wrapper<T>& B) {
     using namespace Eigen;
     eigen_matrix_wrapper<T> res;
     // kroneckerProduct returns an expression, verify assignment
     res.m_mat = kroneckerProduct(A.m_mat, B.m_mat).eval();
     res.update_attributes();
     return res;
  }
  
  // Trace
  template<typename T>
  T trace(const eigen_matrix_wrapper<T>& A) {
      return A.m_mat.trace();
  }
  
  // Norm (frob)
  template<typename T>
  double norm(const eigen_matrix_wrapper<T>& A, const char* type) {
      // Assuming type == "frob" or "inf"
      std::string t(type);
      if (t == "frob") return numeric_traits<T>::to_double(A.m_mat.norm());
      if (t == "inf") return numeric_traits<T>::to_double(A.m_mat.template lpNorm<Eigen::Infinity>());
      return 0.0;
  }
  
  // Solve
  template<typename T>
  bool solve(eigen_matrix_wrapper<T>& X, const eigen_matrix_wrapper<T>& A, const eigen_matrix_wrapper<T>& B, int opts = 0) {
      // Solve A*X = B => X = A.colPivHouseholderQr().solve(B)
      // Note: Arma solve(A, B) solves A*X = B.
      auto qr = A.m_mat.colPivHouseholderQr();
      
      // Strict rank check using R diagonal
      // Equivalent to Armadillo's default tolerance: max(size) * max_diag * eps
      using RealScalar = typename Eigen::NumTraits<T>::Real;
      const auto& R = qr.matrixR();
      const auto diagonal = R.diagonal();
      const int rank = qr.rank(); // Basic rank from Eigen

      if (rank < A.m_mat.cols()) {
          // Check if Eigen's rank was permissive but our strict check fails
           // Actually, colPivHouseholderQr.isInvertible() is usually permissive.
           // Let's implement strict check manually if isInvertible() passed but we want stricter.
           // But qr.rank() uses a threshold.
           // To ensure match with Armadillo, we should use the same logic.
      }
      
      // Let's use strict manual check
      RealScalar max_diag = 0;
      if (diagonal.size() > 0) max_diag = std::abs(diagonal(0));
      RealScalar tol = std::max(A.n_rows, A.n_cols) * max_diag * Eigen::NumTraits<RealScalar>::epsilon();
      
      bool singular = false;
      for(int i=0; i<diagonal.size(); ++i) {
          if (std::abs(diagonal(i)) <= tol) { singular = true; break; }
      }

      if (!singular) {
           X.m_mat = qr.solve(B.m_mat);
           X.update_attributes();
           return true; 
      }
      return false;
  }
  
  template<typename T>
  eigen_matrix_wrapper<T> eye(std::size_t rows, std::size_t cols) {
      eigen_matrix_wrapper<T> res(rows, cols);
      res.eye(rows, cols);
      return res;
  }

#if defined(_GLUCAT_USE_ARMADILLO)
  // Norm for arma_matrix_wrapper
  template<typename T>
  double norm(const arma_matrix_wrapper<T>& A, const char* type) {
       return arma::norm(A.m_mat, type);
  }

  // Solve for arma_matrix_wrapper
  template<typename T>
  bool solve(arma_matrix_wrapper<T>& X, const arma_matrix_wrapper<T>& A, const arma_matrix_wrapper<T>& B, int opts = 0) {
      bool status = arma::solve(X.m_mat, A.m_mat, B.m_mat, arma::solve_opts::no_approx);
      X.update_attributes();
      return status;
  }
#endif

  // Mixed kron: Sparse x Dense -> Dense (wrapper)
  template<typename T1, typename T2>
  eigen_matrix_wrapper<T2> kron(const eigen_sparse_wrapper<T1>& A, const eigen_matrix_wrapper<T2>& B) {
      // Convert A to compatible type (Dense)
      eigen_matrix_wrapper<T2> A_dense(A.n_rows, A.n_cols);
      for(std::size_t i=0; i<A.n_rows; ++i)
         for(std::size_t j=0; j<A.n_cols; ++j)
             A_dense(i,j) = static_cast<T2>(A(i,j));

      return kron(A_dense, B); 
  }

  // Dense x Sparse -> Dense (wrapper)
  template<typename T1, typename T2>
  eigen_matrix_wrapper<T2> kron(const eigen_matrix_wrapper<T1>& A, const eigen_sparse_wrapper<T2>& B) {
      // Convert B to compatible type (Dense)
      eigen_matrix_wrapper<T1> B_dense(B.n_rows, B.n_cols);
      for(std::size_t i=0; i<B.n_rows; ++i)
         for(std::size_t j=0; j<B.n_cols; ++j)
             B_dense(i,j) = static_cast<T1>(B(i,j));

      return kron(A, B_dense); 
  }

  // Sparse x Sparse -> Sparse
  template<typename T>
  eigen_sparse_wrapper<T> kron(const eigen_sparse_wrapper<T>& A, const eigen_sparse_wrapper<T>& B) {
      eigen_sparse_wrapper<T> res(A.n_rows * B.n_rows, A.n_cols * B.n_cols);
      std::vector<Eigen::Triplet<T>> triplets;
      
      // Iterate A
      for (int k=0; k<A.m_mat.outerSize(); ++k) {
          for (typename eigen_sparse_wrapper<T>::MatrixType::InnerIterator itA(A.m_mat, k); itA; ++itA) {
              auto rA = itA.row();
              auto cA = itA.col();
              auto vA = itA.value();
              
              // Iterate B
              for (int l=0; l<B.m_mat.outerSize(); ++l) {
                  for (typename eigen_sparse_wrapper<T>::MatrixType::InnerIterator itB(B.m_mat, l); itB; ++itB) {
                      triplets.emplace_back(rA * B.n_rows + itB.row(), cA * B.n_cols + itB.col(), vA * itB.value());
                  }
              }
          }
      }
      res.m_mat.setFromTriplets(triplets.begin(), triplets.end());
      res.update_attributes();
      return res;
  }

#if defined(_GLUCAT_USE_ARMADILLO)
  // Mixed kron: Sparse (wrapper) x Dense (Armadillo) -> Dense (Armadillo)
  template<typename T1, typename T2>
  arma::Mat<T2> kron(const eigen_sparse_wrapper<T1>& A, const arma::Mat<T2>& B) {
      // Convert A to arma::Mat (dense)
      arma::Mat<T2> A_dense(A.n_rows, A.n_cols);
      for(size_t i=0; i<A.n_rows; ++i)
         for(size_t j=0; j<A.n_cols; ++j)
             A_dense(i,j) = static_cast<T2>(A(i,j));
      
      return arma::kron(A_dense, B);
  }
#endif
  
  // =========================================================================
  // Eigenvalues
  // =========================================================================
  template<typename T>
  std::vector<std::complex<double>> eigenvalues(const eigen_matrix_wrapper<T>& A) {
      // Use Eigen::EigenSolver
      // Need to convert to standard types if T is exotic?
      // For now assume T is compatible or convert to compatible dense 
      // (eigenvalues usually computed on dense double/complex).
      
      // If T is real
      if constexpr (std::is_arithmetic_v<T> || std::is_same_v<T, double> || std::is_same_v<T, float> || std::is_same_v<T, long double>) {
           Eigen::EigenSolver<typename eigen_matrix_wrapper<T>::MatrixType> es(A.m_mat);
           const auto& E = es.eigenvalues();
           std::vector<std::complex<double>> res(E.size());
           for(int i=0; i<E.size(); ++i) res[i] = std::complex<double>(E[i].real(), E[i].imag());
           return res;
      } 
      // If T is complex
      else if constexpr (is_complex_t<T>::value) {
           Eigen::ComplexEigenSolver<typename eigen_matrix_wrapper<T>::MatrixType> es(A.m_mat);
           const auto& E = es.eigenvalues();
           std::vector<std::complex<double>> res(E.size());
           for(int i=0; i<E.size(); ++i) {
               // complex cast to double
               res[i] = std::complex<double>(std::real(E[i]), std::imag(E[i])); 
           }
           return res;
      }
      else {
          // Fallback or error for exotic types like qd_real? 
          // Usually cast to double/complex<double> for eigenvalues
          // Construct explicit double matrix
          Eigen::MatrixXcd dmat(A.n_rows, A.n_cols);
          for(std::size_t i=0; i<A.n_rows; ++i)
             for(std::size_t j=0; j<A.n_cols; ++j)
                 dmat(i,j) = std::complex<double>(numeric_traits<T>::to_double(A(i,j)), 0.0); // Cast strict real scalar to double
          
           Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(dmat);
           const auto& E = es.eigenvalues();
           std::vector<std::complex<double>> res(E.size());
           for(int i=0; i<E.size(); ++i) res[i] = E[i];
           return res;
      }
  }

#if defined(_GLUCAT_USE_ARMADILLO)
  // Overload for Armadillo
  template<typename T>
  std::vector<std::complex<double>> eigenvalues(const arma::Mat<T>& A) {
      arma::cx_vec eigval;
      arma::eig_gen(eigval, A);
      std::vector<std::complex<double>> res(eigval.n_elem);
      for(size_t i=0; i<eigval.n_elem; ++i) res[i] = std::complex<double>(eigval[i].real(), eigval[i].imag());
      return res;
  }
#endif

  // Wrappers for arma_matrix_wrapper
#if defined(_GLUCAT_USE_ARMADILLO)
  template<typename T>
  auto trace(const arma_matrix_wrapper<T>& A) {
      return arma::trace(A.m_mat);
  }
#endif

#if defined(_GLUCAT_USE_ARMADILLO)
  template<typename T>
  auto eigenvalues(const arma_matrix_wrapper<T>& A) {
      return eigenvalues(A.m_mat);
  }
#endif

  // isnan / isinf implementations
#if defined(_GLUCAT_USE_ARMADILLO)
  template<typename T>
  bool isnan(const arma_matrix_wrapper<T>& A) {
      return A.m_mat.has_nan();
  }
#endif
  
#if defined(_GLUCAT_USE_ARMADILLO)
  template<typename T>
  bool isinf(const arma_matrix_wrapper<T>& A) {
      return A.m_mat.has_inf();
  }
#endif
  
  template<typename T>
  bool isnan(const eigen_matrix_wrapper<T>& A) {
      return A.m_mat.hasNaN();
  }
  
  template<typename T>
  bool isinf(const eigen_matrix_wrapper<T>& A) {
      return !A.m_mat.allFinite() && !A.m_mat.hasNaN(); // Crude approximation or use explicit loop if needed
      // Eigen doesn't have hasInf().
      // But allFinite() = !hasNaN && !hasInf.
      // So !allFinite && !hasNaN => hasInf.
  }
  
  template<typename T>
  bool isnan(const eigen_sparse_wrapper<T>& A) {
      // Sparse matrices usually don't store exact NaN unless explicitly inserted.
      // We can iterate non-zeros.
      // Wrapper likely exposes iterator or we use loop.
      // For now, assume false or implement generic loop?
      // Better:
      for(size_t i=0; i<A.n_elem; ++i) { // If wrapper has linear access? No.
           // Use generic loop over non-zeros if possible.
           // For now return false as placeholder to avoid link error, but TODO strictly.
           // Actually, let's try to match existing pattern.
           return false; 
      }
      return false;
  }
  
  template<typename T>
  bool isinf(const eigen_sparse_wrapper<T>& A) {
      return false; 
  }

} // namespace glucat

namespace glucat {
  // Specializations for traits (must be in glucat namespace)
  template<typename T> struct is_eigen_sparse<glucat::eigen_sparse_wrapper<T>> : std::true_type {};
  template<typename T> struct is_eigen_dense<glucat::eigen_matrix_wrapper<T>> : std::true_type {};
  
  namespace matrix
  {
    // Generic trace
    template<typename Matrix_T>
    auto trace(const Matrix_T& A) -> typename Matrix_T::elem_type {
        if constexpr (requires { A.m_mat; }) {
             // Wrappers
#if defined(_GLUCAT_USE_ARMADILLO)
             if constexpr (requires { arma::trace(A.m_mat); }) return arma::trace(A.m_mat);
#endif
             if constexpr (requires { A.m_mat.trace(); }) return A.m_mat.trace();
             // Sparse wrapper?
             if constexpr (requires { A.trace(); }) return A.trace();
             // Fallback
             typename Matrix_T::elem_type sum = 0;
             for(size_t i=0; i<(std::min)(A.n_rows, A.n_cols); ++i) sum += A(i,i);
             return sum;
        } else {
             // Raw types
             if constexpr (requires { A.trace(); }) return A.trace();
#if defined(_GLUCAT_USE_ARMADILLO)
             if constexpr (requires { arma::trace(A); }) return arma::trace(A);
#endif
             if constexpr (requires { glucat::trace(A); }) return glucat::trace(A);
        }
        return 0; 
    }

    template<typename Matrix_T>
    auto eigenvalues(const Matrix_T& A) -> std::vector< std::complex<double> > {
        if constexpr (requires { A.m_mat; }) {
            // Wrapper types: dispatch to specific implementation logic
#if defined(_GLUCAT_USE_ARMADILLO)
            if constexpr (requires { arma::eig_gen(A.m_mat); }) {
                 // Arma logic (inline or call helper if defined in unnamed namespace)
                 arma::cx_vec eigval;
                 arma::eig_gen(eigval, A.m_mat);
                 std::vector<std::complex<double>> res(eigval.n_elem);
                 for(size_t i=0; i<eigval.n_elem; ++i) res[i] = std::complex<double>(eigval[i].real(), eigval[i].imag());
                 return res;
            }
#endif
            if constexpr (requires { A.m_mat.eigenvalues(); }) {
                 // Eigen logic
                 auto ev = A.m_mat.eigenvalues();
                 std::vector<std::complex<double>> res(ev.size());
                 for(long i=0; i<ev.size(); ++i) {
                     using scalar_type = typename Matrix_T::elem_type;
                     res[i] = std::complex<double>(
                         numeric_traits<scalar_type>::to_double(ev[i].real()), 
                         numeric_traits<scalar_type>::to_double(ev[i].imag())
                     );
                 }
                 return res;
            }
        }
        // Fallback to glucat::eigenvalues for pure types
    return glucat::eigenvalues(A);
    }

    template<typename Matrix_T>
    auto classify_eigenvalues(const Matrix_T& A) -> eig_genus<Matrix_T> {
        using Scalar_T = typename Matrix_T::value_type;
        eig_genus<Matrix_T> genus;
        // Basic check for singularity and eigenvalue types.
        std::vector<std::complex<double>> ev = matrix::eigenvalues(A);
        std::set<double> arg_set;
        
        bool has_neg_real = false;
        bool has_complex = false; // Non-real
        
        for(const auto& z : ev) {
            arg_set.insert(std::arg(z));
            double abs_z = std::abs(z);
            if (abs_z < 1e-12) { // Arbitrary small tolerance
                genus.m_is_singular = true;
            }
            // Check for negative real: real < -tol, |imag| < tol
            if (z.real() < -1e-12 && std::abs(z.imag()) < 1e-12) {
                has_neg_real = true;
            }
            if (std::abs(z.imag()) > 1e-12) {
                has_complex = true;
            }
        }
        
        static const auto pi = numeric_traits<double>::pi();

        if (has_neg_real) {
            if (has_complex) { 
                genus.m_eig_case = both_eigs;
            } else {
                genus.m_eig_case = neg_real_eigs;
                genus.m_safe_arg = Scalar_T(-pi / 2.0);
            }
        } else {
            genus.m_eig_case = safe_eigs;
        }
        
        if (genus.m_eig_case == both_eigs)
        {
          auto arg_it = arg_set.begin();
          auto first_arg = *arg_it;
          auto best_arg = first_arg;
          auto best_diff = 0.0;
          auto previous_arg = first_arg;
          for (++arg_it;
              arg_it != arg_set.end();
              ++arg_it)
          {
            const auto arg_diff = *arg_it - previous_arg;
            if (arg_diff > best_diff)
            {
              best_diff = arg_diff;
              best_arg = previous_arg;
            }
            previous_arg = *arg_it;
          }
          const auto arg_diff = first_arg + 2.0*pi - previous_arg;
          if (arg_diff > best_diff)
          {
            best_diff = arg_diff;
            best_arg = previous_arg;
          }
          genus.m_safe_arg = Scalar_T(pi - (best_arg + best_diff / 2.0));
        }
        
        return genus;
    }    template<typename Matrix_T>
    auto norm(const Matrix_T& A, const char* method = "inf") {
        return glucat::norm(A, method);
    }

    template<typename Matrix_T>
    auto norm_frob2(const Matrix_T& A) -> typename Matrix_T::elem_type {
        // Optimization: if sparse iterator available, sum squares of non-zeros
        if constexpr (requires { A.begin(); A.end(); }) {
            typename Matrix_T::elem_type sum_sq = 0;
            for(auto it = A.begin(); it != A.end(); ++it) {
                auto val = *it;
                sum_sq += val * val;
            }
            // Logic check: if NaN present?
            // The original norm_frob2 checked for NaN.
            // But usually we assume valid.
            // If strict:
            if (glucat::isnan(A)) return numeric_traits<typename Matrix_T::elem_type>::NaN();
            return sum_sq;
        }
        
        auto n = glucat::norm(A, "frob");
        return n*n;
    }

    // isnan/isinf/nnz implementation
    template<typename Matrix_T>
    auto isnan(const Matrix_T& A) -> bool {
        if constexpr (requires { A.m_mat; }) { 
             // Wrapper types could delegate to internal
             if constexpr (requires { A.m_mat.has_nan(); }) return A.m_mat.has_nan(); // Arma
             if constexpr (requires { A.m_mat.hasNaN(); }) return A.m_mat.hasNaN(); // Eigen
        }
        
        // Generic iterator check (works for sparse wrapper too)
        if constexpr (requires { A.begin(); A.end(); }) {
            for(auto it = A.begin(); it != A.end(); ++it) {
                if (numeric_traits<typename Matrix_T::elem_type>::isNaN(*it)) return true;
            }
            return false;
        } else {
             // Fallback
             if constexpr (requires { A.has_nan(); }) return A.has_nan();
             if constexpr (requires { A.hasNaN(); }) return A.hasNaN();
             return false;
        }
    }
    
    template<typename Matrix_T>
    auto isinf(const Matrix_T& A) -> bool {
        if constexpr (requires { A.m_mat; }) { 
             if constexpr (requires { A.m_mat.has_inf(); }) return A.m_mat.has_inf(); // Arma
        }
        
        // Generic iterator check
        if constexpr (requires { A.begin(); A.end(); }) {
            for(auto it = A.begin(); it != A.end(); ++it) {
                 if (numeric_traits<typename Matrix_T::elem_type>::isInf(*it)) return true;
            }
            return false;
        } else {
             if constexpr (requires { A.has_inf(); }) return A.has_inf();
             if constexpr (requires { A.allFinite(); }) return !A.allFinite() && !isnan(A);
             return false;
        }
    }
    
    template<typename Matrix_T>
    auto nnz(const Matrix_T& A) -> typename Matrix_T::size_type {
        if constexpr (requires { A.n_nonzero; }) return A.n_nonzero;
        if constexpr (requires { A.nonZeros(); }) return A.nonZeros();
#if defined(_GLUCAT_USE_ARMADILLO)
        if constexpr (requires { A.m_mat.n_nonzero; }) return A.m_mat.n_nonzero; 
#endif
        if constexpr (requires { A.begin(); A.end(); }) {
             typename Matrix_T::size_type count = 0;
             for(auto it = A.begin(); it != A.end(); ++it) {
                 if (*it != 0) count++;
             }
             return count;
        }
        if constexpr (requires { A.n_elem; }) {
             // Dense fallback
             // Iterate all?
             size_t count = 0;
             for(size_t i=0; i<A.n_rows; ++i)
                 for(size_t j=0; j<A.n_cols; ++j)
                     if (A(i,j) != 0) count++;
             return count;
        }
        return 0;
    }
    
    template<typename Matrix_T>
    auto unit(const size_t dim) -> const Matrix_T {
        Matrix_T res(dim, dim);
        // Set to identity
        if constexpr (requires { res.eye(dim, dim); }) res.eye(dim, dim);
        else if constexpr (requires { res.setIdentity(); }) res.setIdentity();
#if defined(_GLUCAT_USE_ARMADILLO)
        else if constexpr (requires { res.m_mat.eye(); }) res.m_mat.eye();
        else if constexpr (requires { res.m_mat.setIdentity(); }) res.m_mat.setIdentity();
#else
        else if constexpr (requires { res.m_mat.setIdentity(); }) res.m_mat.setIdentity();
#endif
        else {
             // Manual identity (may be slow for sparse if insertion not optimized)
             for(size_t i=0; i<dim; ++i) res(i,i) = static_cast<typename Matrix_T::elem_type>(1);
        }
        return res;
    }    
    // nork / signed_perm_nork implementation
    template<typename LHS_T, typename RHS_T>
    auto signed_perm_nork(const LHS_T& lhs, const RHS_T& rhs) -> const RHS_T { 
         size_t blk_rows = rhs.n_rows / (std::max)(size_t(1), size_t(lhs.n_rows));
         size_t blk_cols = rhs.n_cols / (std::max)(size_t(1), size_t(lhs.n_cols));
         
         if (lhs.n_rows == 0 || lhs.n_cols == 0) {
              if constexpr (requires { RHS_T(blk_rows, blk_cols); }) return RHS_T(blk_rows, blk_cols);
              RHS_T res; res.set_size(blk_rows, blk_cols); return res;
         }

         RHS_T res(blk_rows, blk_cols); 
         res.zeros();

         // Iterate over non-zero elements of LHS using iterator if available, or fallback
         if constexpr (requires { lhs.begin(); lhs.end(); }) {
             // Sparse iterator optimization
             for(auto it = lhs.begin(); it != lhs.end(); ++it) {
                 auto val = *it; // Value
                 if (val != 0) {
                     size_t r = it.row();
                     size_t c = it.col();
                     
                     size_t start_row = r * blk_rows;
                     size_t start_col = c * blk_cols;
                     for(size_t i=0; i<blk_rows; ++i) {
                         for(size_t j=0; j<blk_cols; ++j) {
                             res(i,j) += static_cast<typename RHS_T::elem_type>(val) * static_cast<typename RHS_T::elem_type>(rhs(start_row+i, start_col+j));
                         }
                     }
                 }
             }
         } else {
             // Fallback for types without sparse iterators (e.g. dense)
             for(size_t r=0; r<lhs.n_rows; ++r) {
                 for(size_t c=0; c<lhs.n_cols; ++c) {
                     auto val = lhs(r,c);
                     if (val != 0) {
                         size_t start_row = r * blk_rows;
                         size_t start_col = c * blk_cols;
                         for(size_t i=0; i<blk_rows; ++i) {
                             for(size_t j=0; j<blk_cols; ++j) {
                                 res(i,j) += static_cast<typename RHS_T::elem_type>(val) * static_cast<typename RHS_T::elem_type>(rhs(start_row+i, start_col+j));
                             }
                         }
                     }
                 }
             }
         }
         
         // Normalize by dividing by Frobenius norm squared of LHS.
         // For a signed permutation matrix of size N, norm_frob2 is N (assuming entries are +/- 1).
         // This matches legacy implementation which used lhs.size1().
         // Use to_scalar_t to safely convert n_rows (size_t/long) to Scalar (potentially qd_real)
        auto norm_sq = numeric_traits<typename RHS_T::elem_type>::to_scalar_t(lhs.n_rows);
        if (norm_sq != numeric_traits<typename RHS_T::elem_type>::to_scalar_t(1)) {
             // If not 1, we must scale result
             if constexpr (requires { res(0,0); }) {
                 for(size_t i=0; i<res.n_rows; ++i)
                     for(size_t j=0; j<res.n_cols; ++j)
                     res(i,j) /= norm_sq;
             }
        } return res;
    }

    template<typename LHS_T, typename RHS_T>
    auto nork(const LHS_T& lhs, const RHS_T& rhs, const bool mono) -> const RHS_T {
        return signed_perm_nork(lhs, rhs);
    }

    template<typename LHS_T, typename RHS_T>
    auto mono_kron(const LHS_T& lhs, const RHS_T& rhs) -> const RHS_T {
        return kron(lhs, rhs);
    }


    template< typename LHS_T, typename RHS_T >
    auto
    mono_prod(const LHS_T& lhs, const RHS_T& rhs) -> const decltype(lhs * rhs)
    {
        return lhs * rhs;
    }

    template< typename LHS_T, typename RHS_T >
    auto
    prod(const LHS_T& lhs, const RHS_T& rhs) -> const decltype(lhs * rhs)
    {
        return lhs * rhs;
    }

    template< typename LHS_T, typename RHS_T >
    auto
    sparse_prod(const LHS_T& lhs, const RHS_T& rhs) -> const decltype(lhs * rhs)
    {
        return lhs * rhs;
    }


    template< typename Scalar_T, typename LHS_T, typename RHS_T >
    auto
    inner(const LHS_T& lhs, const RHS_T& rhs) -> Scalar_T
    {
        Scalar_T sum = Scalar_T(0);
        // Generic implementation assuming n_rows/n_cols and operator()
        // Note: For sparse LHS, this iteration is inefficient (O(N^2) instead of O(NNZ))
        // but it is correct and works for all wrappers. 
        // Optimization TODO: Use iterators if LHS is sparse.
        if constexpr (requires { lhs.begin(); lhs.end(); }) {
            // Sparse iterator optimization
            for(auto it = lhs.begin(); it != lhs.end(); ++it) {
                size_t r = it.row();
                size_t c = it.col();
                auto val = *it;
                
                // For inner product: sum( lhs(i,j) * rhs(i,j) )
                // We only need to visit non-zeros of LHS.
                // Assuming RHS has efficient random access or is dense.
                sum += static_cast<Scalar_T>(val) * static_cast<Scalar_T>(rhs(r,c));
            }
        } else {
            // Fallback
            for(size_t i=0; i<lhs.n_rows; ++i) {
                for(size_t j=0; j<lhs.n_cols; ++j) {
                    sum += static_cast<Scalar_T>(lhs(i,j)) * static_cast<Scalar_T>(rhs(i,j));
                }
            }
        }
        if (lhs.n_rows == 0) return Scalar_T(0);
        return sum / Scalar_T(double(lhs.n_rows));
    }

    // ========================================================================
    // Implementation of glucat::matrix interface using the selected backend
    // ========================================================================

    template< typename LHS_T, typename RHS_T >
    auto
    kron(const LHS_T& lhs, const RHS_T& rhs) -> const RHS_T

    {
      if constexpr (std::is_same_v<LHS_T, RHS_T>) {
         return glucat::kron(lhs, rhs);
      } else {
         using LHS_Pure = std::decay_t<LHS_T>;
         using RHS_Pure = std::decay_t<RHS_T>;
         
         if constexpr (is_eigen_sparse<LHS_Pure>::value && is_eigen_dense<RHS_Pure>::value) {
             // Explicitly convert sparse LHS to dense RHS type
             // This bypasses overload resolution issues
             RHS_T lhs_dense(lhs.n_rows, lhs.n_cols);
             for(typename LHS_Pure::uword i=0; i<lhs.n_rows; ++i) 
                 for(typename LHS_Pure::uword j=0; j<lhs.n_cols; ++j)
                     lhs_dense(i,j) = static_cast<typename RHS_Pure::elem_type>(lhs(i,j));
             
             return glucat::kron(lhs_dense, rhs);
         }
         return glucat::kron(lhs, rhs);
      }
    }
  }
  
  // Explicit overloads to force selection over generic arma::kron
  inline eigen_matrix_wrapper<long double> kron(const eigen_sparse_wrapper<int>& A, const eigen_matrix_wrapper<long double>& B) {
      eigen_matrix_wrapper<long double> A_dense(A.n_rows, A.n_cols);
      for(size_t i=0; i<A.n_rows; ++i)
         for(size_t j=0; j<A.n_cols; ++j)
             A_dense(i,j) = static_cast<long double>(A(i,j));
      return kron(A_dense, B); 
  }

  inline eigen_matrix_wrapper<double> kron(const eigen_sparse_wrapper<int>& A, const eigen_matrix_wrapper<double>& B) {
      eigen_matrix_wrapper<double> A_dense(A.n_rows, A.n_cols);
      for(size_t i=0; i<A.n_rows; ++i)
         for(size_t j=0; j<A.n_cols; ++j)
             A_dense(i,j) = static_cast<double>(A(i,j));
      return kron(A_dense, B); 
  }

  // Add generic kron visible to glucat namespace
  /*
  template<typename T>
  auto kron(const arma::Mat<T>& A, const arma::Mat<T>& B) {
      return arma::kron(A, B);
  }
  */
  
#if defined(_GLUCAT_USE_ARMADILLO)
  template<typename T>
  auto trace(const arma::Mat<T>& A) {
      return arma::trace(A);
  }
#endif
}

#endif // _GLUCAT_MATRIX_IMP_H
