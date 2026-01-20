
#define _GLUCAT_MATRIX_DEBUG
#include "glucat/glucat.h"

#include "glucat/clifford_algebra_imp.h"
#include "glucat/generation_imp.h"
#include "glucat/errors_imp.h"


#include "glucat/index_set_imp.h"
#include "glucat/framed_multi_imp.h"


#include "glucat/matrix_imp.h"
#include "glucat/matrix_multi_imp.h"


#include <iostream>

using namespace glucat;

int main() {
    using Scalar = float;
    // framed_multi<float,-5,6>
    // p=6, q=5.
    using Multi = framed_multi<Scalar, -5, 6>;
    using MatrixMulti = matrix_multi<Scalar, -5, 6>;

    try {
        index_set< -5, 6 > frame; 
        // Populate frame with some indices? 
        // Default frame is empty? No, framed_multi allows dynamic frames within LO, HI.
        // But matrix_multi constructor uses provided frame.
        
        // Identity 1.22a failed for r=1, s=1.
        // Let's create two vectors.
        Multi a(index_set< -5, 6 >(1), 1.0f); // e_1
        Multi b(index_set< -5, 6 >(2), 1.0f); // e_2
        
        std::cout << "a: " << a << std::endl;
        std::cout << "b: " << b << std::endl;
        
        Multi prod = a * b; // Should be e_1 * e_2 = e_12 (bivector)
        std::cout << "prod (framed): " << prod << std::endl;
        
        // Convert to matrix
        MatrixMulti m(prod); 
        std::cout << "m: " << m << std::endl;
        // If size is large, prod * prod might use matrix. But here we just convert.
        
        // Convert back
        Multi prod_back(m);
        std::cout << "prod_back: " << prod_back << std::endl;
        
        if (prod != prod_back) {
            std::cout << "ROUND TRIP FAILED" << std::endl;
        } else {
            std::cout << "ROUND TRIP PASSED" << std::endl;
        }
        
        // Check grade 2 part
        Multi prod_grade2 = prod(2);
        Multi prod_back_grade2 = prod_back(2);
        
        std::cout << "prod(2): " << prod_grade2 << std::endl;
        std::cout << "prod_back(2): " << prod_back_grade2 << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
    return 0;
}
