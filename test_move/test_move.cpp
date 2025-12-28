#include <iostream>
#include <type_traits>
#include "test/driver.h"
// #include "glucat/index_set_imp.h"
// #include "glucat/framed_multi_imp.h"
// #include "glucat/matrix_multi_imp.h"

template <typename T>
void check_move_semantics(const std::string& name) {
    std::cout << "Checking " << name << ":\n";
    std::cout << "  is_move_constructible: " << std::is_move_constructible<T>::value << "\n";
    std::cout << "  is_move_assignable:    " << std::is_move_assignable<T>::value << "\n";

    if (std::is_move_constructible<T>::value && std::is_move_assignable<T>::value) {
        std::cout << "  [PASS] Move semantics enabled.\n";
    } else {
        std::cout << "  [FAIL] Move semantics NOT enabled.\n";
    }
    std::cout << "----------------------------------------\n";
}

int main() {
    using namespace glucat;

    check_move_semantics<index_set<0, 32>>("index_set<0, 32>");
    check_move_semantics<framed_multi<double, 0, 32>>("framed_multi<double, 0, 32>");
    check_move_semantics<matrix_multi<double, 0, 32>>("matrix_multi<double, 0, 32>");

    return 0;
}
