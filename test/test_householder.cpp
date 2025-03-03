#include <print>
#include<vector>

import math;

void test_householder(){

    jf::matrix::dmat f1(3, 3);
    f1.set_value({2, 5, 1,
                  -4, 5, 0,
                  -1, 2, 1});
    auto [Q, R] = f1.houseHolder();

    std::print("value of Q\n{}\nvalue of R\n{}\n\n", Q, R);
}
void test_eigenvalues(){
    jf::matrix::dmat f2(3, 3);
    f2.set_value({2, 0, 1,
                -1, 4, -1,
                -1, 2, 0});
    auto eigs1 = f2.eigenvalues();
    
    std::print("the eigenvalues of f2 are\n");
    for(const auto& v : eigs1)
        std::print("{} - ", v);
    std::print("\n");
    
    jf::matrix::dmat f3(4, 4);
    f3.set_value({6, 2, 1, 3,
                 2, 5, 2, 1,
                 1, 2, 4, 2,
                3, 1, 2, 7});
    auto eigs2 = f3.eigenvalues();
    
    std::print("the eigenvalues of f3 are\n");
    for(const auto& v : eigs2)
        std::print("{} - ", v);
    
    jf::matrix::dmat f4(2, 2);
    f4.set_value({-6, 3,
                 4, 5});
    auto eigs3 = f4.eigenvalues();
    
    std::print("\nthe eigenvalues of f4 are\n");
    for(const auto& v : eigs3)
        std::print("{} - ", v);
    std::print("\nfim\n");
}
int main(){
    //test_householder();
    test_eigenvalues();
    return 0;
}