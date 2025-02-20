
#include <print>
#include <cmath>
#include <numbers>
// #include <fmt/core.h>
// import Matrix;
import math;

void test_multip(){
    jf::matrix::mat<int> m1{2, 2};
    jf::matrix::mat<int> m2{2, 1};

    m1(0, 0) = 1; m1(0, 1) = 2;
    m1(1, 0) = 5; m1(1, 1) = -3;

    m2(0, 0) = 4;
    m2(1, 0) = 8;
    jf::matrix::mat<int> m3 = m1 * m2;

    std::print("{} - {}\n", m3(0, 0), m3(1, 0));
    // std::print("value of m3 = \n{}\n\n", m3);
}

void test_better(){
    jf::matrix::imat m1{3, 3};
    jf::matrix::imat m2{3, 1};

    m1.set_value({1, 2, 5,
                  5, 8, -3,
                  0, 2, 9});

    m2.set_value({4, 8, 7});

    jf::matrix::imat m3 = m1 * m2;

    std::print("{} - {} - {}\n", m3(0, 0), m3(1, 0), m3(2, 0));
    std::print("{}\n", m3);
    // std::print("value of m3 = \n{}\n\n", m3);
}
void test_base_points()
{
    auto [rho, phi1, z] = jf::math::point_cart_cyl(-2, 6, 3);

    auto [r, theta, phi2] = jf::math::point_cart_spher(-2, 6, 3);

    std::println("values: rho = {}, phi1 = {}, phi2 = {}, z = {}, r = {}, theta = {}.", rho, phi1, phi2, z, r, theta);

    std::print("cos(45) = {}, acos(0.707) = {}\ncos(180) = {}, cos(pi) = {}\n\n", 
            std::cos(45*std::numbers::pi/180), std::acos(0.707*std::numbers::pi/180), std::cos(180*std::numbers::pi/180), std::cos(std::numbers::pi));
}

void test_base_change1()
{
    auto [r, theta, phi] = jf::math::point_cart_spher(-2, 6, 3);

    auto [Bx, By, Bz] = jf::math::cart_cyl(6, 1, 0, phi);

    std::println("values: Bx = {}, By = {}, Bz = {}.", Bx, By, Bz);

    auto [Ax, Ay, Az] = jf::math::cart_spher(6, 1, 0, theta, phi);
    std::println("values: Ax = {}, Ay = {}, Az = {}.", Ax, Ay, Az);

    auto [x, y, z] = jf::var::variables<3>();

    auto [Cx, Cy, Cz] = jf::math::cart_cyl({y, x + z, 0*z}, {-2, 6, 3});

    std::println("values: Cx = {}, Cy = {}, Cz = {}.", Cx, Cy, Cz);

    // auto ax = jf::math::sqrt(x^2 + y^2) / jf::math::sqrt(x^2 + y^2 + z^2);
    auto ax = jf::math::sqrt(x*x + y*y) / jf::math::sqrt(x*x + y*y + z*z);
    auto ay = 0*x;
    // auto az = y*z / jf::math::sqrt(x^2 + y^2 + z^2);
    auto az = -1*y*z / jf::math::sqrt(x*x + y*y + z*z);

    auto [Dx, Dy, Dz] = jf::math::cart_cyl({ax, ay, az}, {0, -4, 3});

    std::println("values: Dx = {}, Dy = {}, Dz = {}.", Dx, Dy, Dz);
}

void test_base_change2()
{
    auto [r, theta, phi] = jf::math::point_cart_spher(-3, 4, 0);

    auto [Bx, By, Bz] = jf::math::spher_cart(10/r, r*std::cos(theta), 1, theta, phi);

    std::println("values: Bx = {}, By = {}, Bz = {}.", Bx, By, Bz);

    auto [r2, theta2, phi2] = jf::math::point_cyl_spher(5, std::numbers::pi/2, -2);
    auto [Ax, Ay, Az] = jf::math::spher_cyl(10/r2, r2*std::cos(theta2), 1, theta2);
    std::println("values: Ax = {}, Ay = {}, Az = {}.", Ax, Ay, Az);

}
int main()
{
    std::println("test matrix.");
    // test_multip();
    // test_better()/* ; */
   test_base_change1();
   // test_base_change2();
    // test_base_points();
    return 0;
}
