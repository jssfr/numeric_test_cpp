
#include "headers/std_input.hpp"
#include <fmt/core.h>
#include <numbers>
import math;
void test_bissec_falsepos(){

    auto a = 0.5;
    auto b = 1.0;

    auto [x] = jf::var::variables<1>();
    auto fx = x * x + jf::math::ln(x);

    auto result_bissec = jf::math::bissec(a, b, fx);
    auto result_falpos = jf::math::false_position(a, b, fx);

    fmt::println("the bissec for 0.05         = {}.", result_bissec);
    fmt::println("the false_position for 0.05 = {}.", result_falpos);
}

void test_fixed_point(){
    //

    auto [x] = jf::var::variables<1>();

    auto fx0 = jf::math::cbrt(x + 1);

    auto result = jf::math::fixed_point(fx0, 1.0);

    fmt::println("the fixed_point for cbrt(x + 1) = {}.", result);
    
    auto fx = x*x*x - x - 1;
    auto result_newton = jf::math::newton_method(fx, 1.0);

    fmt::println("the newton method for x^3 - x - 1 = {}.", result_newton);

    //  function for a pendulum:
    //  80 + 90*cos(pi/3 * t) = d
    //  if d = 10 cm 
    //  70 + 90 * cos(pi/3 * t) = 0
    //  what is the value of t for t0 = 4s and |f(x)| < 0.002?
    auto [t] = jf::var::variables<1>();
    auto pendulum_eq = 70.0 + 90.0*jf::math::cos(std::numbers::pi_v<double>/3.0 * t);

    auto r_pendulum = jf::math::newton_method(pendulum_eq, 4.0, 0.002);

    fmt::println("the newton method for pendulum_eq = {}.", r_pendulum);
}

void test_secant_method(){
    //volume of a sphere
    //(pi * h^2 * (3 * R - h) ) / 3 = v
    //if R = 1, v = 0.5
    //what is the value of h?
    //given h0 = 0.25 and h1 = 0.5
    //|f(x)| < 0.002

    auto [h] = jf::var::variables<1>();

    auto func_sphere = std::numbers::pi_v<double> / 3 * h*h*(3*1 - h) - 0.5;

    auto r_sphere = jf::math::secant_method(func_sphere, 0.25, 0.5, 0.002);
    fmt::println("the secant method for func_sphere = {}.", r_sphere);

    // x^3 - x - 1
    // x0 = 1.2 and x1 = 1.5
    // |f(x)| < 0.02
    auto [x] = jf::var::variables<1>();

    auto fx = x*x*x - x - 1;

    auto r_fx = jf::math::secant_method(fx, 1.2, 1.5, 0.02);
    fmt::println("the secant method for x^3 - x - 1 = {}.", r_fx);
}

auto main() -> int
{
    // test_bissec_falsepos();
    // test_fixed_point();
    test_secant_method();
    return 0;
}
