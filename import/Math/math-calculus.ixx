module;

#include<cmath>
#include<limits>

export module math:Calculus;
import Config;
import :Variable;

namespace{
    // constexpr auto h1 = 2.0001e-4;
    constexpr double h1 = 0.001;
}

export namespace jf::math{
    inline constexpr auto diff = []<jf::types::floating_c Number>(auto func, Number x)
    {
        auto h = std::sqrt(std::numeric_limits<double>::epsilon());
        return (func(x + h) - func(x - h)) / (2.0 * h);
    };

    template<class FuncType, jf::types::floating_c Number>
    auto five_point_stencil(FuncType&& func, Number x){
        auto fn1 = func(x + 2.0 * h1);
        auto fn2 = func(x + h1);
        auto fn3 = func(x - h1);
        auto fn4 = func(x - 2.0 * h1);

        return (8.0*(fn2 -fn3) + fn4 - fn1) / (12.0 * h1);
    }
} //namespace jf::math

