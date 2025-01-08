module;

#include<cmath>
#include<limits>

export module math:calculus;
import Config;
import :variable;

namespace{
    auto h = std::sqrt(std::numeric_limits<double>::epsilon());
}

namespace jf::math{
    template<class FuncTypes, jf::types::floating_c Number>
    auto diff(FuncTypes&& func, Number x)
    {
        return (func(x + h) - func(x - h)) / (2.0 * h);
    }

    template<class FuncTypes>
    auto diff(FuncTypes&& func)
    {
        return [func](auto arg){ return (func(arg + h) - func(arg - h)) / (2.0 * h);};
    }
} //namespace jf::math

