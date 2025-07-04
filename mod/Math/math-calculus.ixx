module;

#ifndef USING_IMPORT_STD_MOD
#  include <cmath>
#  include <limits>
#endif

export module math:Calculus;
#ifdef USING_IMPORT_STD_MOD
import std;
#endif
import Config;
import :Variable;

namespace ncnt
{
// constexpr auto h1 = 2.0001e-4;

struct constants
{
  constexpr static double h1 = 0.001;
  constexpr static double h2 = 0.001;
  constexpr static double h3 = 0.01;
  constexpr static double h4 = 0.01;
} cnt;

template<auto Order> auto delta()
{
  // constants cnt;
  if (Order == 1) return cnt.h1;
  if (Order == 2) return cnt.h2;
  if (Order == 3) return cnt.h3;
  if (Order == 4) return cnt.h4;
}
}  // namespace ncnt

export namespace jf::math
{

inline constexpr auto diff = []<jf::types::floating_c Number>(auto func, Number x)
{
  auto h = 1.25e-6;
  return (func(x + h) - func(x - h)) / (2.0 * h);
};

// f p s backward
//  ...

//  f p s forward
//  ...

// f p s central difference method
template<auto DerivOrder = 1, class FuncType, jf::types::floating_c Number>
auto five_point_stencil(FuncType&& func, Number x)
{
  auto h = ncnt::delta<DerivOrder>();
  double fn[5];
  fn[0] = func(x);
  fn[1] = func(x + 2 * h);
  fn[2] = func(x + h);
  fn[3] = func(x - h);
  fn[4] = func(x - 2 * h);

  if (DerivOrder == 1) return (8 * (fn[2] - fn[3]) + fn[4] - fn[1]) / (12 * h);
  if (DerivOrder == 2)
    return (16 * (fn[2] + fn[3]) - 30 * fn[0] - (fn[1] + fn[4])) / (12 * std::pow(h, 2.0));
  if (DerivOrder == 3) return (2 * (fn[3] - fn[2]) + fn[1] - fn[4]) / (2 * std::pow(h, 3.0));
  if (DerivOrder == 4) return (6 * fn[0] + fn[1] + fn[4] - 4 * (fn[2] + fn[3])) / std::pow(h, 4.0);

  return static_cast<Number>(0);  // seven point needed
}

template<auto DerivOrder = 1, class FuncType, jf::types::floating_c Number>
auto diff_stencil(FuncType&& func, Number x)
{
  if (DerivOrder < 5) return five_point_stencil<DerivOrder>(func, x);

  return static_cast<Number>(0);
}

}  // namespace jf::math
