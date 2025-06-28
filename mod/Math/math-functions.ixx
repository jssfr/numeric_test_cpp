module;

#ifndef USING_IMPORT_STD_MOD
#  include <cmath>
#endif

export module math:Functions;
#ifdef USING_IMPORT_STD_MOD
import std;
#endif
import Config;
import :Variable;

namespace jf::var
{
// using namespace jf::var;
template<class FuncType> inline constexpr auto ln(FuncType&& func)
{
  return [=]<typename T, typename... ts>(T x0, ts... args) { return std::log(func(x0, args...)); };
}
template<class FuncType> inline constexpr auto cbrt(FuncType&& func)
{
  return [=]<typename T, typename... ts>(T x0, ts... args) { return std::cbrt(func(x0, args...)); };
}

template<class FuncType> inline constexpr auto sqrt(FuncType&& func)
{
  return [=]<typename T, typename... ts>(T x0, ts... args) { return std::sqrt(func(x0, args...)); };
}

template<class FuncType> inline constexpr auto sin(FuncType&& func)
{
  return [=]<typename T, typename... ts>(T x0, ts... args) { return std::sin(func(x0, args...)); };
}

template<class FuncType> inline constexpr auto cos(FuncType&& func)
{
  return [=]<typename T, typename... ts>(T x0, ts... args) { return std::cos(func(x0, args...)); };
}

}  // namespace jf::var

export namespace jf::math
{

template<class FuncType> inline constexpr auto ln(FuncType&& func)
{
  return jf::var::ln(func);
}

template<class FuncType> inline constexpr auto cbrt(FuncType&& func)
{
  return jf::var::cbrt(func);
}

template<class FuncType> inline constexpr auto sqrt(FuncType&& func)
{
  return jf::var::sqrt(func);
}

template<class FuncType> inline constexpr auto sin(FuncType&& func)
{
  return jf::var::sin(func);
}

template<class FuncType> inline constexpr auto cos(FuncType&& func)
{
  return jf::var::cos(func);
}

}  // namespace jf::math
