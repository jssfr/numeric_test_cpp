
#ifndef USING_IMPORT_STD_MOD
  #include "headers/std_input.hpp"
  #include <fmt/core.h>
  #include <print>
  #include <format>
  #include <numbers>
#endif

#ifdef USING_IMPORT_STD_MOD
  import std;
#endif

import math;

namespace fft = jf::math::fft;
// specialization for vector values in std::format/std::print
template<>
struct std::formatter<std::vector<std::complex<double>>>
{
    template<class ParseContext>
    constexpr ParseContext::iterator parse(ParseContext& ctx)
    {
        return ctx.begin();
    }

    template<class FmtContext>
    FmtContext::iterator format(const std::vector<std::complex<double>>& s, FmtContext& ctx) const
    {
        auto out = ctx.out();
        if(s.empty()) { std::format_to(out, "[  ]"); }
        else{
            std::format_to(out, "[ ");
            auto size_s = s.size();
            int inter{};
            for(const auto& cplx: s){
                std::format_to(out,"({0}, {1})", cplx.real(), cplx.imag());
                ++inter;
                (inter%2 != 0)? std::format_to(out, ", "): (inter==(int)size_s)?std::format_to(out, " ]"): std::format_to(out, "\n");
            } 
        }
        return out;
    }
};



void test_fft_evaluate(){

    // std::vector vec1{ 1.0, 5.0, 7.0, 1.0, 2.0, 3.0 };
    std::vector<std::complex<double>> vec1{ 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    auto answer = fft::evaluate(vec1);

    std::print("the fft for\n{}\n\n", vec1);
    std::print("= {}\n", answer);
}

    // template<typename PolynomialType>
    // auto fft_interpolate(const PolynomialType& p){
    //     using element_t = typename PolynomialType::value_type;
    //     std::size_t size = p.size();
    //
    //     if(size == 1){ return p;}
    //
    //     auto v_2pi_n = 2*std::numbers::pi_v<element_t> / (element_t)size;
    //
    //     auto w = [v_2pi_n](auto j){ return std::exp(- v_2pi_n * j); };
    //
    //     auto [Pe, Po] = fft::get_even_odd(p);
    //     auto Ye = fft_interpolate(Pe);
    //     auto Yo = fft_interpolate(Po);
    //
    //     PolynomialType y(size);
    //     auto n_2 = size / 2;
    //
    //     for(std::size_t j{}; j < n_2; ++j)
    //     {
    //         auto w_j = w(j) * Yo[j];
    //
    //         y[j] = Ye[j] + w_j;
    //         y[j + n_2] = Ye[j] - w_j;
    //     }
    //     auto v_1_n = 1.0 / (element_t)size;
    //     // for (auto& val : y) {
    //     //     val *= v_1_n;
    //     // }
    //     auto ifft = std::views::all(y) | std::views::transform([v_1_n](auto val){return val * v_1_n;})
    //                | std::ranges::to<std::vector>();
    //     return ifft;
    // } // end of fft interpolate

// void test_fft_interpolate(){
//     fmt::println("test_fft_evaluate called");
//
//     std::vector vec1{ 1.0, 5.0, 7.0, 1.0, 2.0, 3.0 };
//
//     auto answer = fft::fft_evaluate(vec1);
//
//     fmt::println("the fft for");
//     for (const auto& vec : vec1) {
//     fmt::print("{} ", vec);
//     }
//     fmt::println("=");
//    for (const auto& vec : answer) {
//     fmt::print("{} ", vec);
//    } 
//     fmt::println("");
//
//     auto ifft = fft::fft_interpolate(answer);
//
//     fmt::println("fft interpolate =");
//     for (const auto& vec : ifft) {
//         fmt::print("{} ", vec);
//     }
//     fmt::println("");
//
// }
void test_fft_derivative(){

    auto sin_f = [](auto x){ return std::sin(x); };
    auto cos_f = [](auto x){ return std::cos(x); };

    auto x = std::numbers::template pi_v<double> / 6;

    auto n_sin_7 = fft::derivative<7, 64>(sin_f, x);
    auto a_sin_7 = -cos_f(x);

    std::print("Numerical = {}\nAnalitic = {}\n\n", n_sin_7, a_sin_7);

}

void test_fft_derivative2(){

    auto sqrt_f = [](auto x){ return std::pow(x, 1.0/2); };
    auto sqrt_f3 = [](auto x){ return 1.0/2 * -1.0/2 * -3.0/2 * std::pow(x, -5.0/2 ); };

    auto x = 5.0;
    auto n_sqrt_3 = fft::derivative<3, 32>(sqrt_f, x);
    auto a_sqrt_3 = sqrt_f3(x);
    auto stencil = jf::math::diff_stencil<3>(sqrt_f, x);
    std::print("Numerical = {}\nAnalitic  = {}\nstencil   = {}\n\n", n_sqrt_3, a_sqrt_3, stencil);

}


/// jacobian matrix

// template<typename ContainerType, typename... FuncTypes, typename... ArgTypes, std::size_t Size = sizeof...(ArgTypes)>
// auto create_jacobian_matrix(const std::tuple<FuncTypes...>& functions, const std::tuple<ArgTypes...>& args, bool with_result = false)
// {
//     auto column_count = with_result? Size+1 : Size;
//
//     if constexpr(matrix::is_vector_v<ContainerType>){
//         ContainerType jacobian(Size * column_count);
//     }
//
//     jf::types::lambda_seq<Size>([&functions, &args, &jacobian](auto i){
//             auto& func = std::get<i.Index>(functions);
//
//             jf::types::lambda_seq<Size>([&](auto j){
//                     auto f_fixed = fix_arguments<j.Index>(func, args);
//
//                     auto pos = Size * i.Index + j.Index;
//
//                     jacobian[pos] = jf::math::diff_stencil(f_fixed, std::get<j.Index>(args));
//
//             });
//     });
//
//     return jacobian;
// }



auto main() -> int
{
   // test_fft_evaluate();
   // test_fft_interpolate();
   // test_fft_derivative();
   test_fft_derivative2();
    return 0;
}
