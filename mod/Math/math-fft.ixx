module;
// cmath limits complex etc.
#include "std_input.hpp"
#include<numbers>

export module math:Fft;
import :ncrnpr;
import Config;

// @link https://youtu.be/h7apO7q16V0?si=sm4suPYBzwqJESdC
// @link https://en.wikipedia.org/wiki/Fast_Fourier_transform

// @link https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
// @link https://www.youtube.com/@HomoSiliconiens

namespace jf::math::fft {

    template<typename PolynomialType>
   auto get_even_odd(const PolynomialType& poly)
   {
       std::size_t size = poly.size();
       auto n_2 = size / 2;

       PolynomialType even_poly;
       PolynomialType odd_poly;

       even_poly.reserve(n_2);
       odd_poly.reserve(n_2);

       for (std::size_t i{}; i < size; ++i) {
            if(i % 2 == 0){ even_poly.emplace_back(poly[i]);}
            else {odd_poly.emplace_back(poly[i]);}
       }
        
       return std::pair {even_poly, odd_poly};
   }

    export template<typename PolynomialType>
    constexpr auto evaluate(const PolynomialType& poly){
        using cmplx = typename PolynomialType::value_type;
        using element_t = typename cmplx::value_type;
        std::size_t size = poly.size();

        if(size == 1){ return poly;}

        auto v_2pi_n = 2*std::numbers::template pi_v<element_t> / static_cast<element_t>(size);

        auto w = [v_2pi_n](auto j){ return std::exp( std::complex<element_t>{ element_t{}, v_2pi_n * j});};

        auto [Pe, Po] = get_even_odd(poly);
        auto Ye = evaluate(Pe);
        auto Yo = evaluate(Po);

        PolynomialType y(size);
        auto n_2 = size / 2;

        for(std::size_t j{}; j < n_2; ++j)
        {
            auto w_j = w(j) * Yo[j];

            y[j] = Ye[j] + w_j;
            y[j + n_2] = Ye[j] - w_j;
        }

        return y;
    } // end of fft evaluate
    
    template<typename PolynomialType>
    auto fft_interpolate_impl(const PolynomialType& poly){
        using compx_t   = typename PolynomialType::value_type;
        using element_t = typename compx_t::value_type;

        std::size_t size = poly.size();

        if(size == 1){ return poly; }

        auto v_2pi_n = 2*std::numbers::template pi_v<element_t> / static_cast<element_t>(size);

        auto w = [v_2pi_n](auto j){ return std::exp(compx_t{element_t{}, - v_2pi_n * j} ); };

        auto [Pe, Po] = get_even_odd(poly);
        auto Ye = fft_interpolate_impl(Pe);
        auto Yo = fft_interpolate_impl(Po);

        PolynomialType y(size);
        auto n_2 = size / 2;

        for(std::size_t j{}; j < n_2; ++j)
        {
            auto w_j = w(j) * Yo[j];

            y[j] = Ye[j] + w_j;
            y[j + n_2] = Ye[j] - w_j;
        }

        return y;
    } // end of fft interpolate
    
    export template<typename PolynomialType>
    auto interpolate(const PolynomialType& poly){
    
        using type_t = typename PolynomialType::value_type;
        auto P = fft_interpolate_impl(poly);

        auto v_1_n = 1.0 / static_cast<type_t>(poly.size());
        std::ranges::for_each(P, [v_1_n](auto& val){val *= v_1_n;});
       return P;
    }

    // fft derivative
    template<std::size_t id, int Order, int N,
           typename FuncType, jf::types::number_c ArgType, jf::types::number_c... ArgTypes,
           auto count = jf::types::param_count_v<FuncType>>
    auto fft_derivative_impl(FuncType&& func, ArgType arg, ArgTypes... args)
    {
        using Cplx = std::complex<ArgType>;

        auto w = [](auto j)
        {
            return std::exp(Cplx{ArgType{}, 2 * std::numbers::template pi_v<ArgType> / N * j});
        };

        auto ff = [&func, &w, argpack = std::tuple{arg, args...}](auto x)
        {
            auto set_arg = [&]<auto i>(jf::types::sequence<i>) -> Cplx
            {
                if constexpr(i == id)
                    return std::get<i>(argpack) + x;
                else
                    return std::get<i>(argpack);
            };

            return jf::types::lambda_seq<count>(
                    [&]<auto...i>(jf::types::sequence<i...>){
                       return func(set_arg(jf::types::sequence<i>{} ...));
                        }
                    );
        };

        std::vector<Cplx> ValuesOf_ff(N);

        for(int j{}; j < N; ++j)
        {
            ValuesOf_ff[j] = ff( w(j) );
        }

        auto coeffs = interpolate(ValuesOf_ff);

        return coeffs[Order] * jf::ncrnpr::fact(static_cast<ArgType>(Order));
    }

    template<std::size_t id, std::size_t Order, int N, typename FuncType>
        auto fft_derivative(FuncType func)
        {
            return [func]<typename ArgType, typename... ArgTypes>
                (ArgType arg, ArgTypes... args)
            {
                using f_return_type = decltype(func(arg, args...));

                auto df = fft_derivative_impl<id, Order, N>(func, arg, args...);

                if constexpr (jf::types::complex_c<f_return_type> || jf::types::complex_c<ArgType> )
                    return df;
                else
                    return df.real();
            };
        }

        export template<std::size_t Order = 1, int N = 32, std::size_t id = 0,
            typename FuncType, typename... ParamTypes>
            requires (sizeof...(ParamTypes) == jf::types::param_count_v<FuncType>)
        auto derivative(FuncType func, ParamTypes... params)
        {
            auto df = fft_derivative<id, Order, N>(func);
            return df(params...);
        }


} // namespace jf::math::fft
