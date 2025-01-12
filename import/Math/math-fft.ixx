module;
// cmath limits complex etc.
#include "std_input.hpp"

export module math:Fft;
// import :Config;

// @link https://youtu.be/h7apO7q16V0?si=sm4suPYBzwqJESdC
// @link https://en.wikipedia.org/wiki/Fast_Fourier_transform

// @link https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm

 /*! TODO:  
  *  \todo ajust for complex
  *  \todo implement cooley turkey fft algorithm
 */
namespace jf::math {

    template<typename PolynomialType>
   auto get_even_odd(const PolynomialType& p)
   {
       std::size_t size = p.size();
       auto n_2 = size / 2;

       PolynomialType even_poly, odd_poly;
       even_poly.reserve(n_2);
       odd_poly.reserve(n_2);

       for (std::size_t i{}; i < size; ++i) {
            if(i % 2 == 0){ even_poly.emplace_back(p[i]);}
            else {odd_poly.emplace_back(p[i]);}
       }
        
       return std::pair {even_poly, odd_poly};
   }

    export template<typename PolynomialType>
    auto fft_evaluate(const PolynomialType& p){
        using element_t = typename PolynomialType::value_type;
        std::size_t size = p.size();

        if(size == 1){ return p;}

        auto v_2pi_n = 2*std::numbers::pi_v<element_t> / (element_t)size;

        auto w = [v_2pi_n](auto j){ return std::exp( v_2pi_n * j); };

        auto [Pe, Po] = get_even_odd(p);
        auto Ye = fft_evaluate(Pe);
        auto Yo = fft_evaluate(Po);

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

    export template<typename PolynomialType>
    auto fft_interpolate(const PolynomialType& p){
        using element_t = typename PolynomialType::value_type;
        std::size_t size = p.size();

        if(size == 1){ return p;}

        auto v_2pi_n = 2*std::numbers::pi_v<element_t> / (element_t)size;

        auto w = [v_2pi_n](auto j){ return std::exp(- v_2pi_n * j); };

        auto [Pe, Po] = get_even_odd(p);
        auto Ye = fft_interpolate(Pe);
        auto Yo = fft_interpolate(Po);

        PolynomialType y(size);
        auto n_2 = size / 2;

        for(std::size_t j{}; j < n_2; ++j)
        {
            auto w_j = w(j) * Yo[j];

            y[j] = Ye[j] + w_j;
            y[j + n_2] = Ye[j] - w_j;
        }
        auto v_1_n = 1.0 / (element_t)size;

        std::ranges::for_each(y, [v_1_n](auto& val){val *= v_1_n;});

        return y;
    } // end of fft interpolate
      
      // cooley turkey algorithm

} // namespace jf::math
