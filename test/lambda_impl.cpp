#include "std_input.hpp"
#include <print>
#include <iostream>
//#include <fmt/core.h>

import Config;


// namespace jf {
//     namespace hidden{
//         template<typename T>
//             struct st_is_tuple_or_array{
//                 static constexpr auto value = false;
//             };
//         template<typename ...Types>
//             struct st_is_tuple_or_array<std::tuple<Types...>>
//             {
//                 static constexpr auto value = true;
//             };
//         template<typename T, size_t Nm>
//             struct st_is_tuple_or_array<std::array<T, Nm>>
//             {
//                 static constexpr auto value = true;
//             };
//     } // namespace hidden
//     template<typename T>
//         constexpr auto is_tuple_or_array_v = hidden::st_is_tuple_or_array<std::remove_cvref_t<T> >::value;
//     template<typename T>
//         concept tuple_or_array_c = is_tuple_or_array_v<T>;
//     
// } // namespace jf



template< typename TypeContainer, auto  endVal = std::tuple_size_v<TypeContainer>>
auto lambda_seq( auto&& func) -> decltype(auto)
{
    return func(jf::types::make_sequence<decltype(endVal), 0, endVal, 1>{});
    
}

template< auto  endVal>
auto lambda_seq( auto&& func) -> decltype(auto)
{
    return func(jf::types::make_sequence<decltype(endVal), 0, endVal, 1>{});
    
}
template< auto startVal, auto  endVal, auto stepVal>
auto lambda_seq( auto&& func) -> decltype(auto)
{
    return func(jf::types::make_sequence<decltype(endVal), startVal, endVal, stepVal>{});
    
}

auto process_lambdas_add(auto&& ...args) -> decltype(auto)
{
    auto process_arg = [](auto& arg){
        using arg_t = std::remove_cvref_t<decltype(arg)>;

        if constexpr(jf::types::tuple_or_array_c<arg_t>)
        {
            constexpr size_t N = std::tuple_size_v<arg_t>;
            auto process = [&arg]<auto ...i>(std::index_sequence<i...>){
                return process_lambdas_add(std::get<i>(arg)...);
            };

            return lambda_seq<arg_t>(process);
            // return process(std::make_index_sequence<N>{});
        }
        else {
            return arg;
        }
    };

    return (process_arg(args) + ...);
}

auto test_process_lambdas_add() -> void
{

    auto result = process_lambdas_add(1, 2, 3, std::array{4, 5, 6}, 7, std::tuple{8, 9});

    std::print("value of result = {0}\n\n", result);
    //std::cout<<result<<'\n';
}
//////////////////////
auto process_lambdas2(auto&& point, auto&&...funcs) -> decltype(auto){
    constexpr auto size = sizeof...(funcs);
    
    auto process_funcs = 
        [&point, fgroup = std::forward_as_tuple(funcs...)]<auto ...i>(std::index_sequence<i...> sequ)
    {
        auto call = [&point](auto&& func){
            if constexpr(jf::types::arithmetic_c<decltype(func)>)
                return func;
            else
             return std::apply(func, point);
        };
        return std::tuple{call(std::get<i>(fgroup)) ... };
    };
    
    return lambda_seq<size>(process_funcs);
}

void test_process_lambdas2()
{
    auto func1 = [](auto x, auto y, auto z){
        return (x + y) / 2.0 * z * z;
    };
    // auto func2 = [](auto x, auto y, auto z){
    //     return (x - y) / 2.0 * z;
    // };
    auto func2 = 4.5;
    auto func3 = [](auto x, auto y, auto z){
        return (x - z) / 4.0 * y;
    };
    auto func4 = [](auto x, auto y, auto z){
        return (y + z) / 5.0 * x;
    };
    
    auto [re1, re2, re3, re4] = process_lambdas2(std::tuple{1.0, 2.0, 3.0}, func1, func2, func3, func4);

    std::print("result0 = {0}\nresult1 = {1}\nresult2 = {2}\nresult3 = {3}\n\n", re1, re2, re3, re4);
}

auto process_lambdas_as_func_arg_pairs(auto ...args) -> decltype(auto)
    requires (sizeof...(args) % 2 == 0)
{
    constexpr auto N = sizeof...(args);

    auto process_func = [arguments = std::make_tuple(args...)]<auto ...i>(std::index_sequence<i...> seq){
        return std::tuple{ std::apply(std::get<i>(arguments), std::get<i+1>(arguments)) ... };
    };

    return lambda_seq<0, N, 2>(process_func);
}
auto test_process_lambdas_as_func_arg_pairs() -> void
{
    auto func1 = [](auto x, auto y, auto z){
        return (x + y) / 2.0 * z * z;
    };
    auto func2 = [](auto x, auto y, auto z){
        return (x - y) / 2.0 * z;
    };
    auto func3 = [](auto x, auto y, auto z){
        return (x - z) / 4.0 * y;
    };
    auto func4 = [](auto x, auto y, auto z){
        return (y + z) / 5.0 * x;
    };
    
    auto [re1, re2, re3, re4] = process_lambdas_as_func_arg_pairs(
             func1, std::tuple{1.0, 2.0, 3.0}, func2, std::tuple{4.0, 2.0, 3.0}, func3, std::tuple{1.0, 2.5, 3.0},
             func4, std::tuple{1.0, 2.0, 3.2});

    std::print("point as tuple:\nresult0 = {0}\nresult1 = {1}\nresult2 = {2}\nresult3 = {3}\n", re1, re2, re3, re4);
    
    std::println("point as array:");
    auto [ra1, ra2, ra3, ra4] = process_lambdas_as_func_arg_pairs(
             func1, std::array{1.0, 2.0, 3.0}, func2, std::array{4.0, 2.0, 3.0}, func3, std::array{1.0, 2.5, 3.0},
             func4, std::array{1.0, 2.0, 3.2});

 
    std::print("result0 = {0}\nresult1 = {1}\nresult2 = {2}\nresult3 = {3}\n\n", ra1, ra2, ra3, ra4);
   std::print("result for lambda_seq(func1, sequence<1, 5, 4>) = {0}\n\n",jf::types::lambda_seq(func1, jf::types::sequence<1, 5, 4>{}) );
   auto lamb = [&func1]<auto ... i>(jf::types::sequence<i...>){
       return std::apply(func1, std::tuple{i...});
   };
   auto enx = jf::types::lambda_seq(lamb, jf::types::sequence<1, 5, 4>{});

   std::print("Result for lamb <1, 5, 4> = {0}\n\n", enx);
}

auto main() -> int
{
    std::println("process lambda impl");
    test_process_lambdas_add();
    std::println("process lambda 2");
    test_process_lambdas2();
    std::println("process lambda as arg pairs");
    test_process_lambdas_as_func_arg_pairs();
    return 0;
}
