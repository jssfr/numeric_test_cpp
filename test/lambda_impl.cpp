#include "std_input.hpp"
#include <print>
import Config;


namespace jf {
    namespace hidden{
        template<typename T>
            struct st_is_tuple_or_array{
                static constexpr auto value = false;
            };
        template<typename ...Types>
            struct st_is_tuple_or_array<std::tuple<Types...>>
            {
                static constexpr auto value = true;
            };
        template<typename T, size_t Nm>
            struct st_is_tuple_or_array<std::array<T, Nm>>
            {
                static constexpr auto value = true;
            };
    } // namespace hidden
    template<typename T>
        constexpr auto is_tuple_or_array_v = hidden::st_is_tuple_or_array<jf::types::remove_cv_ref_t<T> >::value;
    template<typename T>
        concept tuple_or_array_c = is_tuple_or_array_v<T>;
    
} // namespace jf



// template< typename TypeContainer, auto N = std::tuple_size_v<TypeContainer>, typename Indx = std::make_index_sequence<N>>
// auto lambda_seq( auto&& func) -> decltype(auto)
// {
//     return func(Indx{});
//     
// }


auto process_lambdas_add(auto&& ...args) -> decltype(auto)
{
    auto process_arg = [](auto& arg){
        using arg_t = jf::types::remove_cv_ref_t<decltype(arg)>;

        if constexpr(jf::tuple_or_array_c<arg_t>)
        {
            constexpr size_t N = std::tuple_size_v<arg_t>;
            auto process = [&arg]<auto ...i>(std::index_sequence<i...>){
                return process_lambdas_add(std::get<i>(arg)...);
            };

            //return lambda_seq<arg_t>(process);
            return process(std::make_index_sequence<N>{});
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

    std::println("value of result = {}", result);
}

auto main() -> int
{
    std::println("process lambda impl");
    test_process_lambdas_add();
    return 0;
}
