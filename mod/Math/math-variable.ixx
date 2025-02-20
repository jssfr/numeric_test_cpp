module;

#include "std_input.hpp" 
#include<tuple>
export module math:Variable;

import Config;

export namespace jf::var {

    // creating lambdas 
    // @link https://www.youtube.com/@HomoSiliconiens
template <auto N, auto n>
auto create_oprvar() noexcept {
    auto f = [](auto... args)
        requires(sizeof...(args) == N)
    { return std::get<n>(std::forward_as_tuple(args...)); };

    return 1 * f;
}

template <auto N, auto... n>
auto variables() noexcept {
    if constexpr (sizeof...(n) == 0) {
        return []<typename T, T... k>(std::integer_sequence<T, k...>) {
            return std::tuple{create_oprvar<N, k>()...};
        }(std::make_integer_sequence<int, N>{});
    } else {
        return std::tuple{create_oprvar<N, n>()...};
    }
}

//======================
// for convienience
// just a test, better implementation is needed
template<typename T, typename Arg, typename... Args>
auto validate(T&& func, const Arg& x1, const Args&... x2) {
    return func(x1, x2...);
}

// operator for operations with lambdas (lambda1 + lambda2)

template <class F, class G>
auto operator+(F&& f, G&& g) noexcept {
    return [=]<typename T, typename... ts>(T x0, ts... args) {
        return (validate(f, x0, args...) + validate(g, x0, args...));
    };
}

template <class F, jf::types::number_c Gtype>
auto operator+(F&& f, Gtype&& g) noexcept {
    return [=]<typename T, typename... ts>(T x0, ts... args) {
        return (validate(f, x0, args...) + g);
    };
}

template <jf::types::number_c Ftype, class G>
auto operator+(Ftype&& f, G&& g) noexcept {
    return [=]<typename T, typename... ts>(T x0, ts... args) {
        return (f + validate(g, x0, args...));
    };
}

// operator for operations with lambdas (lambda1 - lambda2)
template <class F, class G>
auto operator-(F&& f, G&& g) noexcept {
    return [=]<typename T, typename... ts>(T x0, ts... args) {
        return (validate(f, x0, args...) - validate(g, x0, args...));
    };
}

template <class F, jf::types::number_c Gtype>
auto operator-(F&& f, Gtype&& g) noexcept {
    return [=]<typename T, typename... ts>(T x0, ts... args) {
        return (validate(f, x0, args...) - g);
    };
}

template <jf::types::number_c Ftype, class G>
auto operator-(Ftype&& f, G&& g) noexcept {
    return [=]<typename T, typename... ts>(T x0, ts... args) {
        return (f - validate(g, x0, args...));
    };
}

// operator for operations with lambdas (lambda1 * lambda2)
template <class F, class G>
auto operator*(F&& f, G&& g) noexcept {
    return [=]<typename T, typename... ts>(T x0, ts... args) {
        return (validate(f, x0, args...) * validate(g, x0, args...));
    };
}

template <class F, jf::types::number_c Gtype>
auto operator*(F&& f, Gtype&& g) noexcept {
    return [=]<typename T, typename... ts>(T x0, ts... args) {
        return (validate(f, x0, args...) * g);
    };
}

template <jf::types::number_c Ftype, class G>
auto operator*(Ftype&& f, G&& g) noexcept {
    return [=]<typename T, typename... ts>(T x0, ts... args) {
        return (f * validate(g, x0, args...));
    };
}
// operator for operations with lambdas (lambda1 / lambda2)
template <class F, class G>
auto operator/(F&& f, G&& g) noexcept {
    return [=]<typename T, typename... ts>(T x0, ts... args) {
        return (validate(f, x0, args...) / validate(g, x0, args...));
    };
}

template <class F, jf::types::number_c Gtype>
auto operator/(F&& f, Gtype&& g) noexcept {
    return [=]<typename T, typename... ts>(T x0, ts... args) {
        return (validate(f, x0, args...) / g);
    };
}

template <jf::types::number_c Ftype, class G>
auto operator/(Ftype&& f, G&& g) noexcept {
    return [=]<typename T, typename... ts>(T x0, ts... args) {
        return (f / validate(g, x0, args...));
    };
}

// template <class F, jf::types::number_c Gtype>
// auto operator^(F&& f, Gtype&& g) noexcept {
//     return [=]<typename T, typename... ts>(T x0, ts... args) {
//         return std::pow(validate(f, x0, args...), g);
//     };
// }
// operator to add some value to a function or function operations (lambdas)
template <typename FuncType, typename NATType>
auto operator,(FuncType&& f, NATType&& args) noexcept {
    return std::apply(f, args);
}
}  // namespace jf::var
