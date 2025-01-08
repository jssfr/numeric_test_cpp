

#ifndef _JF_NCRNPR_HPP
#define _JF_NCRNPR_HPP

#ifdef _MSVC_LANG
#if _MSVC_LANG < 201703L
#error This library requires c++17 or later
#endif
#else
#if __cplusplus < 201703L
#error This library requires c++17 or later
#endif
#endif  // _MSVC_LANG

/*
    INCLUDES HERE
*/

#include "std_input.hpp"

#ifdef max
#define _JF_NCRNPR_MAX_DEFINED
#pragma push_macro("max")
#undef max
#endif

#ifdef min
#define _JF_NCRNPR_MIN_DEFINED
#pragma push_macro("min")
#undef min
#endif

namespace jf {
namespace ncrnpr {
template <typename T>
T gcd(T a, T b) {
    if (b > a) std::swap(a, b);

    T r;
    while ((r = a % b) != 0) {
        a = b;
        b = r;
    }
    return b;
}

// least common multiple
template <typename T>
T lcm(T a, T b) {
    T G = gcd(a, b);
    return (a / G) * b;
}

template <typename T>
class UnsignedRational {
   public:
    using value_t = std::make_unsigned_t<T>;

   private:
    value_t m_n{0};  // numerator
    value_t m_d{1};  // denominator

    void reduce() {
        value_t G = gcd(this->m_n, this->m_d);

        this->m_n /= G;
        this->m_d /= G;
    }

   public:
    UnsignedRational() {}
    UnsignedRational(T n, T d) : m_n(n), m_d(d) { reduce(); }
    UnsignedRational(T n, T d, bool) : m_n(n), m_d(d) {}

    UnsignedRational(const UnsignedRational& rhs) = default;
    UnsignedRational& operator=(const UnsignedRational& rhs) = default;

    UnsignedRational(UnsignedRational&& rhs) = default;
    UnsignedRational& operator=(UnsignedRational&& rhs) = default;

    ~UnsignedRational() = default;

    template <typename C>
    operator C() const {
        return C((double)this->m_n / (double)this->m_d);
    }

    void operator*=(const UnsignedRational& r) {
        UnsignedRational r1(this->m_n, r.m_d);
        UnsignedRational r2(r.m_n, this->m_d);

        this->m_n = r1.m_n * r2.m_n;
        this->m_d = r1.m_d * r2.m_d;
    }

    void operator+=(const UnsignedRational& r) {
        T n = this->m_n, d = this->m_d;
        T L = lcm(m_d, r.m_d);

        this->m_n = n * L / d + r.m_n * (L / r.m_d);
        this->m_d = L;

        reduce();
    }
};

template <typename T>
T fact(T n) {
    T factorial{n};
    for (auto i = n - 1; i > 0; --i) factorial *= i;
    return factorial;
}

template <typename T>
T nCr(T n, T r) {
    if (r > n)
        return 0;
    else if (r == 0 || n == r)
        return 1;
    else if (r == (n - 1) || r == 1)
        return n;
    else {
        if (r > (n - r)) r = n - r;

        UnsignedRational<T> rlt{n, r};

        for (T i = 1; i < r; ++i) rlt *= UnsignedRational<T>(n - i, r - i);

        return (T)rlt;
    }
}

template <typename T>
T nPr(T n, T r) {
    if (r > n)
        return 0;
    else if (r == 0)
        return 1;
    else {
        T count{n};
        do {
            if (r == 1)
                break;
            else if (r == n || r == (n - 1)) {
                count = fact(n);
                break;

            } else {
                --n;
                --r;
                count *= n;
            }
        } while (true);

        return count;
    }
}

template <typename T>
T combi_double(T n, T r) {
    if (r > n)
        return 0;
    else if (r == 0 || n == r)
        return 1;
    else if (r == (n - 1) || r == 1)
        return n;
    else {
        if (r > (n - r)) r = n - r;

        double rlt{n, r};

        for (T i = 1; i < r; ++i) rlt *= double(n - i) / double(r - i);

        return T(rlt);
    }
}

// enum combination

template <class set_of_sets_t = std::vector<std::vector<int>>>
void enum_comb(set_of_sets_t&& SS, size_t n, size_t r, size_t mth) {
    std::vector<int> S(n);

    std::generate(S.begin(), S.end(), [n = size_t{}]() mutable { return ++n; });
    std::vector<int> R;
    do {
        if (r == 0 && !R.empty()) {
            SS.emplace_back(R);
            break;
        } else if (r == S.size()) {
            R.insert(R.end(), std::make_move_iterator(S.begin()), std::make_move_iterator(S.end()));
            SS.emplace_back(R);
            break;
        }
        // nCr = n-1Cr-1 + n-1Cr
        auto n_1_c_r_1 = nCr(S.size() - 1, r - 1);
        if (mth < n_1_c_r_1) {
            R.push_back(S.front());
            S.erase(S.begin());
            --r;
        } else {
            S.erase(S.begin());

            mth -= n_1_c_r_1;
        }

    } while (true);
}

template <class set_t = std::vector<int>>
void enum_permu(set_t& R, set_t S, size_t r, size_t mth) {
    // https://youtu.be/gixbV9dMlRg?si=yXEF7Jw8IICofgAg
    do {
        if (r == 0 || S.empty()) { break; }

        auto n = S.size();

        auto n_1_p_r_1 = jf::ncrnpr::nPr(n - 1, r - 1);

        auto quotient = mth / n_1_p_r_1;  // index of S

        R.push_back(S[quotient]);

        S.erase(S.begin() + quotient);

        mth %= n_1_p_r_1;
        --r;

    } while (true);
}
auto enum_permu(size_t n, size_t r, size_t mth) {
    std::vector<int> S{};

    std::generate_n(std::back_inserter(S), n, [count = 0]() mutable { return count++; });

    std::vector<int> R{};
    do {
        if (r == 0 || S.empty()) { break; }

        auto n = S.size();

        auto n_1_p_r_1 = jf::ncrnpr::nPr(n - 1, r - 1);

        auto quotient = mth / n_1_p_r_1;  // index of S

        R.push_back(S[quotient]);

        S.erase(S.begin() + quotient);

        mth %= n_1_p_r_1;
        --r;

    } while (true);

    return R;
}
auto enum_permu_remainder(size_t n, size_t r, size_t mth) {
    std::vector<int> S{};

    std::generate_n(std::back_inserter(S), n, [count = 0]() mutable { return count++; });

    std::vector<int> R{};
    do {
        if (r == 0 || S.empty()) {
            if (!S.empty()) {
                for (auto i = 0; i < S.size(); ++i) R.push_back(S[i]);
                break;
            }
            break;
        }

        auto n = S.size();

        auto n_1_p_r_1 = jf::ncrnpr::nPr(n - 1, r - 1);

        auto quotient = mth / n_1_p_r_1;  // index of S

        R.push_back(S[quotient]);

        S.erase(S.begin() + quotient);

        mth %= n_1_p_r_1;
        --r;

    } while (true);

    return R;
}

template <template <typename, typename...> class ContainerType, typename ElementType,
          typename... Types>
int sgn(ContainerType<ElementType, Types...>& cntr) {
    // return +1 if even inversion, -1 if odd inversion
    // computers the inversion of a permutation
    ElementType count = 0, cntr_i;
    ElementType size = cntr.size();
    for (ElementType i = 1; i < size; ++i) {
        cntr_i = cntr[i];

        for (ElementType j = 0; j < i; ++j)
            if (cntr_i < cntr[j]) ++count;
    }

    return (count % 2 ? -1 : 1);
}

}  // namespace ncrnpr

}  // namespace jf

#ifdef _JF_NCRNPR_MAX_DEFINED
#pragma pop_macro("max")
#undef _JF_NCRNPR_MAX_DEFINED
#endif

#ifdef _JF_NCRNPR_MIN_DEFINED
#pragma pop_macro("min")
#undef _JF_NCRNPR_MIN_DEFINED
#endif

#endif  // _JF_NCRNPR_HPP
