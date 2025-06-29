module;

#ifndef USING_IMPORT_STD_MOD
#include "std_input.hpp"
#endif

export module math:ncrnpr;

#ifdef USING_IMPORT_STD_MOD
import std;
#endif

namespace jf::ncrnpr {

/// @brief greater common divisor
export template <typename T> 
T gcd(T a, T b) {
  if (b > a)
    std::swap(a, b);

  T r;
  while ((r = a % b) != 0) {
    a = b;
    b = r;
  }
  return b;
}

/// @brief least common multiple
export template <typename T> 
T lcm(T a, T b) {
  T G = gcd(a, b);
  return (a / G) * b;
}

template <typename T> class Rational {
public:
  using value_t = std::make_unsigned_t<T>;

private:
  value_t m_n; // numerator
  value_t m_d; // denominator

  void reduce() {
    value_t G = gcd(this->m_n, this->m_d);

    this->m_n /= G;
    this->m_d /= G;
  }

public:
  Rational() : m_n(0), m_d(1) {}
  Rational(T n, T d) : m_n(n), m_d(d) { reduce(); }
  Rational(T n, T d, bool) : m_n(n), m_d(d) {}

  Rational(Rational &&rhs) = default;
  Rational(const Rational &rhs) = default;

  auto operator=(const Rational &rhs) -> Rational & = default;

  auto operator=(Rational &&rhs) -> Rational & = default;

  ~Rational() = default;

  template <typename C> operator C() const {
    return C((double)this->m_n / (double)this->m_d);
  }

  auto operator*=(const Rational &rat) -> Rational & {
    Rational rat1(this->m_n, rat.m_d);
    Rational rat2(rat.m_n, this->m_d);

    this->m_n = rat1.m_n * rat2.m_n;
    this->m_d = rat1.m_d * rat2.m_d;

    return *this;
  }

  auto operator+=(const Rational &rat) -> Rational & {
    T num = this->m_n, den = this->m_d;
    T Least = lcm(m_d, rat.m_d);

    this->m_n = num * Least / den + rat.m_n * (Least / rat.m_d);
    this->m_d = Least;

    reduce();
    return *this;
  }
};

export template <typename T> T fact(T n) {
  T factorial{n};
  if (n == 0)
    return T{1};
  for (auto i = n - 1; i > 0; --i)
    factorial *= i;
  return factorial;
}

/// @return comutation of n r
export template <typename T> 
T nCr(T n, T r) {
  if (r > n)
    return 0;
  if (r == 0 || n == r)
    return 1;
  if (r == (n - 1) || r == 1)
    return n;
  
  if (r > (n - r))
     r = n - r;

  Rational<T> rlt{n, r};

  for (T i = 1; i < r; ++i)
     rlt *= Rational<T>(n - i, r - i);

  return (T)rlt;
  
}

/// @return permutation n r
export template <typename T> 
T nPr(T n, T r) {
  if (r > n)
    return 0;
  if (r == 0)
    return 1;
  
  T count{n};
  do {
    if (r == 1)
      break;
    if (r == n || r == (n - 1)) {
        count = fact(n);
        break;

    } 
    --n;
    --r;
    count *= n;
    
  } while (true);

  return count;
  
}

// export template <typename T> 
// T combi_double(T n, T r) {
//   if (r > n)
//     return 0;
//   else if (r == 0 || n == r)
//     return 1;
//   else if (r == (n - 1) || r == 1)
//     return n;
//   else {
//     if (r > (n - r))
//       r = n - r;
//
//     double rlt{n, r};
//
//     for (T i = 1; i < r; ++i)
//       rlt *= double(n - i) / double(r - i);
//
//     return T(rlt);
//   }
// }

// enum combination

template <class set_of_sets_t = std::vector<std::vector<int>>>
void enum_comb(set_of_sets_t &&SS, std::size_t n, std::size_t r, std::size_t mth) {
  std::vector<int> S(n);

  std::generate(S.begin(), S.end(),
                [n = std::size_t{}]() mutable { return ++n; });
  std::vector<int> R;
  do {
    if (r == 0 && !R.empty()) {
      SS.emplace_back(R);
      break;
    }  
    if (r == S.size()) {
      R.insert(R.end(), std::make_move_iterator(S.begin()),
               std::make_move_iterator(S.end()));
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

export template <class set_t = std::vector<int>>
void enum_permu(set_t &R, set_t S, std::size_t r, std::size_t mth) {
  // https://youtu.be/gixbV9dMlRg?si=yXEF7Jw8IICofgAg
  do {
    if (r == 0 || S.empty()) {
      break;
    }

    auto n = S.size();

    auto n_1_p_r_1 = jf::ncrnpr::nPr(n - 1, r - 1);

    auto quotient = mth / n_1_p_r_1; // index of S

    R.push_back(S[quotient]);

    S.erase(S.begin() + quotient);

    mth %= n_1_p_r_1;
    --r;

  } while (true);
}
export auto enum_permu(std::size_t n, std::size_t r, std::size_t mth) {
  std::vector<int> S{};

  std::generate_n(std::back_inserter(S), n,
                  [count = 0]() mutable { return count++; });

  std::vector<int> R{};
  do {
    if (r == 0 || S.empty()) {
      break;
    }

    auto n = S.size();

    auto n_1_p_r_1 = jf::ncrnpr::nPr(n - 1, r - 1);

    auto quotient = mth / n_1_p_r_1; // index of S

    R.push_back(S[quotient]);

    S.erase(S.begin() + quotient);

    mth %= n_1_p_r_1;
    --r;

  } while (true);

  return R;
}
export auto enum_permu_remainder(std::size_t n, std::size_t r, std::size_t mth) {
  std::vector<int> S{};

  std::generate_n(std::back_inserter(S), n,
                  [count = 0]() mutable { return count++; });

  std::vector<int> R{};
  do {
    if (r == 0 || S.empty()) {
      if (!S.empty()) {
        for (std::size_t i{}; i < S.size(); ++i)
          R.push_back(S[i]);
        break;
      }
      break;
    }

    auto n = S.size();

    auto n_1_p_r_1 = jf::ncrnpr::nPr(n - 1, r - 1);

    auto quotient = mth / n_1_p_r_1; // index of S

    R.push_back(S[quotient]);

    S.erase(S.begin() + quotient);

    mth %= n_1_p_r_1;
    --r;

  } while (true);

  return R;
}

export template <template <typename, typename...> class ContainerType, typename ElementType,
       typename... Types>
/// @return +1 if even inversion, -1 if odd inversion
/// @brief computers the inversion of a permutation
int sgn(ContainerType<ElementType, Types...> &cntr) {
  ElementType count{}, cntr_i{};
  ElementType size = cntr.size();
  for (auto i = ElementType{1}; i < size; ++i) {
    cntr_i = cntr[i];

    for (ElementType j{}; j < i; ++j)
      if (cntr_i < cntr[j])
        ++count;
  }

  return (count % 2 ? -1 : 1);
}

} // namespace jf::ncrnpr
