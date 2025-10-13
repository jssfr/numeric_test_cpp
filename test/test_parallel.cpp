
#ifndef USING_IMPORT_STD_MOD
#  include <chrono>
#  include <execution>
#  include <functional>
#  include <numeric>
#  include <print>
#  include <thread>
#  include <utility>
#  include <vector>
#endif

#define SEQ std::execution::seq
#define PAR std::execution::par

#ifdef USING_IMPORT_STD_MOD
import std;
#endif
import parallel;

void test_parallel_for()
{
  jf::par::parallel_for(0, 10,
      [](auto i) { std::print("Processando {} na thread {}\n", i, std::this_thread::get_id()); });

  std::println("");

  jf::par::parallel_for(jf::par::blocked_range { std::size_t {}, static_cast<std::size_t>(10) },
      [](auto& sub_range)
      {
        for (auto i = sub_range.begin(); i < sub_range.end(); ++i) {
            std::print("Na sub range. Processando {} na thread {}\n", i,
                std::this_thread::get_id());
          }
      });
}

void test_parallel_for_each()
{
  std::vector<int> data = { 1, 2, 3, 4, 5 };
  jf::par::parallel_for_each(data.begin(), data.end(),
      [](int& x)
      {
        x *= 2;  // Multiplica cada elemento por 2
      });

  std::print("value of data = {}\n", data);
}

void test_parallel_reduce()
{
  std::vector<int> numbers = { 1, 2, 3, 4, 5, 6, 7, 8 };
  int sum1 = jf::par::parallel_reduce(numbers.begin(), numbers.end(), int {},
      [](auto i, auto a) { return i + a; });

  auto handle = [&](auto& range, auto det)
  {
    decltype(det) sum {};
    for (auto j = range.begin(); j != range.end(); ++j) {
        sum += j;
      }
    return sum;
  };

  auto sum_up = [](auto left_det, auto right_det) { return left_det + right_det; };

  int sum2 = jf::par::parallel_reduce(jf::par::blocked_range { std::size_t {},
                                          static_cast<std::size_t>(9) },
      std::size_t {}, handle, sum_up);
  std::print("Soma: {} == {}\n", sum1, sum2);  // Deve imprimir 36

  auto eval = [](auto fun)
  {
    const auto t1 = std::chrono::high_resolution_clock::now();
    const auto [name, result] = fun();
    const auto t2 = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double, std::milli> ms = t2 - t1;
    std::print("{}, sum: {}     time: {} ms\n", name, result, ms.count());
  };

  {
    const std::vector<double> v(100'000'007, 0.1);

    eval(
        [&v]
        {
          return std::pair { "jf::par::parallel_reduce (double)",
            jf::par::parallel_reduce(v.cbegin(), v.cend(), 0.0,
                [](auto i, auto a) { return i + a; }) };
        });  //[](auto i, auto a){ return  i + a;}
    eval(
        [&v]
        {
          return std::pair { "std::reduce (seq, double)", std::reduce(SEQ, v.cbegin(), v.cend()) };
        });
    eval(
        [&v]
        {
          return std::pair { "std::reduce (par, double)", std::reduce(PAR, v.cbegin(), v.cend()) };
        });
  }

  {
    const std::vector<long> v(100'000'007, 1);

    eval(
        [&v]
        {
          return std::pair { "jf::par::parallel_reduce (long)",
            jf::par::parallel_reduce(v.cbegin(), v.cend(), 0l,
                [](auto i, auto a) { return i + a; }) };
        });
    eval(
        [&v]
        {
          return std::pair { "std::reduce (seq, long)", std::reduce(SEQ, v.cbegin(), v.cend()) };
        });
    eval(
        [&v]
        {
          return std::pair { "std::reduce (par, long)", std::reduce(PAR, v.cbegin(), v.cend()) };
        });
  }
}

int main()
{
  // test_parallel_for();

  // test_parallel_for_each();

  test_parallel_reduce();
  return 0;
}
