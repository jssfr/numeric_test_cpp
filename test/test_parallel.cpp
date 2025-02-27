#include<print>
#include <vector>
#include<thread>

import parallel;

void test_parallel_for(){
    jf::par::parallel_for(0, 10, [](auto i) {
        std::print("Processando {} na thread {}\n", i, std::this_thread::get_id());
    });

    std::println("");

    jf::par::parallel_for(jf::par::blocked_range{size_t{}, static_cast<size_t>(10)}, [](auto& sub_range) {
        for (int i = sub_range.begin(); i < sub_range.end(); ++i) {
            std::print("Na sub range. Processando {} na thread {}\n", i, std::this_thread::get_id());
        }
    });
    
}

void test_parallel_for_each(){
    std::vector<int> data = {1, 2, 3, 4, 5};
    jf::par::parallel_for_each(data.begin(), data.end(), [](int& x) {
    x *= 2; // Multiplica cada elemento por 2
    });

    std::print("value of data = {}\n", data);

}

void test_parallel_reduce()
{
    std::vector<int> numbers = {1, 2, 3, 4, 5, 6, 7, 8};
    int sum1 = jf::par::parallel_reduce(numbers.begin(), numbers.end(), int{}, [](auto i, auto a){ return  i + a;});
        
    auto handle = [&](auto& range, auto det) {
        decltype(det) sum{};
        for (auto j = range.begin(); j != range.end(); ++j) {
            sum += j;
        }
        return sum;
    };

    auto sum_up = [](auto left_det, auto right_det) { return left_det + right_det; };

    int sum2 = jf::par::parallel_reduce(jf::par::blocked_range{size_t{}, static_cast<size_t>(9)}, size_t{}, handle, sum_up);
    std::print("Soma: {} == {}\n", sum1, sum2 ); // Deve imprimir 36


}

int main()
{

    // test_parallel_for();

    // test_parallel_for_each();

    test_parallel_reduce();
    return 0;
}