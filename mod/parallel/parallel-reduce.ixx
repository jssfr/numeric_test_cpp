module;

#include <vector>
#include <future>
#include <numeric>
#include <algorithm>
#include <functional>

export module parallel:Reduce;
import :Blocked_range;

namespace jf::par
{
   
    template <typename Iterator, typename T, typename Func>
    T parallel_reduce_impl(Iterator begin, Iterator end, T init, Func&& func) {
        const unsigned num_threads = std::thread::hardware_concurrency();
        auto dist = std::distance(begin, end - 1);
        auto chunk_size = (dist + num_threads - 1) / num_threads;
        std::vector<std::future<T>> futures;

        for (unsigned i = 0; i < num_threads; ++i) {
            auto chunk_begin = std::next(begin, i * chunk_size);
            auto chunk_end = (i == num_threads - 1) ? end : std::next(chunk_begin, chunk_size);

            if (chunk_begin != chunk_end) {
                futures.push_back(std::async(std::launch::async, [=, &func]() {
                    return std::accumulate(chunk_begin, chunk_end, T{}, func);
                }));
            }
        }

        for (auto& fut : futures) {
            init = func(init, fut.get());
        }

        return init;
    }
    
    template <typename T, class Func, class Proj>
    T parallel_reduce_impl(blocked_range range, T init, Func&& func, Proj&& proj) {
        const unsigned num_threads = std::thread::hardware_concurrency();

        auto chunk_size = (range.begin() + range.end() + num_threads - 1) / num_threads;
        std::vector<std::future<T>> futures;

        for (unsigned i = 0; i < num_threads; ++i) {
            size_t chunk_begin = range.begin() + i * chunk_size;
            size_t chunk_end = std::min(chunk_begin + chunk_size, range.end()); // (i == num_threads - 1) ? range.end() : chunk_begin + chunk_size;

            if (chunk_begin != chunk_end) {
                blocked_range subrange{chunk_begin, chunk_end};
                futures.push_back(std::async(std::launch::async, [=, &func]() {
                    return func(subrange, init); // std::accumulate(chunk_begin, chunk_end, T{}, func);
                }));
            }
        }

        for (auto& fut : futures) {
            init = proj(init, fut.get());
        }

        return init;
    }

    export template <typename Iterator, typename T, typename Func>
    T parallel_reduce(Iterator begin, Iterator end, T init, Func&& func) {
        return parallel_reduce_impl(begin, end, init, func);
    }

    export template <typename T, class Func, class Proj>
    T parallel_reduce(blocked_range&& range, T init, Func&& func, Proj&& proj) {
        return parallel_reduce_impl(range, init, func, proj);
    }
}