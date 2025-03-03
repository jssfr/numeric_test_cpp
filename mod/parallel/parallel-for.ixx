module;

#include <algorithm>
#include <execution>
#include <thread>
#include <functional>
#include <vector>

export module parallel:For;
import :Blocked_range;

namespace jf::par{
    
    template<typename Index, class Func>
    void parallel_for_impl(Index start, Index end, Index step, Func&& func)
    {
        const unsigned num_threads = std::thread::hardware_concurrency();

        std::vector<std::jthread> threads;

        Index chunk_size = (end - start + num_threads - 1) / num_threads;

        for (unsigned i = 0; i < num_threads; ++i) {
            Index chunk_start = start + i * chunk_size; //i=0: 0 ; 
            Index chunk_end = std::min(chunk_start + chunk_size, end);// i=0: 1;

            if (chunk_start < chunk_end) {
                threads.emplace_back([=, &func]() {
                    for (Index j = chunk_start; j < chunk_end; j += step) {
                        func(j);
                    }
                });
            }

        }
    }

    template<class Func>
    void parallel_for_impl(blocked_range range, Func&& func)
    {
        const unsigned num_threads = std::thread::hardware_concurrency();

        std::vector<std::jthread> threads;

        size_t chunk_size = (range.end() - range.begin() + num_threads - 1) / num_threads;

        for (unsigned i = 0; i < num_threads; ++i) {
            size_t chunk_start = range.begin() + i * chunk_size;
            size_t chunk_end = std::min(chunk_start + chunk_size, range.end());

            if (chunk_start < chunk_end) {

                blocked_range subrange{chunk_start, chunk_end};

                threads.emplace_back([=, &func]() {
                           func(subrange);
                });
            }

        }
    }

    export template<typename Index, class Func>
    void parallel_for(Index start, Index end, Func&& func)
    {
        parallel_for_impl(start, end, Index{1}, func);
    }
    
    export template<typename Index, class Func>
    void parallel_for(Index start, Index end, Index step, Func&& func)
    {
        parallel_for_impl(start, end, step, func);
    }
    
    export template<class Func>
    void parallel_for(blocked_range&& range, Func&& func)
    {
        parallel_for_impl(range, func);
    }
    
}