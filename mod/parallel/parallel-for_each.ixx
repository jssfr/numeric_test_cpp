module;

#ifndef USING_IMPORT_STD_MOD
  #include <vector>
  #include <thread>
  #include <iterator>
#endif

export module parallel:For_each;
#ifdef USING_IMPORT_STD_MOD
  import std;
#endif

namespace jf::par
{
    template <typename Iterator, typename Func>
    void parallel_for_each_impl(Iterator begin, Iterator end, Func&& func) {
        const unsigned num_threads = std::thread::hardware_concurrency();
        std::vector<std::jthread> threads;
        auto dist = std::distance(begin, end - 1);
        auto chunk_size = (dist + num_threads - 1) / num_threads;

        for (unsigned i = 0; i < num_threads; ++i) {
            auto chunk_begin = std::next(begin, i * chunk_size);
            auto chunk_end = (i == num_threads - 1) ? end : std::next(chunk_begin, chunk_size);

            if (chunk_begin != chunk_end) {
                threads.emplace_back([=, &func]() {
                    for (auto it = chunk_begin; it != chunk_end; ++it) {
                        func(*it);
                    }
                });
            }
        }
    }
    export template <typename Iterator, typename Func>
    void parallel_for_each(Iterator begin, Iterator end, Func&& func) {
        parallel_for_each_impl(begin, end, func);
    }
}
