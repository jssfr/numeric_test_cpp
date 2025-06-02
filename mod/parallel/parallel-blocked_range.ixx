module;
#ifndef USING_IMPORT_STD_MOD
  #include<algorithm>
#endif

export module parallel:Blocked_range;
#ifdef USING_IMPORT_STD_MOD
  import std;
#endif

namespace jf::par   
{
    export class blocked_range
    {
        std::size_t m_begin, m_end;

        public:

        blocked_range(std::size_t begin, std::size_t end): m_begin{begin}, m_end{end} {}
        blocked_range(blocked_range&) = default;
        blocked_range(const blocked_range&) = default;
        blocked_range& operator=(const blocked_range&) = default;

        std::size_t begin()const { return m_begin;}
        std::size_t end()const { return m_end;}
    };
} // namespace jf::par
