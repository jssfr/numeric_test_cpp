module;
#include<algorithm>

export module parallel:Blocked_range;

namespace jf::par   
{
    export class blocked_range
    {
        size_t m_begin, m_end;

        public:

        blocked_range(size_t begin, size_t end): m_begin{begin}, m_end{end} {}
        blocked_range(blocked_range&) = default;
        blocked_range(const blocked_range&) = default;
        blocked_range& operator=(const blocked_range&) = default;

        size_t begin()const { return m_begin;}
        size_t end()const { return m_end;}
    };
} // namespace jf::par
