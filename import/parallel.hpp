#ifndef PARALLEL_HPP_
#define PARALLEL_HPP_

#include "types_value.hpp"
#include <atomic>
#include <algorithm>
#include <execution>


/*
     *@brief implements set operations
*/
namespace jf::parallel
{

template<typename ElementType, typename Operation>
inline void atomic_operation(std::atomic<ElementType>& atom, Operation&& opr)
{
    ElementType old_value = atom;
    ElementType new_value = opr(old_value);

    while (!atom.compare_exchange_strong(old_value, new_value))
    {
        old_value = atom;
        new_value = opr(old_value);
    }

}

}// end of namespace jf::parallel 


#endif
