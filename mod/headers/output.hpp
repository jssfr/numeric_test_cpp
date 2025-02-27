

#ifndef _JF_OUTPUT_HPP
#define _JF_OUTPUT_HPP

#ifdef _MSVC_LANG

    #if _MSVC_LANG < 201703L
        #error This library requires c++17 or later
    #endif
	#define WIN32_LEAN_AND_MEAN      // Exclude rarely-used stuff from Windows headers
	#include <windows.h>

#else
    #if __cplusplus < 201703L
        #error This library requires c++17 or later
    #endif
#endif // _MSVC_LANG

/*
    INCLUDES HERE
*/
#ifdef max
    #undef max
#endif

#ifdef min
    #undef min
#endif
#include<iostream>

#ifdef USING_TBBLIB
#include <oneapi/tbb.h>
#include<oneapi/tbb/partitioner.h>
#include <oneapi/tbb/concurrent_vector.h>
#include <oneapi/tbb/tbb_allocator.h>
#include <oneapi/tbb/scalable_allocator.h>


#endif
#ifndef max
    #define max
#endif

#ifndef min
    #define min
#endif

#ifdef USING_FMTLIB
#include<fmt/core.h>
#include<fmt/color.h>
#include<fmt/format.h>
#include<fmt/ranges.h>
#endif


#endif // _JF_OUTPUT_HPP
