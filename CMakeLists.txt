
cmake_minimum_required(VERSION 4.0.3)
##### For use of 'import std'

#VERSION 4.0.2 
#set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "a9e1cf81-9932-4810-974b-6eccaf14e457")

#VERSION 4.0.3
set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")

set(CMAKE_CXX_MODULE_STD 1)

#######
project(
    c_numeric
    VERSION 0.2.5
    DESCRIPTION "calc numeric libs"
    LANGUAGES CXX
)

# For use of 'import std' in place of '#include' directives
set(IMPORT_STD_ENABLE ON)

set(PROJECT_VERSION 4.0.2)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(MODULE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/mod)
set(THIRDY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdy)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_C_STANDARD 23)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
# Downloads CPM if not already downloaded. CPM is an easy-to-use package manager nicely integrated with CMake.
include(cmake/cpm.cmake)

# This commands downloads AND configures fmt. It sets up some variables.
# include(cmake/thirdparty.cmake)
# if (MSVC)
#     set(TBB_DIR "C:/bin/tbb-2022/oneapi-tbb-2022.0.0-win/oneapi-tbb-2022.0.0/lib/cmake/tbb")  # Caminho onde está o TBBConfig.cmake
#     #include_directories("C:/bin/TBB/include")
#     #link_directories("C:/bin/TBB/lib")
# else()
#     set(TBB_DIR "C:/msys64/uctr64/lib/cmake/TBB")  # Caminho para a versão do MSYS2
# endif()

# find_package(TBB REQUIRED)

set(TBB_ENABLE OFF)

# if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

# if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#     find_package(TBB REQUIRED )
#
#     if (TBB_FOUND)
#         message(STATUS "TBB encontrado: ${TBB_DIR}")
#         set(TBB_ENABLE ON)
#     else()
#         message(STATUS "TBB não encontrado.")
#     endif()
# endif()

# ---- Declare library ----

add_subdirectory(mod)
#add_subdirectory(source)
add_subdirectory(test)

