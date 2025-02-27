cmake_minimum_required(VERSION 3.31.3)

project(
    c_numeric
    VERSION 0.1.2
    DESCRIPTION "calc numeric libs"
    LANGUAGES CXX
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(MODULE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/mod)
set(THIRDY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdy)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_C_STANDARD 20)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
# Downloads CPM if not already downloaded. CPM is an easy-to-use package manager nicely integrated with CMake.
include(cmake/cpm.cmake)

find_package(TBB REQUIRED)
if (TBB_FOUND)
        message(STATUS "TBB encontrado: ${TBB_DIR}")
        set(TBB_ENABLE ON)
    else()
        message(STATUS "TBB não encontrado.")
        set(TBB_ENABLE OFF)
    endif()

# ---- Declare library ----

add_subdirectory(mod)
#add_subdirectory(source)
add_subdirectory(test)

