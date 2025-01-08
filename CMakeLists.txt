cmake_minimum_required(VERSION 3.31.3)

project(
    c_numeric
    VERSION 0.1.2
    DESCRIPTION "calc numeric libs"
    LANGUAGES CXX
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(IMPORT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/import)
set(THIRDY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/thirdy)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_C_STANDARD 20)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
include(cmake/project-is-top-level.cmake)
include(cmake/variables.cmake)
# Downloads CPM if not already downloaded. CPM is an easy-to-use package manager nicely integrated with CMake.
include(cmake/cpm.cmake)

# This commands downloads AND configures fmt. It sets up some variables.
include(cmake/thirdparty.cmake)

find_package(TBB REQUIRED)

# ---- Declare library ----

add_subdirectory(import)
#add_subdirectory(source)
add_subdirectory(test)


# ---- Declare executable ----

# add_executable(c_numeric_exe source/main.cpp)
# add_executable(c_numeric::exe ALIAS c_numeric_exe)
#
# set_property(TARGET c_numeric_exe PROPERTY OUTPUT_NAME c_numeric)
#
# target_compile_features(c_numeric_exe PRIVATE cxx_std_20)
#
# target_link_libraries(c_numeric_exe PRIVATE c_numeric_lib)

