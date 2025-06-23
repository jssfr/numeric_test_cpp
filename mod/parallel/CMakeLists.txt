cmake_minimum_required(VERSION ${PROJECT_VERSION})
project(parallel_modules)

add_library(parallel SHARED)
add_library(PARALLEL::PARALLEL ALIAS parallel)
target_sources(parallel
    PUBLIC
    FILE_SET CXX_MODULES FILES
    parallel.ixx
    parallel-blocked_range.ixx
    parallel-for.ixx
    parallel-for_each.ixx
    parallel-reduce.ixx
)

target_link_libraries(parallel PUBLIC CONFIG::CONFIG)

target_include_directories(parallel PUBLIC ${MODULE_DIR}/headers)

if(IMPORT_STD_ENABLE)
    add_compile_definitions(USING_IMPORT_STD_MOD)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set_target_properties(parallel PROPERTIES
        WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
endif()
