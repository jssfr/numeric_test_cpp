add_library(parallel SHARED)
add_library(PARALLEL::PARALLEL ALIAS parallel)
target_sources(parallel
    PUBLIC
    FILE_SET CXX_MODULES FILES
    parallel.ixx
    parallel-blocked_range.ixx;
    parallel-for.ixx
    parallel-for_each.ixx
    parallel-reduce.ixx
)

#add_compile_definitions(CONTAINERS )

#target_link_libraries(parallel PUBLIC CONFIG::CONFIG)

#target_include_directories(parallel PUBLIC ${MODULE_DIR}/headers)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set_target_properties(parallel PROPERTIES
        WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
endif()
