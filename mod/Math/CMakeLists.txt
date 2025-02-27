add_library(math SHARED)
add_library(MATH::MATH ALIAS math)
target_sources(math
    PUBLIC
    FILE_SET CXX_MODULES FILES
    math.ixx
    math-ncrnpr.ixx
    math-numerical.ixx
    math-compose.ixx
    math-variable.ixx
    math-functions.ixx
    math-calculus.ixx
    math-fft.ixx
    matrix2d.ixx
    math-base.ixx
)

add_compile_definitions(CONTAINERS
		    #USING_FMTLIB
                    )

target_link_libraries(math PUBLIC CONFIG::CONFIG)
target_link_libraries(math PUBLIC PARALLEL::PARALLEL)

if(TBB_ENABLE)
    add_compile_definitions(USING_TBBLIB)
    target_link_libraries(math PUBLIC TBB::tbb)
    target_link_libraries(math PUBLIC TBB::tbbmalloc) # for allocation
endif()

target_include_directories(math PUBLIC ${MODULE_DIR}/headers)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set_target_properties(math PROPERTIES
        WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
endif()
