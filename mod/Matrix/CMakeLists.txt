add_library(matrix SHARED)
add_library(MATRIX::MATRIX ALIAS matrix)
target_sources(matrix
    PUBLIC
    FILE_SET all_my_modules TYPE CXX_MODULES FILES
    matrix.ixx
    matrix2d.ixx
)

add_compile_definitions(CONTAINERS
		    # USING_FMTLIB
                    USING_TBBLIB
)

target_link_libraries(matrix MATH::MATH)
target_link_libraries(matrix TBB::tbb)
target_link_libraries(matrix TBB::tbbmalloc) # for allocation
target_include_directories(matrix PUBLIC ${IMPORT_DIR}/headers)
