# CMakeLists.txt

# add_library(config)
# add_library(CONFIG::CONFIG ALIAS config)
#
# target_sources(config
#     PUBLIC
#     FILE_SET all_my_modules TYPE CXX_MODULES FILES
#     ./Config/Config.ixx
#     ./Config/Config-Name.ixx
#     ./Config/Config-Types.ixx
# )
# target_include_directories(config PUBLIC ${IMPORT_DIR}/headers)
#
#
# add_library(math)
# add_library(MATH::MATH ALIAS math)
# target_sources(math
#     PUBLIC
#     FILE_SET all_my_modules TYPE CXX_MODULES FILES
#     ./Math/math.ixx
#     ./Math/math-ncrnpr.ixx
# )
#
#
# target_link_libraries(math TBB::tbb)
# target_link_libraries(math TBB::tbbmalloc) # for allocation
# target_include_directories(math PUBLIC ${IMPORT_DIR}/headers)
#
#
# add_library(matrix)
# add_library(MATRIX::MATRIX ALIAS matrix)
# target_sources(matrix
#     PUBLIC
#     FILE_SET all_my_modules TYPE CXX_MODULES FILES
#     ./Matrix/matrix.ixx
#     ./Matrix/matrix2d.ixx
# )
#
# add_compile_definitions(CONTAINERS
# 		    USING_FMTLIB
#                     USING_TBBLIB
# )
#
# target_link_libraries(matrix fmt CONFIG::CONFIG MATH::MATH)
# target_link_libraries(matrix TBB::tbb)
# target_link_libraries(matrix TBB::tbbmalloc) # for allocation
# target_include_directories(matrix PUBLIC ${IMPORT_DIR}/headers)
#
add_subdirectory(Config)
add_subdirectory(Math)
add_subdirectory(Matrix)

add_library(import INTERFACE)
target_link_libraries(import
	INTERFACE
	CONFIG::CONFIG
	MATH::MATH
	MATRIX::MATRIX
)
