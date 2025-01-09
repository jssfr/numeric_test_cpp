# CMakeLists.txt
add_subdirectory(Config)
add_subdirectory(Math)
# add_subdirectory(Matrix)

add_library(import INTERFACE)
target_link_libraries(import
	INTERFACE
	CONFIG::CONFIG
	MATH::MATH
    #	MATRIX::MATRIX
)
