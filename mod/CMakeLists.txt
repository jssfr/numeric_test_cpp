# CMakeLists.txt
add_subdirectory(Config)
add_subdirectory(Math)
#add_subdirectory(Matrix)

add_subdirectory(parallel)

add_library(mod INTERFACE)
target_link_libraries(mod
	INTERFACE
	CONFIG::CONFIG
	MATH::MATH
# 	 MATRIX::MATRIX
)

    target_link_libraries(mod
        INTERFACE
        PARALLEL::PARALLEL
    )

