## ---- Tests ----
# set(SOURCE main.cpp)
# set(EXE test01)
set(SOURCE fft_test.cpp)
set(EXE test02)

# set(SOURCE lambda_impl.cpp)
# set(EXE lambda_impl)
# set(SOURCE numeric2.cpp)
# set(EXE numeric2)
# set(SOURCE like-t.cpp)
# set(EXE like_t)
# set(SOURCE testmatrix.cpp)
# set(EXE mat)
# set(SOURCE test_parallel.cpp)
# set(EXE par1)
# set(SOURCE test_householder.cpp)
# set(EXE house)
add_executable(${EXE} ${SOURCE})
target_link_libraries(${EXE}
	PRIVATE
        mod
       # fmt
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
    # for use of std::print() 
    target_link_libraries(${EXE}
        PRIVATE
	stdc++exp
        )
endif()


target_include_directories(
    ${EXE}
     PUBLIC
     ${MODULE_DIR}
 )

add_compile_definitions(CONTAINERS)

if(IMPORT_STD_ENABLE)
    add_compile_definitions(USING_IMPORT_STD_MOD)
endif()


if(TBB_ENABLE)
    add_compile_definitions(USING_TBBLIB)
endif()
