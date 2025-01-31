## ---- Tests ----
# set(SOURCE main.cpp)
# set(EXE test01)
# set(SOURCE fft_test.cpp)
# set(EXE test02)

# set(SOURCE lambda_impl.cpp)
# set(EXE lambda_impl)
set(SOURCE numeric2.cpp)
set(EXE numeric2)

add_executable(${EXE} ${SOURCE})
target_link_libraries(${EXE}
		PRIVATE
        mod
        fmt
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # for use of std::print() 
    target_link_libraries(${EXE}
        PRIVATE
        stdc++exp)
endif()


target_include_directories(
    ${EXE}
     PUBLIC
     ${MODULE_DIR}
 )

add_compile_definitions(CONTAINERS
                    USING_TBBLIB
)

