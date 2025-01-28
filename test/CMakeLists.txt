## ---- Tests ----
# set(SOURCE main.cpp)
# set(EXE test01)
# set(SOURCE fft_test.cpp)
# set(EXE test02)

set(SOURCE lambda_impl.cpp)
set(EXE lambda_impl)

add_executable(${EXE} ${SOURCE})
target_link_libraries(${EXE}
		PRIVATE
        mod
        fmt
        stdc++exp
)

target_include_directories(
    ${EXE}
     PUBLIC
     ${MODULE_DIR}
 )

add_compile_definitions(CONTAINERS
                    USING_TBBLIB
)

