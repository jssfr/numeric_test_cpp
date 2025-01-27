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
        import
        fmt
        stdc++exp
)

target_include_directories(
    ${EXE}
     PUBLIC
     ${IMPORT_DIR}
 )

add_compile_definitions(CONTAINERS
                    USING_TBBLIB
)

