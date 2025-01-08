## ---- Tests ----
set(SOURCE main.cpp)
set(EXE test01)

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

