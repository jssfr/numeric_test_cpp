install(
    TARGETS c_numeric_exe
    RUNTIME COMPONENT c_numeric_Runtime
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
