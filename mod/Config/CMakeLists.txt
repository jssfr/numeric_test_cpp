add_library(config SHARED)
add_library(CONFIG::CONFIG ALIAS config)

target_sources(config
    PUBLIC
    FILE_SET CXX_MODULES FILES
    Config.ixx
    Config-Name.ixx
    Config-Types.ixx
)
add_compile_definitions(CONTAINERS
)
target_include_directories(config PUBLIC ${MODULE_DIR}/headers)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set_target_properties(config PROPERTIES
        WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
endif()
