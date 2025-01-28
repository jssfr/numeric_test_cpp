add_library(config)
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


