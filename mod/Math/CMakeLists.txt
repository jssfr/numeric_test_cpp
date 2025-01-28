add_library(math)
add_library(MATH::MATH ALIAS math)
target_sources(math
    PUBLIC
    FILE_SET CXX_MODULES FILES
    math.ixx
    math-ncrnpr.ixx
    math-numerical.ixx
    math-compose.ixx
    math-variable.ixx
    math-functions.ixx
    math-calculus.ixx
    math-fft.ixx
)

add_compile_definitions(CONTAINERS
)

target_link_libraries(math CONFIG::CONFIG)
target_include_directories(math PUBLIC ${MODULE_DIR}/headers)
