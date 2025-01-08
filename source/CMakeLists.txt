# object from classes

add_library(
    windowry_lib OBJECT
)

#target_link_libraries(windowry_lib
#    PUBLIC
#    planes
#)
target_sources(windowry_lib
    PUBLIC
    FILE_SET all_my_modules TYPE CXX_MODULES FILES
   ./../import/BasicPlane.Figures-Point.ixx
    ./../import/BasicPlane.Figures-Rectangle.ixx
    BasicPlane.Figures-Rectangle.cpp
)

# add_executable(SDLtest01)
# target_sources(SDLtest01 PRIVATE main.cpp )
#
 # target_include_directories(
 #     windowry_lib ${warning_guard}
 #     PUBLIC
 #     "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/import>"
 # )
#
# target_link_libraries(SDLtest01
# 		PRIVATE
#  		vendor
#         $<TARGET_OBJECTS:sdl_windowry_lib>
# )


