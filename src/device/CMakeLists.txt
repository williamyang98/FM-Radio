cmake_minimum_required(VERSION 3.21)

if(NOT DEFINED RTLSDR_LIBS)
message(FATAL_ERROR "RTLSDR_LIBS must be defined")
endif()

set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
set(ROOT_DIR ${CMAKE_SOURCE_DIR}/src)

add_library(device_lib STATIC
    ${SRC_DIR}/device.cpp
    ${SRC_DIR}/device_selector.cpp)
set_target_properties(device_lib PROPERTIES CXX_STANDARD 17)
target_include_directories(device_lib PRIVATE ${SRC_DIR} ${ROOT_DIR})
target_link_libraries(device_lib PRIVATE ${RTLSDR_LIBS} fmt)

if(MSVC)
target_compile_options(device_lib PRIVATE "/MP")
endif()