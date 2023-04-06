cmake_minimum_required(VERSION 3.21)
project(rds_decoder)

set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
set(ROOT_DIR ${CMAKE_SOURCE_DIR}/src)

add_library(rds_decoder STATIC
    ${SRC_DIR}/crc10.cpp
    ${SRC_DIR}/rds_group_sync.cpp
    ${SRC_DIR}/rds_decoder.cpp
    ${SRC_DIR}/rds_database_decoder_handler.cpp)
set_target_properties(rds_decoder PROPERTIES CXX_STANDARD 17)
target_include_directories(rds_decoder PRIVATE ${SRC_DIR} ${ROOT_DIR})

if(MSVC)
target_compile_options(rds_decoder PRIVATE "/MP")
endif()