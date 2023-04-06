cmake_minimum_required(VERSION 3.21)

if(NOT DEFINED RTLSDR_LIBS)
    message(FATAL_ERROR "RTLSDR_LIBS must be defined")
endif()

set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
set(ROOT_DIR ${CMAKE_SOURCE_DIR}/src)

add_library(app_gui STATIC
    ${SRC_DIR}/imgui_skeleton.cpp
    ${SRC_DIR}/imgui_extensions.cpp
    ${SRC_DIR}/render_app.cpp
    ${SRC_DIR}/render_portaudio_controls.cpp
    ${SRC_DIR}/render_fm_demod.cpp
    ${SRC_DIR}/render_bpsk_sync.cpp
    ${SRC_DIR}/render_rds_database.cpp)
set_target_properties(app_gui PROPERTIES CXX_STANDARD 17)
target_include_directories(app_gui PRIVATE ${SRC_DIR} ${ROOT_DIR})
target_link_libraries(app_gui PRIVATE 
    imgui implot
    audio_lib fm_demod rds_decoder
    getopt fmt)

add_library(device_gui STATIC
    ${SRC_DIR}/render_device_selector.cpp
    ${SRC_DIR}/imgui_extensions.cpp)
set_target_properties(device_gui PROPERTIES CXX_STANDARD 17)
target_include_directories(device_gui PRIVATE ${SRC_DIR} ${ROOT_DIR})
target_link_libraries(device_gui PRIVATE ${RTLSDR_LIBS} device_lib imgui fmt)

if(MSVC)
target_compile_options(device_gui PRIVATE "/MP")
endif()