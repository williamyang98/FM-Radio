cmake_minimum_required(VERSION 3.21)

set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
set(ROOT_DIR ${CMAKE_SOURCE_DIR}/src)

add_library(fm_demod STATIC
    ${SRC_DIR}/fm_demod.cpp
    ${SRC_DIR}/broadcast_fm_demod.cpp
    ${SRC_DIR}/bpsk_synchroniser.cpp
    ${SRC_DIR}/pll_mixer.cpp
    ${SRC_DIR}/local_oscillator.cpp
    ${SRC_DIR}/trigger_cooldown.cpp
    ${SRC_DIR}/zero_crossing_detector.cpp
    ${SRC_DIR}/ted_clock.cpp)
set_target_properties(fm_demod PROPERTIES CXX_STANDARD 17)
target_include_directories(fm_demod PRIVATE ${SRC_DIR} ${ROOT_DIR})
target_link_libraries(fm_demod PRIVATE dsp_lib)

if(MSVC)
target_compile_options(fm_demod PRIVATE "/MP")
endif()