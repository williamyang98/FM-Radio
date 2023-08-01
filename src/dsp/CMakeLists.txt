cmake_minimum_required(VERSION 3.21)

set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
set(ROOT_DIR ${CMAKE_SOURCE_DIR}/src)

if(NOT DEFINED FFTW3_LIBS)
    message(FATAL_ERROR "FFTW3_LIBS must be defined")
endif()

add_library(dsp_lib STATIC
    ${SRC_DIR}/filter_designer.cpp
    ${SRC_DIR}/calculate_fft.cpp
    ${SRC_DIR}/calculate_fft_mag.cpp
    ${SRC_DIR}/simd/apply_harmonic_pll.cpp
    ${SRC_DIR}/simd/c32_f32_cum_mul.cpp
    ${SRC_DIR}/simd/f32_cum_mul.cpp
)
set_target_properties(dsp_lib PROPERTIES CXX_STANDARD 17)
target_include_directories(dsp_lib PRIVATE ${SRC_DIR} ${ROOT_DIR})
target_link_libraries(dsp_lib PRIVATE ${FFTW3_LIBS})