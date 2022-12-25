cmake_minimum_required(VERSION 3.10)
project(simple_sdr_app)

if(MSVC)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:fast /arch:AVX2")
# NOTE: If we are on a platform without AVX2
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fp:fast")
add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
else()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -mfma -ffast-math")
# NOTE: If we are compiling with ssse3 not avx2+fma (affects pll and viterbi decoder)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mssse3 -ffast-math")
endif()

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(imgui REQUIRED)
find_package(implot REQUIRED)
find_package(fmt REQUIRED)

# Windows + MSYS = Install rtlsdr from submodule
# Ubuntu = Install rtlsdr from package manager
if (WIN32) 
set(rtlsdr_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/librtlsdr)
find_package(rtlsdr CONFIG REQUIRED)
set(RTLSDR_LIBS rtlsdr)
else()
find_package(PkgConfig REQUIRED)
pkg_check_modules(rtlsdr REQUIRED IMPORTED_TARGET librtlsdr)
set(RTLSDR_LIBS PkgConfig::rtlsdr)
endif()

# MSVC = vcpkg package manager
# MSYS2 + Ubuntu = package manager
if(MSVC)
find_package(portaudio CONFIG REQUIRED)
find_package(FFTW3f CONFIG REQUIRED)
set(PORTAUDIO_LIBS portaudio)
set(FFTW3_LIBS FFTW3::fftw3f)
set(EXTRA_LIBS "")
else()
find_package(PkgConfig REQUIRED)
pkg_check_modules(portaudio REQUIRED IMPORTED_TARGET portaudio-2.0)
pkg_check_modules(fftw3f REQUIRED IMPORTED_TARGET fftw3f)
set(PORTAUDIO_LIBS PkgConfig::portaudio)
set(FFTW3_LIBS PkgConfig::fftw3f)
set(EXTRA_LIBS pthread)
endif()

# Our core libraries
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(DSP_DIR ${SRC_DIR}/dsp)
add_library(dsp_lib STATIC
    ${DSP_DIR}/filter_designer.cpp
    ${DSP_DIR}/calculate_fft.cpp
    ${DSP_DIR}/calculate_fft_mag.cpp
    ${DSP_DIR}/hilbert_fir_filter.cpp)
set_target_properties(dsp_lib PROPERTIES CXX_STANDARD 17)
target_include_directories(dsp_lib PRIVATE ${DSP_LIB} ${SRC_DIR})
target_link_libraries(dsp_lib PRIVATE ${FFTW3_LIBS})

set(DEMOD_DIR ${SRC_DIR}/demod)
add_library(demod_lib STATIC
    ${DEMOD_DIR}/fm_demod.cpp
    ${DEMOD_DIR}/broadcast_fm_demod.cpp
    ${DEMOD_DIR}/bpsk_synchroniser.cpp
    ${DEMOD_DIR}/pll_mixer.cpp
    ${DEMOD_DIR}/local_oscillator.cpp
    ${DEMOD_DIR}/trigger_cooldown.cpp
    ${DEMOD_DIR}/zero_crossing_detector.cpp
    ${DEMOD_DIR}/ted_clock.cpp)
set_target_properties(demod_lib PROPERTIES CXX_STANDARD 17)
target_include_directories(demod_lib PRIVATE ${DEMOD_LIB} ${SRC_DIR})
target_link_libraries(demod_lib PRIVATE dsp_lib)

set(RDS_DIR ${SRC_DIR}/rds_decoder)
add_library(rds_lib STATIC
    ${RDS_DIR}/rds_group_sync.cpp
    ${RDS_DIR}/rds_decoder.cpp
    ${RDS_DIR}/rds_database_decoder_handler.cpp)
set_target_properties(rds_lib PROPERTIES CXX_STANDARD 17)
target_include_directories(rds_lib PRIVATE ${RDS_DIR} ${SRC_DIR})

set(AUDIO_DIR ${SRC_DIR}/audio)
add_library(audio_lib STATIC 
    ${AUDIO_DIR}/resampled_pcm_player.cpp
    ${AUDIO_DIR}/audio_mixer.cpp
    ${AUDIO_DIR}/portaudio_output.cpp
    ${AUDIO_DIR}/portaudio_utility.cpp)
set_target_properties(audio_lib PROPERTIES CXX_STANDARD 17)
target_include_directories(audio_lib PRIVATE ${AUDIO_DIR} ${SRC_DIR})
target_link_libraries(audio_lib PRIVATE ${PORTAUDIO_LIBS})

set(GETOPT_DIR ${SRC_DIR}/utility/getopt)
add_library(getopt STATIC ${GETOPT_DIR}/getopt.c)
target_include_directories(getopt PRIVATE ${GETOPT_DIR})

set(GUI_DIR ${SRC_DIR}/gui)
add_library(gui_lib STATIC
    ${GUI_DIR}/imgui_skeleton.cpp
    ${GUI_DIR}/imgui_extensions.cpp
    ${GUI_DIR}/render_app.cpp
    ${GUI_DIR}/render_portaudio_controls.cpp
    ${GUI_DIR}/render_fm_demod.cpp
    ${GUI_DIR}/render_bpsk_sync.cpp
    ${GUI_DIR}/render_rds_database.cpp)
set_target_properties(gui_lib PROPERTIES CXX_STANDARD 17)
target_include_directories(gui_lib PRIVATE ${GUI_DIR} ${SRC_DIR})
target_link_libraries(gui_lib PRIVATE 
    imgui implot
    audio_lib demod_lib rds_lib
    getopt fmt
    ${EXTRA_LIBS})

# Our application
add_executable(main 
    ${SRC_DIR}/main.cpp 
    ${SRC_DIR}/app.cpp)
target_include_directories(main PRIVATE ${SRC_DIR})
set_target_properties(main PROPERTIES CXX_STANDARD 17)
target_link_libraries(main PRIVATE 
    imgui implot gui_lib
    ${RTLSDR_LIBS} ${FFTW3_LIBS} audio_lib dsp_lib demod_lib rds_lib
    getopt fmt
    ${EXTRA_LIBS})

if(WIN32)
    # Copy dlls for prebuilt librtlsdr
    if(CMAKE_SIZEOF_VOID_P EQUAL 8)
        add_custom_command(TARGET main 
            POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${rtlsdr_DIR}/bin/x64/
            $<TARGET_FILE_DIR:main>)
    elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
        message(FATAL_ERROR "x86 rtlsdr is unsupported")
    endif()
endif()

if(MSVC)
target_compile_options(getopt       PRIVATE "/MP")
target_compile_options(demod_lib    PRIVATE "/MP")
target_compile_options(rds_lib      PRIVATE "/MP")
target_compile_options(audio_lib    PRIVATE "/MP")
target_compile_options(main         PRIVATE "/MP")
endif()
