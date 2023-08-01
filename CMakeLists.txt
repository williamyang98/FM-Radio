cmake_minimum_required(VERSION 3.21)
project(simple_sdr_app)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(imgui REQUIRED)
find_package(implot REQUIRED)
find_package(fmt REQUIRED)

# Windows + MSYS = Install rtlsdr from submodule
# Ubuntu = Install rtlsdr from package manager
if(WIN32) 
    set(rtlsdr_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/librtlsdr)
    find_package(rtlsdr CONFIG REQUIRED)
    set(RTLSDR_LIBS rtlsdr::rtlsdr)
else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(rtlsdr REQUIRED IMPORTED_TARGET librtlsdr)
    set(RTLSDR_LIBS PkgConfig::rtlsdr)
endif()

if(WIN32)
    find_package(portaudio CONFIG REQUIRED)
    find_package(FFTW3f CONFIG REQUIRED)
    set(PORTAUDIO_LIBS portaudio)
    set(FFTW3_LIBS FFTW3::fftw3f)
else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(portaudio REQUIRED IMPORTED_TARGET portaudio-2.0)
    pkg_check_modules(fftw3f REQUIRED IMPORTED_TARGET fftw3f)
    set(PORTAUDIO_LIBS PkgConfig::portaudio)
    set(FFTW3_LIBS PkgConfig::fftw3f)
endif()

# for posix threads
if(NOT WIN32) 
    find_package(Threads REQUIRED)
    add_compile_options(-pthread) 
    link_libraries(Threads::Threads)
endif()

# simd compile options
if (NOT ${CMAKE_SYSTEM_PROCESSOR} STREQUAL "aarch64")
    if(MSVC)
        add_compile_definitions(_SILENCE_NONFLOATING_COMPLEX_DEPRECATION_WARNING)
        add_compile_definitions(_CRT_SECURE_NO_WARNINGS)
        # for MSVC multi processor compilation
        add_compile_options(/MP)
        add_compile_options(/fp:fast /arch:AVX2)
        # add_compile_options(/fp:fast /arch:AVX)
        # add_compile_options(/fp:fast)
    else()
        add_compile_options(-march=native -ffast-math)
    endif()
else()
    add_compile_options(-ffast-math)
endif()

add_subdirectory(${CMAKE_SOURCE_DIR}/src)
