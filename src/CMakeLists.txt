cmake_minimum_required(VERSION 3.21)

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/dsp)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/fm_demod)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/rds_decoder)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/audio)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/device)
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/gui)

if (NOT DEFINED RTLSDR_LIBS)
    message(FATAL_ERROR "RTLSDR_LIBS must be defined")
endif()

set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR})

set(GETOPT_DIR ${SRC_DIR}/getopt)
add_library(getopt STATIC ${GETOPT_DIR}/getopt.c)
target_include_directories(getopt PRIVATE ${GETOPT_DIR})

function(init_target target)
    target_compile_features(${target} PRIVATE cxx_std_17)
    if(MSVC)
    target_compile_options(${target} PRIVATE "/MP")
    else()
    target_link_libraries(${target} PRIVATE pthread)
    endif()
endfunction()

function(install_dlls target)
    if (WIN32)
        add_custom_command(
            TARGET ${target}
            POST_BUILD 
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            $<TARGET_RUNTIME_DLLS:${target}> 
            $<TARGET_FILE_DIR:${target}>
            COMMAND_EXPAND_LISTS
        )
    endif()
endfunction()

add_library(app_lib STATIC ${SRC_DIR}/app.cpp)
target_include_directories(app_lib PRIVATE ${SRC_DIR})
target_link_libraries(app_lib PRIVATE dsp_lib fm_demod rds_decoder)
init_target(app_lib)

add_library(fm_scraper STATIC ${SRC_DIR}/fm_scraper.cpp)
target_include_directories(fm_scraper PRIVATE ${SRC_DIR})
init_target(fm_scraper)

# Utility applications
add_executable(rtl_sdr ${SRC_DIR}/rtl_sdr.cpp)
target_include_directories(rtl_sdr PRIVATE ${SRC_DIR})
init_target(rtl_sdr)
target_link_libraries(rtl_sdr PRIVATE ${RTLSDR_LIBS} getopt)
install_dlls(rtl_sdr)

add_executable(fm_demod_benchmark ${SRC_DIR}/fm_demod_benchmark.cpp)
target_include_directories(fm_demod_benchmark PRIVATE ${SRC_DIR})
target_link_libraries(fm_demod_benchmark PRIVATE 
    getopt app_lib)
init_target(fm_demod_benchmark)

add_executable(fm_demod_scraper ${SRC_DIR}/fm_demod_scraper.cpp)
target_include_directories(fm_demod_scraper PRIVATE ${SRC_DIR})
target_link_libraries(fm_demod_scraper PRIVATE 
    getopt dsp_lib fm_demod fm_scraper)
init_target(fm_demod_scraper)

add_executable(rds_decode ${SRC_DIR}/rds_decode.cpp)
target_include_directories(rds_decode PRIVATE ${SRC_DIR})
target_link_libraries(rds_decode PRIVATE 
    getopt rds_decoder)
init_target(rds_decode)

# GUI applications
add_executable(fm_demod_no_tuner ${SRC_DIR}/fm_demod_no_tuner.cpp)
target_include_directories(fm_demod_no_tuner PRIVATE ${SRC_DIR})
target_link_libraries(fm_demod_no_tuner PRIVATE 
    getopt
    imgui implot
    app_lib app_gui
    audio_lib ${PORTAUDIO_LIBS})
init_target(fm_demod_no_tuner)
install_dlls(fm_demod_no_tuner)

add_executable(fm_demod_tuner ${SRC_DIR}/fm_demod_tuner.cpp)
target_include_directories(fm_demod_tuner PRIVATE ${SRC_DIR})
target_link_libraries(fm_demod_tuner PRIVATE 
    getopt
    imgui implot
    app_lib app_gui
    audio_lib ${PORTAUDIO_LIBS}
    device_lib device_gui ${RTLSDR_LIBS})
init_target(fm_demod_tuner)
install_dlls(fm_demod_tuner)
