cmake_minimum_required(VERSION 3.21)

if(NOT DEFINED PORTAUDIO_LIBS)
    message(FATAL_ERROR "PORTAUDIO_LIBS must be defined")
endif()

set(SRC_DIR ${CMAKE_CURRENT_LIST_DIR})
set(ROOT_DIR ${CMAKE_SOURCE_DIR}/src)

add_library(audio_lib STATIC 
    ${SRC_DIR}/resampled_pcm_player.cpp
    ${SRC_DIR}/audio_mixer.cpp
    ${SRC_DIR}/portaudio_output.cpp
    ${SRC_DIR}/portaudio_utility.cpp)
set_target_properties(audio_lib PROPERTIES CXX_STANDARD 17)
target_include_directories(audio_lib PRIVATE ${SRC_DIR} ${ROOT_DIR})
target_link_libraries(audio_lib PRIVATE ${PORTAUDIO_LIBS})

if(MSVC)
target_compile_options(audio_lib PRIVATE "/MP")
endif()