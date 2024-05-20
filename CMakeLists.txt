cmake_minimum_required(VERSION 3.4...3.18)
project(tp_3)

add_subdirectory(matplotplusplus)
add_subdirectory(pybind11)

include_directories (${AudioFile_SOURCE_DIR})

add_definitions (-DPROJECT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")
file (COPY test-audio.wav DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

pybind11_add_module(tp_3 example.cpp)

target_link_libraries(tp_3 PUBLIC matplot)
set_property(TARGET matplot PROPERTY POSITION_INDEPENDENT_CODE ON)