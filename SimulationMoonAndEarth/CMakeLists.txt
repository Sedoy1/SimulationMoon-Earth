cmake_minimum_required(VERSION 3.20)
project(untitled11)

set(CMAKE_CXX_STANDARD 14)
find_package(SFML 2.5 COMPONENTS graphics audio REQUIRED)
add_executable(untitled11 main.cpp)
target_link_libraries(untitled11 sfml-graphics sfml-audio)
