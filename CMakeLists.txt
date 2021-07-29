cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(GIFS VERSION 0.2.0 LANGUAGES C CXX)

include(CTest)
enable_Testing()

include_directories(include)

add_subdirectory(gifs_src)
