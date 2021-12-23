cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(GIFS VERSION 0.2.0 LANGUAGES C CXX)

include(CTest)
enable_Testing()

if(${CMAKE_VERSION} VERSION_LESS "3.13")
    set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
    set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
endif()

include_directories(include)

add_subdirectory(gifs_src)
