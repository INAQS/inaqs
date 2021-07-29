find_package(Armadillo 10)
if(${ARMADILLO_FOUND})
    include_directories(${ARMADILLO_INCLUDE_DIRS})
    link_libraries(${ARMADILLO_LIBRARIES})
else()
    include_directories(../armadillo/include)
endif()

find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
link_libraries(${HDF5_LIBRARIES})

add_library(gifs SHARED
    bomd.cpp
    bomd_rescale.cpp
    configreader.cpp
    conversion.cpp
    decoherence_afssh.cpp
    electronic.cpp
    fssh.cpp
    gifs.cpp
    gifs_implementation.cpp
    linkatoms.cpp
    properties.cpp
    qm_interface.cpp
    qm_model.cpp
    qm_qchem.cpp
    util.cpp
)

target_compile_options(gifs PRIVATE -DARMA_DONT_USE_WRAPPER -DARMA_USE_HDF5)
target_compile_options(gifs PRIVATE -Werror -Wall -Wextra)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set_target_properties(gifs PROPERTIES VERSION ${PROJECT_VERSION})
set_target_properties(gifs PROPERTIES PUBLIC_HEADER ../include/gifs.hpp)

install(TARGETS gifs
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
