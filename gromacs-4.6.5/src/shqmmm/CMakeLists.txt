set(CMAKE_POSITION_INDEPENDENT_CODE 1)

enable_language(CXX)

find_package(Armadillo)
if(${ARMADILLO_FOUND})
    include_directories(${ARMADILLO_INCLUDE_DIRS})
    link_libraries(${ARMADILLO_LIBRARIES})
else()
    include_directories(${GMX_GIFS}/armadillo/include)
endif()

set(GIFS_SOURCES
    "gmxgifs.c"
    "eshake.cpp"
)

add_library(gifsgmx
    SHARED
    ${GIFS_SOURCES}
)

                              
target_include_directories(gifsgmx
    PRIVATE
    ${GMX_GIFS}/include
)

target_link_libraries(gifsgmx
    PRIVATE
    ${GMX_GIFS}/lib/libinaqs.so
)

install(TARGETS gifsgmx 
        DESTINATION ${LIB_INSTALL_DIR} 
        COMPONENT libraries)
