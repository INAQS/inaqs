set(CMAKE_POSITION_INDEPENDENT_CODE 1)
#
file(GLOB GIFS_SOURCES *.c)

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
    ${GMX_GIFS}/lib/libgifs.so
)

install(TARGETS gifsgmx 
        DESTINATION ${LIB_INSTALL_DIR} 
        COMPONENT libraries)
