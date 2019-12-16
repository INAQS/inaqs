set(CMAKE_POSITION_INDEPENDENT_CODE 1)
#
file(GLOB GIFS_SOURCES *.c)

add_library(gifsgmx
    SHARED
    ${GIFS_SOURCES}
)

#set_property(TARGETS pyqmmmgromacs 
#    PROPERTY 
#      POSITION_INDEPENDENT_CODE 1
#)

target_include_directories(gifsgmx
    PRIVATE
)

target_link_libraries(gifsgmx
    PRIVATE
)

install(TARGETS gifsgmx 
        DESTINATION ${LIB_INSTALL_DIR} 
        COMPONENT libraries)
