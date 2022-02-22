find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(test-inaqs
    main.cpp
    util.cpp
    electronic.cpp)

target_link_libraries(test-inaqs inaqs)

if(${CMAKE_VERSION} VERSION_LESS "3.13")
    target_link_libraries(test-inaqs ${BLAS_LIBRARIES} ${BLAS_LINKER_FLAGS})
    target_link_libraries(test-inaqs ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS})
    target_link_libraries(test-inaqs -fopenmp)
else()
    target_link_options(test-inaqs PRIVATE ${BLAS_LIBRARIES} ${BLAS_LINKER_FLAGS})
    target_link_options(test-inaqs PRIVATE ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS})
endif()

target_include_directories(test-inaqs PRIVATE ..)

add_test(NAME all-tests COMMAND test-inaqs --success)
