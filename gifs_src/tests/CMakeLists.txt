find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(test-gifs
    main.cpp
    electronic.cpp)

target_link_libraries(test-gifs gifs)

if(${CMAKE_VERSION} VERSION_LESS "3.13")
    target_link_libraries(test-gifs ${BLAS_LIBRARIES} ${BLAS_LINKER_FLAGS})
    target_link_libraries(test-gifs ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS})
    target_link_libraries(test-gifs -fopenmp)
else()
    target_link_options(test-gifs PRIVATE ${BLAS_LIBRARIES} ${BLAS_LINKER_FLAGS})
    target_link_options(test-gifs PRIVATE ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS})
endif()

target_include_directories(test-gifs PRIVATE ..)

add_test(NAME all-tests COMMAND test-gifs --success)
