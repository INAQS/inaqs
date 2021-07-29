find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(test-gifs
    main.cpp
    electronic.cpp)

target_link_libraries(test-gifs gifs)
target_link_options(test-gifs PRIVATE ${BLAS_LIBRARIES} ${BLAS_LINKER_FLAGS})
target_link_options(test-gifs PRIVATE ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS})

target_include_directories(test-gifs PRIVATE ..)

add_test(NAME all-tests COMMAND test-gifs --success)
