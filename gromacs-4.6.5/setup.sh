mkdir -p build/
mkdir -p install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/../install  -DGMX_BUILD_OWN_FFTW=ON
