current=$PWD
echo $current
mkdir -p build
cd build
cmake .. -DGMX_DOUBLE=off -DGMX_CPU_ACCELERATION=AVX_256 -DGMX_OPENMP=ON -DGMX_GIFS=$current/../ -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX=$current/install
make -j 64
make -j 64 install
