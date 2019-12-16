mkdir -p build/
mkdir -p install
cd build
CC=icc cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/../install -DGMX_GIFS=/data/menger/gromacs/gifs -DGMX_BUILD_OWN_FFTW=ON
