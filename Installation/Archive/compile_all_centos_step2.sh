cmake .. -DGMX_FFT_LIBRARY=fftw3 -DFFTWF_LIBRARY=/opt/fftw3/3.3.9/float-sse2-avx2/lib/libfftw3f.so -DFFTWF_INCLUDE_DIR=/opt/fftw3/3.3.9/float-sse2-avx2/include -DGMX_GPU=ON -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda -DCMAKE_INSTALL_PREFIX=/opt/gromacs/2020.6
make -j32
make install
cd ~
mkdir -p install && mv gromacs-2020.6.tar.gz install && rm -rf gromacs-2020.6

tar zxvf gcc-11.1.0.tar.gz
cd gcc-11.1.0
mkdir build
cd build
../configure --enable-language=c,c++ --disable-multilib --prefix=/opt/gcc/11.1.0
make -j32
make install
cd ~
mkdir -p install && mv gcc-11.1.0.tar.gz install && rm -rf gcc-11.1.0
