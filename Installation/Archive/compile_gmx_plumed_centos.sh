tar zxvf fftw-3.3.9.tar.gz
cd fftw-3.3.9
./configure --enable-float --enable-sse2 --enable-avx2 --enable-shared --prefix=/opt/fftw3/3.3.9/float-sse2-avx2
make -j32
make install -j32
cd ~
mkdir -p install && mv fftw-3.3.9.tar.gz install && rm -rf fftw-3.3.9

tar zxvf mpich-3.4.2.tar.gz
cd mpich-3.4.2/
./configure --prefix=/opt/mpich/3.4.2-ch3 --with-device=ch3 --disable-fortran
make -j32
make install -j32
cd ~
mkdir -p install && mv mpich-3.4.2.tar.gz install && rm -rf mpich-3.4.2

export PATH=/opt/mpich/3.4.2/bin:$PATH
export LD_LIBRARY_PATH=/opt/mpich/3.4.2/lib/:$LD_LIBRARY_PATH
tar zxvf plumed-2.6.3.tgz
cd plumed-2.6.3/
./configure --prefix=/opt/plumed/2.6.3
make -j32
make install -j32
cd ~
mkdir -p install && mv plumed-2.6.3.tgz install && rm -rf plumed-2.6.3

export PATH=/opt/plumed/2.6.3/bin:$PATH
export LD_LIBRARY_PATH=/opt/plumed/2.6.3/lib/:$LD_LIBRARY_PATH
export PLUMED_KERNEL=/opt/plumed/2.6.3/lib/libplumedKernel.so
tar xvf gromacs-2019.6.tar.gz
cd gromacs-2019.6/
/opt/plumed/2.6.3/bin/plumed patch -p
scl enable devtoolset-7 bash
# Comment out compute_30 35 37 50

mkdir build
cd build
cmake .. -DGMX_FFT_LIBRARY=fftw3 -DFFTWF_LIBRARY=/opt/fftw3/3.3.9/float-sse2-avx2/lib/libfftw3f.so -DFFTWF_INCLUDE_DIR=/opt/fftw3/3.3.9/float-sse2-avx2/include -DGMX_GPU=ON -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda -DCMAKE_INSTALL_PREFIX=/opt/gromacs/2019.6-plumed -DGMX_MPI=ON -DMPI_CXX_COMPILER=/opt/mpich/3.4.2/bin/mpicxx -DMPI_C_COMPILER=/opt/mpich/3.4.2/bin/mpicc
make -j32
make install -j32
cd ~
mkdir -p install && mv gromacs-2019.6.tar.gz install && rm -rf gromacs-2019.6
