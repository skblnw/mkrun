apt install libssl-dev
apt install dracut-core

tar zxvf CMake-3.20.2.tar.gz
cd CMake-3.20.2
./bootstrap
make -j16
make install
cd ~
mkdir -p install && mv CMake-3.20.2.tar.gz install && rm -r CMake-3.20.2

tar zxvf fftw-3.3.9.tar.gz
cd fftw-3.3.9
./configure --enable-float --enable-sse2 --enable-avx2 --enable-shared --prefix=/opt/fftw3/3.3.9/float-sse2-avx2
make -j16
make install
cd ~
mkdir -p install && mv fftw-3.3.9.tar.gz install && rm -r fftw-3.3.9

# tar xvf NAMD_3.0alpha9_Linux-x86_64-multicore-CUDA.tar
# mkdir /opt/namd
# mv NAMD_3.0alpha9_Linux-x86_64-multicore-CUDA /opt/namd/3.0alpha9
# cd ~
# mkdir -p install && mv NAMD_3.0alpha9_Linux-x86_64-multicore-CUDA.tar install

# tar xvf vmd-1.9.4a51.bin.LINUXAMD64.opengl.tar
# cd vmd-1.9.4a51
# sed -i 's/^\$install_bin_dir=.*/\$install_bin_dir=\"\/opt\/vmd\/1.9.4\"\;/g' configure
# sed -i 's/^\$install_library_dir=.*/\$install_library_dir=\"\/opt\/vmd\/1.9.4\/lib\"\;/g' configure
# ./configure
# cd src
# make install
# cd ~
# mkdir -p install && mv vmd-1.9.4a51.bin.LINUXAMD64.opengl.tar install && rm -r vmd-1.9.4a51

# tar zxvf gromacs-2020.6.tar.gz
# cd gromacs-2020.6
# mkdir build
# cd build
# cmake .. -DGMX_FFT_LIBRARY=fftw3 -DFFTWF_LIBRARY=/opt/fftw3/3.3.9/float-sse2-avx2/lib/libfftw3f.so -DFFTWF_INCLUDE_DIR=/opt/fftw3/3.3.9/float-sse2-avx2/include -DGMX_GPU=ON -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda -DCMAKE_INSTALL_PREFIX=/opt/gromacs/2020.6
# make -j32
# make install
# cd ~
# mkdir -p install && mv gromacs-2020.6.tar.gz install && rm -r gromacs-2020.6

tar zxvf gromacs-2023.tar.gz
cd gromacs-2023
mkdir build
cd build
cmake .. -DGMX_FFT_LIBRARY=fftw3 -DFFTWF_LIBRARY=/opt/fftw3/3.3.9/float-sse2-avx2/lib/libfftw3f.so -DFFTWF_INCLUDE_DIR=/opt/fftw3/3.3.9/float-sse2-avx2/include -DGMX_GPU=CUDA -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda -DCMAKE_INSTALL_PREFIX=/opt/gromacs/2023
make -j16
make install -j16
cd ~
mkdir -p install && mv gromacs-2023.tar.gz install/ && rm -r gromacs-2023
