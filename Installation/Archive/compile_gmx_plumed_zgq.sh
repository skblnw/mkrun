#!/bin/bash

tar -zfxv ucx-1.12.1.tar.gz
cd ucx-1.12.1
./contrib/configure-release --prefix=/public/home/zgq/software/package/ucx-1.12.1/install --with-cuda=/usr/local/cuda-11.4
make -j8
make install
echo "export PATH=/public/home/zgq/software/package/ucx-1.12/bin:$PATH" >> ~/.bashrc

cd ..
tar -zxvf openmpi-4.1.4.tar.gz
cd openmpi-4.1.4
mkdir build
cd build
../configure --prefix=/public/home/zgq/software/package/openmpi/v4.1.4 --with-cuda=/usr/local/cuda-11.4 --with-ucx=/public/home/zgq/software/package/ucx-1.12.1/install

make all
make -j8 install
echo "export PATH=/public/home/zgq/software/package/openmpi/v4.1.4/bin:$PATH" >> ~/.bashrc

#Than install plumed2-2.8.0 and gromacs2022 with MPI-version

use devtools to upgrade gcc

source /opt/rh/devtoolset-7/enable
(gmx 2022版本都需要gcc很高的级别)

download FFTW
GROMACS 2018依赖于快速傅立叶变换库FFTW 3.3.8，可以在http://www.fftw.org/fftw-3.3.8.tar.gz下载。将其压缩包解压，进入此目录后运行
./configure --prefix=/home/zgq/software/fftw338 --enable-sse2 --enable-avx --enable-float --enable-shared
以上语句代表FFTW将被安装到/sob/fftw338目录。如果你的CPU相对较新，支持AVX2指令集，可再加上--enable-avx2选项以获得更好性能。

然后运行make -j install开始编译，过一会儿编译完毕后，就出现了/sob/fftw338目录。然后可以把FFTW的解压目录和压缩包删掉了

gmx 按照教程
tar xfz gromacs-2021.4.tar.gz
mkdir build
cd build
cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX=~/software/GMX2022/2022.3 -DGMX_FFT_LIBRARY=fftw3 -DFFTWF_LIBRARY=~/software/fftw338/lib/libfftw3f.so -DFFTWF_INCLUDE_DIR=~/software/fftw338/include/ -DREGRESSIONTEST_PATH=~/software/regressiontests-2022.3
make -j 16
make -j 16 check
make -j 16 install

plumed download
./configure --disable-mpi --prefix=/home/zgq/software/plumed/
make -j 16
make -j 16 install
进入gromacs2022.3文件夹（源代码所在地）
plumed patch -p

then recompile the gromacs

make -j 32
make -j 32 check
make -j 32 install
