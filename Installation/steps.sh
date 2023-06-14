#!/bin/bash

# Set default values for software versions and installation paths
namd_version="3.0alpha9"
vmd_version="1.9.4a51"
cmake_version="3.20.2"
gdrcopy_version="2.3"
ucx_version="1.12.1"
openmpi_version="4.1.4"
plumed_version="2.8.2"
fftw_version="3.3.9"
gromacs_version="2023.1"
install_base="/opt"

# Set default values for execution control
install_namd=false
install_vmd=false
install_cmake=false
install_gdrcopy=false
install_ucx=false
install_openmpi=false
install_fftw=false
install_plumed=false
install_gromacs=false
install_gromacs_mpi=false
install_gromacs_double=false

if $ubuntu; then
  apt install libssl-dev
  apt install dracut-core
elif $centos; then
  yum install centos-release-scl
  yum install devtoolset-7-gcc*
  yum install openssl-devel
  scl enable devtoolset-7 bash
fi

check_files() {
  local files=("$@")
  for file in "${files[@]}"; do
    if [ ! -f "$file" ]; then
      echo "Error: $file does not exist"
      exit 1
    fi
  done
}

check_directories() {
  local directories=("$@")
  for directory in "${directories[@]}"; do
    if [ ! -d "$directory" ]; then
      echo "Error: $directory does not exist"
      exit 1
    fi
  done
}

# Install NAMD
if $install_namd; then
  tar xvf "NAMD_${namd_version}_Linux-x86_64-multicore-CUDA.tar"
  mkdir "${install_base}/namd"
  mv "NAMD_${namd_version}_Linux-x86_64-multicore-CUDA" "${install_base}/namd/${namd_version}"
fi 

# Install VMD
if $install_vmd; then
  tar xvf vmd-${vmd_version}.bin.LINUXAMD64.opengl.tar
  cd vmd-${vmd_version}
  sed -i "s/^\$install_bin_dir=.*/\$install_bin_dir=\"${install_base}\/vmd\/${vmd_version}\";/g" configure
  sed -i "s/^\$install_library_dir=.*/\$install_library_dir=\"${install_base}\/vmd\/${vmd_version}\/lib\";/g" configure
  ./configure
  cd src
  make -j8 install
  cd ~
  rm -r vmd-${vmd_version}
fi

# Install CMake
if $install_cmake; then
  files=("CMake-${cmake_version}.tar.gz")
  check_files "${files[@]}"

  tar zxvf "CMake-${cmake_version}.tar.gz"
  cd "CMake-${cmake_version}"
  ./bootstrap
  make -j8
  make -j8 install
  cd ~
  rm -r "CMake-${cmake_version}"
fi

# Install gdrcopy
if $install_gdrcopy; then 
  files=("gdrcopy-${gdrcopy_version}.tar.gz")
  check_files "${files[@]}"
  directories=("/usr/local/cuda")
  check_directories "${directories[@]}"

  # apt install check libsubunit0 libsubunit-dev
  tar zfx "gdrcopy-${gdrcopy_version}.tar.gz"
  cd "gdrcopy-${gdrcopy_version}"
  make CUDA="/usr/local/cuda" all install
  ./insmod.sh
  cd ~
  rm -r "gdrcopy-${gdrcopy_version}"
fi

# Install UCX
if $install_ucx; then
  files=("ucx-${ucx_version}.tar.gz")
  check_files "${files[@]}"
  directories=("/usr/local/cuda")
  check_directories "${directories[@]}"

  tar zfxv "ucx-${ucx_version}.tar.gz"
  cd "ucx-${ucx_version}"
  #./autogen.sh
  ./contrib/configure-release --prefix="${install_base}/ucx/${ucx_version}" --with-cuda="/usr/local/cuda"
  make -j8
  make -j8 install
  cd ~
  rm -r "ucx-${ucx_version}"
fi

# Install OpenMPI
if $install_openmpi; then
  files=("openmpi-${openmpi_version}.tar.gz")
  check_files "${files[@]}"
  directories=("/usr/local/cuda" "/opt/ucx/${ucx_version}")
  check_directories "${directories[@]}"

  tar zfxv "openmpi-${openmpi_version}.tar.gz"
  cd "openmpi-${openmpi_version}"
  mkdir build
  cd build
  ../configure --prefix="${install_base}/openmpi/${openmpi_version}" --with-cuda="/usr/local/cuda" --with-ucx="${install_base}/ucx/${ucx_version}"
  make -j8 all
  make -j8 install
  cd ~
  rm -r "openmpi-${openmpi_version}"
fi

# Install FFTW
if $install_fftw; then
  files=("fftw-${fftw_version}.tar.gz")
  check_files "${files[@]}"

  tar zxvf "fftw-${fftw_version}.tar.gz"
  cd "fftw-${fftw_version}"
  ./configure --enable-float --enable-sse2 --enable-avx2 --enable-shared --prefix=/opt/fftw3/${fftw_version}/float-sse2-avx2
  make -j8
  make -j8 install
  cd ~
  rm -r "fftw-${fftw_version}"
fi

# Install FFTW (MPI)
if $install_fftw_mpi; then
  files=("fftw-${fftw_version}.tar.gz")
  check_files "${files[@]}"

  tar zxvf "fftw-${fftw_version}.tar.gz"
  cd "fftw-${fftw_version}"
  ./configure --enable-float --enable-sse2 --enable-avx2 --enable-threads --enable-openmp --enable-mpi --enable-shared \
              MPICC= \
              MPICXX= \
              --prefix=/opt/fftw3/${fftw_version}/float-sse2-avx2-omp-mpi
  make -j8
  make -j8 install
  cd ~
  rm -r "fftw-${fftw_version}"
fi

# Install Plumed
if $install_plumed; then
  files=("plumed-${plumed_version}.tar" "gromacs-${gromacs_version}.tar")
  check_files "${files[@]}"
  directories=("/usr/local/cuda" "/opt/openmpi/${openmpi_version}")
  check_directories "${directories[@]}"

  tar zfxv "plumed-${plumed_version}.tar"
  cd "plumed-${plumed_version}"
  ./configure --prefix=/opt/plumed/${plumed_version}
  make -j8
  make -j8 install
  cd ~
  rm -r "plumed-${plumed_version}"

  export PATH=/opt/plumed/${plumed_version}/bin:$PATH
  export LD_LIBRARY_PATH=/opt/plumed/${plumed_version}/lib/:$LD_LIBRARY_PATH
  export PLUMED_KERNEL=/opt/plumed/${plumed_version}/lib/libplumedKernel.so

  tar xvf "gromacs-${gromacs_version}.tar"
  cd "gromacs-${gromacs_version}"
  plumed patch -p
  # Comment out compute_30 35 37 50

  mkdir build
  cd build
  cmake .. -DGMX_FFT_LIBRARY=fftw3 \
            -DFFTWF_LIBRARY=/opt/fftw3/${fftw_version}/float-sse2-avx2/lib/libfftw3f.so \
            -DFFTWF_INCLUDE_DIR=/opt/fftw3/${fftw_version}/float-sse2-avx2/include \
            -DGMX_GPU=CUDA \
            -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
            -DCMAKE_INSTALL_PREFIX=/opt/gromacs/${gromacs_version}-plumed \
            -DGMX_MPI=ON \
            -DMPI_CXX_COMPILER=/opt/openmpi/${openmpi_version}/bin/mpicxx \
            -DMPI_C_COMPILER=/opt/openmpi/${openmpi_version}/bin/mpicc
  make -j8
  make -j8 install
  cd ~
  rm -r "gromacs-${gromacs_version}"
fi

# Install Gromacs
if $install_gromacs; then
  files=("gromacs-${gromacs_version}.tar.gz")
  check_files "${files[@]}"
  directories=("/usr/local/cuda" "/opt/fftw3/${fftw_version}")
  check_directories "${directories[@]}"

  tar xvf "gromacs-${gromacs_version}.tar.gz"
  cd "gromacs-${gromacs_version}"
  mkdir build
  cd build
  cmake .. -DGMX_FFT_LIBRARY=fftw3 \
            -DFFTWF_LIBRARY=/opt/fftw3/${fftw_version}/float-sse2-avx2/lib/libfftw3f.so \
            -DFFTWF_INCLUDE_DIR=/opt/fftw3/${fftw_version}/float-sse2-avx2/include \
            -DGMX_GPU=CUDA \
            -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
            -DCMAKE_INSTALL_PREFIX=/opt/gromacs/${gromacs_version} \
            -DCMAKE_CXX_COMPILER= \
            -DCMAKE_C_COMPILER= 
  make -j8
  make -j8 install
  cd ~
  rm -r "gromacs-${gromacs_version}"
fi

# Install Gromacs (MPI version)
if $install_gromacs_mpi; then
  files=("gromacs-${gromacs_version}.tar.gz")
  check_files "${files[@]}"
  directories=("/usr/local/cuda" "/opt/fftw3/${fftw_version}" "/opt/openmpi/${openmpi_version}")
  check_directories "${directories[@]}"

  tar xvf "gromacs-${gromacs_version}.tar.gz"
  cd "gromacs-${gromacs_version}"
  mkdir build
  cd build
  cmake .. -DGMX_FFT_LIBRARY=fftw3 \
            -DFFTWF_LIBRARY=/opt/fftw3/${fftw_version}/float-sse2-avx2/lib/libfftw3f.so \
            -DFFTWF_INCLUDE_DIR=/opt/fftw3/${fftw_version}/float-sse2-avx2/include \
            -DGMX_GPU=CUDA \
            -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
            -DCMAKE_INSTALL_PREFIX=/opt/gromacs/${gromacs_version} \
            -DGMX_MPI=ON \
            -DMPI_CXX_COMPILER=/opt/openmpi/${openmpi_version}/bin/mpicxx \
            -DMPI_C_COMPILER=/opt/openmpi/${openmpi_version}/bin/mpicc
  make -j8
  make -j8 install
  cd ~
  rm -r "gromacs-${gromacs_version}"
fi

# Install Gromacs (double precision)
if $install_gromacs_double; then
  files=("gromacs-${gromacs_version}.tar.gz" "fftw-${fftw_version}.tar.gz")
  check_files "${files[@]}"

  tar zxvf "fftw-${fftw_version}.tar.gz"
  cd "fftw-${fftw_version}"
  ./configure --enable-sse2 --enable-avx2 --enable-shared --prefix=/opt/fftw3/${fftw_version}/double-sse2-avx2
  make -j8
  make -j8 install
  cd ~
  rm -r "fftw-${fftw_version}"

  tar xvf "gromacs-${gromacs_version}.tar.gz"
  cd "gromacs-${gromacs_version}"
  mkdir build
  cd build
  cmake ..  -DGMX_DOUBLE=ON \
            -DGMX_FFT_LIBRARY=fftw3 \
            -DFFTW_LIBRARY=/opt/fftw3/${fftw_version}/double-sse2-avx2/lib/libfftw3.so \
            -DFFTW_INCLUDE_DIR=/opt/fftw3/${fftw_version}/double-sse2-avx2/include \
            -DCMAKE_INSTALL_PREFIX=/opt/gromacs/${gromacs_version}
  make -j8
  make -j8 install
  cd ~
  rm -r "gromacs-${gromacs_version}"
fi
