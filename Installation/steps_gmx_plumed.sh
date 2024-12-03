#!/bin/bash
# GCC must >= 9.5.0
# installation sequence:
## 1. gdrcopy
## 2. ucx
## 3. openmpi
## 4. plumed
## 5. gmx

# Set default values for software versions and installation paths
cmake_version="3.20.2"
gdrcopy_version="2.3.1"
ucx_version="1.12.1"
openmpi_version="4.1.4"
plumed_version="2.9.0"
fftw_version="3.3.10"
gromacs_version="2022.5"
install_base="/opt"

# Set default values for execution control
install_cmake=false
install_gdrcopy=false
install_ucx=false
install_openmpi=false
install_fftw_mpi=false
install_plumed=false

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
