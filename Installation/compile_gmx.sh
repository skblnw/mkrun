#!/bin/bash

for ii in 4; do
  cd gromacs-2020.$ii
  rm -rf build
  mkdir build
  cd build
  cmake .. -DGMX_DOUBLE=ON -DCMAKE_INSTALL_PREFIX=/opt/gromacs/2020.$ii -DGMX_BUILD_OWN_FFTW=ON
  make -j24
  make install
  cd ../..
done
