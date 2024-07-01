cd NAMD_2.14_Source
tar xf charm-6.10.2.tar
cd charm-6.10.2
./build charm++ multicore-linux-x86_64 --with-production
cd ..
wget http://www.ks.uiuc.edu/Research/namd/libraries/fftw-linux-x86_64.tar.gz
tar xzf fftw-linux-x86_64.tar.gz
mv linux-x86_64 fftw
wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.5.9-linux-x86_64.tar.gz
wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.5.9-linux-x86_64-threaded.tar.gz
tar xzf tcl8.5.9-linux-x86_64.tar.gz
tar xzf tcl8.5.9-linux-x86_64-threaded.tar.gz
mv tcl8.5.9-linux-x86_64 tcl
mv tcl8.5.9-linux-x86_64-threaded tcl-threaded
cp arch/Linux-x86_64-g++.arch arch/Linux-x86_64-g++.arch.BAK
cat > arch/Linux-x86_64-g++.arch <<'EOF'
NAMD_ARCH = Linux-x86_64
CHARMARCH = multicore-linux-x86_64

CXX = g++ -m64 -std=c++0x
CXXOPTS = -O3 -fexpensive-optimizations -ffast-math -no-pie
CC = gcc -m64
COPTS = -O3 -fexpensive-optimizations -ffast-math
EOF
./config Linux-x86_64-g++ --charm-arch multicore-linux-x86_64 --with-cuda
cd Linux-x86_64-g++
#make -j4
cd ..
mv Linux-x86_64-g++ 2.14