#!/bin/bash

MAINDIR=`pwd`

wget http://www.fftw.org/fftw-3.3.3.tar.gz

tar -xzvf fftw-3.3.3.tar.gz

cd fftw-3.3.3


FFTWPATH=$MAINDIR/fftw-3.3.3/build

./configure CFLAGS="-fPIC" --enable-threads --prefix=$FFTWPATH --enable-float

make
make install
make clean

./configure CFLAGS="-fPIC" --enable-threads --prefix=$FFTWPATH 

make 
make install
make clean

rm ../fftw-3.3.3.tar.gz

cd $MAINDIR

wget -O InsightToolkit-4.4.1.tar.gz http://sourceforge.net/projects/itk/files/itk/4.4/InsightToolkit-4.4.1.tar.gz/download/

tar -xzvf InsightToolkit-4.4.1.tar.gz


mkdir ./InsightToolkit-4.4.1/build

ITKPATH=$MAINDIR/InsightToolkit-4.4.1/build
cd $ITKPATH/

make clean

cmake \
-D BUILD_DOCUMENTATION:BOOL=OFF \
-D BUILD_TESTING:BOOL=ON \
-D BUILD_EXAMPLES:BOOL=OFF \
-D CMAKE_CXX_FLAGS:STRING="-fPIC" \
-D ITK_USE_FFTWF:BOOL=ON \
-D FFTWF_LIB:FILEPATH="$FFTWPATH/lib/libfftw3f.a" \
-D FFTWF_THREADS_LIB:FILEPATH="$FFTWPATH/lib/libfftw3f_threads.a" \
-D ITK_USE_FFTWD:BOOL=ON \
-D FFTWD_LIB:FILEPATH="$FFTWPATH/lib/libfftw3.a" \
-D FFTWD_THREADS_LIB:FILEPATH="$FFTWPATH/lib/libfftw3_threads.a" \
-D FFTW_INCLUDE_PATH:FILEPATH="$FFTWPATH/include" \
-D ITK_USE_REVIEW:BOOL= OFF \
-D ITK_USE_SYSTEM_FFTW:BOOL=ON \
-D BUILD_SHARED_LIBS:BOOL=OFF \
../

make -j 10

rm InsightToolkit-4.4.1.tar.gz

cd $MAINDIR

LCClogDemonsPATH=$MAINDIR/LCClogDemonsV1.0/build
mkdir $LCClogDemonsPATH
cd $LCClogDemonsPATH

cmake \
-D CMAKE_CXX_FLAGS:STRING="-fPIC" \
-D ITK_DIR:FILEPATH="$ITKPATH" \
../src

make -j 10
cd $MAINDIR



