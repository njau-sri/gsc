#!/bin/bash

rm -rf gsc-$1
mkdir gsc-$1

if [ $1 == "lnx64" ]; then

    g++ *.cpp -o gsc-$1/gsc -s -O2 -std=c++11 -static -llapack -lrefblas -lgfortran -lquadmath

elif [ $1 == "win32" ]; then

    i686-w64-mingw32-g++ *.cpp -o gsc-$1/gsc.exe -s -O2 -std=c++11 -static -llapack -lrefblas -lgfortran -lquadmath

elif [ $1 == "win64" ]; then

    x86_64-w64-mingw32-g++ *.cpp -o gsc-$1/gsc.exe -s -O2 -std=c++11 -static -llapack -lrefblas -lgfortran -lquadmath

fi
