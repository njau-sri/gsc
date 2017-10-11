#!/bin/bash

g++ *.cpp -o gsc -s -O2 -std=c++11 -static -llapack -lblas -lgfortran -lquadmath
