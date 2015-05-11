#!/bin/bash
rm main10
g++ -std=c++11 minimization10_write_determinant.cpp -o main10 -O3 -I/usr/local/include  -I/usr/local/include -larmadillo
./main10