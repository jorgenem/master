#!/bin/bash
rm main9
g++ -std=c++11 minimization9_fit_mass_differences.cpp -o main9 -O3 -I/usr/local/include  -I/usr/local/include -larmadillo
./main9