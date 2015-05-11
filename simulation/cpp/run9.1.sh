#!/bin/bash
rm main9.1
g++ -std=c++11 minimization9.1_fit_mass_differences_with_dilepton_edge.cpp -o main9.1 -O3 -I/usr/local/include  -I/usr/local/include -larmadillo
./main9.1