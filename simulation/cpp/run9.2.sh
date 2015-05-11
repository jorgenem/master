#!/bin/bash
rm main9.2
g++ -std=c++11 minimization9.2_fit_mass_differences_with_dilepton_edge_several_SP.cpp -o main9.2 -O3 -I/usr/local/include  -I/usr/local/include -larmadillo
./main9.2