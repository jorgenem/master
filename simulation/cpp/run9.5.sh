#!/bin/bash
rm main9.5
g++ -std=c++11 minimization9.5_fit_mass_differences_with_dilepton_edge_multiple_combos_summed-subdetAcut.cpp -o main9.5 -O3 -I/usr/local/include  -I/usr/local/include -larmadillo
./main9.5