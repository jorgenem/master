#!/bin/bash
rm main9.4
g++ -std=c++11 minimization9.4_fit_mass_differences_with_dilepton_edge_multiple_combos_summed-detAcut.cpp -o main9.4 -O3 -I/usr/local/include  -I/usr/local/include -larmadillo
./main9.4