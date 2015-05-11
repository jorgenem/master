#!/bin/bash
rm main11.1
g++ -std=c++11 minimization11.1_MD-fit_all_combos_separately_with_dilepton_edge.cpp -o main11.1 -O3 -I/usr/local/include  -I/usr/local/include -larmadillo
./main11.1