#!/bin/bash
rm main11
g++ -std=c++11 minimization11_MD-fit_all_combos_separately.cpp -o main11 -O3 -I/usr/local/include  -I/usr/local/include -larmadillo
./main11