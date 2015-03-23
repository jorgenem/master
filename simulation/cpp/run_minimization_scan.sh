#!/bin/bash
rm main_ms
g++ -std=c++11 minimization6.1_combinatorics-count_avoid_unphysical_regions_minimization_scan.cpp -o main_ms -O3 -I/usr/local/include  -I/usr/local/include -larmadillo
./main_ms $1 $2 