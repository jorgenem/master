#!/bin/bash
rm main9.2.1
g++ -std=c++11 minimization9.2.1-combosum-several_SP.cpp -o main9.2.1 -O3 -I/usr/local/include  -I/usr/local/include -larmadillo
./main9.2.1