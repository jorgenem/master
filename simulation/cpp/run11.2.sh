#!/bin/bash
rm main11.2
g++ -std=c++11 minimization11.2_MD-fit_event-pairing-preminimization.cpp -o main11.2 -O3 -I/usr/local/include  -I/usr/local/include -larmadillo
./main11.2