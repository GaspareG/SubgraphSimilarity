#!/bin/bash

echo "Compiling..."
g++ -O3 -Wall -Wextra -pedantic -fopenmp -fconcepts --std=c++17 -march=native -o nSimGram-BloomFilter-multilabeled nSimGram-BloomFilter-multilabeled.cpp
echo "Done!"

echo
./nSimGram-BloomFilter-multilabeled --bruteforce --baseline --fsample -e 10 -r 100 -t 4 -Q 3 -Z 32 -H 1 --verbose
