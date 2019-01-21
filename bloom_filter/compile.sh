#!/bin/bash

echo "Compiling..."
g++ -O3 -Wall -Wextra -pedantic -fopenmp --std=c++17 -march=native -o nSimGram-BloomFilter nSimGram-BloomFilter.cpp
echo "Done!"

echo "Testing..."
cat NetInf.txt | ./nSimGram-BloomFilter --bruteforce --baseline --fsample -e 10 -r 5000 -A 0 -B 1 -t 1 --verbose
cat NetInf.txt | ./nSimGram-BloomFilter --bruteforce --baseline --fsample -e 10 -r 5000 -A 0 -B 1 -t 1 --first --verbose
