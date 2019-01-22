#!/bin/bash

echo "Compiling..."
g++ -O3 -Wall -Wextra -pedantic -fopenmp -fconcepts --std=c++17 -march=native -o nSimGram-BloomFilter nSimGram-BloomFilter.cpp
echo "Done!"

echo "Testing. H=1 Z=64"
cat NetInf.txt | ./nSimGram-BloomFilter --bruteforce --baseline --fsample -e 10 -r 50000 -A 0 -B 1 -t 4 -Q 6 -Z 64 -H 1
echo "Testing. H=2 Z=64"
cat NetInf.txt | ./nSimGram-BloomFilter --bruteforce --baseline --fsample -e 10 -r 50000 -A 0 -B 1 -t 4 -Q 6 -Z 64 -H 2
echo "Testing. H=4 Z=64"
cat NetInf.txt | ./nSimGram-BloomFilter --bruteforce --baseline --fsample -e 10 -r 50000 -A 0 -B 1 -t 4 -Q 6 -Z 64 -H 4
echo "Testing. H=8 Z=64"
cat NetInf.txt | ./nSimGram-BloomFilter --bruteforce --baseline --fsample -e 10 -r 50000 -A 0 -B 1 -t 4 -Q 6 -Z 64 -H 8
echo "Testing. H=16 Z=64"
cat NetInf.txt | ./nSimGram-BloomFilter --bruteforce --baseline --fsample -e 10 -r 50000 -A 0 -B 1 -t 4 -Q 6 -Z 64 -H 16
echo "Testing. H=32 Z=64"
cat NetInf.txt | ./nSimGram-BloomFilter --bruteforce --baseline --fsample -e 10 -r 50000 -A 0 -B 1 -t 4 -Q 6 -Z 64 -H 32
#cat NetInf.txt | ./nSimGram-BloomFilter --bruteforce --baseline --fsample -e 10 -r 5000 -A 0 -B 1 -t 1 --verbose
#cat NetInf.txt | ./nSimGram-BloomFilter --bruteforce --baseline --fsample -e 10 -r 5000 -A 0 -B 1 -t 1 --verbose
#cat NetInf.txt | ./nSimGram-BloomFilter --bruteforce --baseline --fsample -e 10 -r 5000 -A 0 -B 1 -t 1 --first --verbose
