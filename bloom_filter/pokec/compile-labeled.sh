#!/bin/bash

echo "Compiling..."
g++ -O3 -Wall -Wextra -pedantic -fopenmp -fconcepts --std=c++17 -march=native -o nSimGram-BloomFilter-multilabeled nSimGram-BloomFilter-multilabeled.cpp
echo "Done!"

echo
echo
echo "30322 - 38063"
./nSimGram-BloomFilter-multilabeled --bruteforce -e 10 -r 10000 -A 30322 -B 38063 -t 4 -Q 3 -Z 64 -H 24 --verbose

echo
echo
echo "30322 - 41249"
./nSimGram-BloomFilter-multilabeled --bruteforce -e 10 -r 10000 -A 30322 -B 41249 -t 4 -Q 3 -Z 64 -H 24 --verbose

echo
echo
echo "30322 - 60611"
./nSimGram-BloomFilter-multilabeled --bruteforce -e 10 -r 10000 -A 30322 -B 60611 -t 4 -Q 3 -Z 64 -H 24 --verbose

echo
echo
echo "30322 - 77017"
./nSimGram-BloomFilter-multilabeled --bruteforce -e 10 -r 10000 -A 30322 -B 77017 -t 4 -Q 3 -Z 64 -H 24 --verbose

echo
echo
echo "30322 - 77637"
./nSimGram-BloomFilter-multilabeled --bruteforce -e 10 -r 10000 -A 30322 -B 77637 -t 4 -Q 3 -Z 64 -H 24 --verbose

########################################

echo
echo
echo "41249 - 30322"
./nSimGram-BloomFilter-multilabeled --bruteforce -e 10 -r 10000 -A 41249 -B 30322 -t 4 -Q 3 -Z 64 -H 24 --verbose

echo
echo
echo "41249 - 38063"
./nSimGram-BloomFilter-multilabeled --bruteforce -e 10 -r 10000 -A 41249 -B 38063 -t 4 -Q 3 -Z 64 -H 24 --verbose

echo
echo
echo "41249 - 60611"
./nSimGram-BloomFilter-multilabeled --bruteforce -e 10 -r 10000 -A 41249 -B 60611 -t 4 -Q 3 -Z 64 -H 24 --verbose

echo
echo
echo "41249 - 77017"
./nSimGram-BloomFilter-multilabeled --bruteforce -e 10 -r 10000 -A 41249 -B 77017 -t 4 -Q 3 -Z 64 -H 24 --verbose

echo
echo
echo "41249 - 77637"
./nSimGram-BloomFilter-multilabeled --bruteforce -e 10 -r 10000 -A 41249 -B 77637 -t 4 -Q 3 -Z 64 -H 24 --verbose


#30322 Giuseppe--Di--Battista
#38063 J.--Ian--Munro
#41249 Jeffrey--Scott--Vitter
#60611 Maurizio--Patrignani
#77017 Robert--Endre--Tarjan
#77637 Roberto--Grossi

