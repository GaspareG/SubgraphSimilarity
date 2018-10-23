#!/bin/bash

g++ -O3 -Wall -Wextra -pedantic -fopenmp --std=c++17 -march=native -o nSimGram-BloomFilter nSimGram-BloomFilter.cpp
