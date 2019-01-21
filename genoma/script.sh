#!/bin/bash

ulimit -s unlimited

g++ -Wall -g -Wextra -pedantic -pipe -march=native --std=c++17 -Ofast -fopenmp -o graph graph_fast_int.cpp graph.hpp
#g++ -Wall -g -Wextra -pedantic -pipe -march=native --std=c++17 -Ofast -fopenmp -o graph graph.cpp

#perf stat ./graph 12 8 100 < reads_cleaned.txt
#perf stat -B -e cache-references,cache-misses,cycles,instructions,branches,faults,migrations ./graph 12 8 100 < reads_cleaned.txt

#./graph 12 6 10 < reads_cleaned.txt

cp reads_final.txt /tmp/reads_final.txt
./graph < /tmp/reads_final.txt
