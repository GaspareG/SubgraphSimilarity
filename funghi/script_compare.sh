#!/bin/bash

g++ -O3 -Wall -Wextra --pedantic -march=native -fopenmp -lm --std=c++17 -static -pipe -o compare compare.cpp

FUNGHI=(ERR1308675 ERR1308682 ERR1308727 ERR1308789 ERR1308797 ERR1308906 ERR1308934 ERR1308982 ERR1308996 ERR1309034 ERR1309255 ERR1309259 ERR1309341 ERR1309458)

RESULT="results/"
SAMPLING="sampling/"

for i in "${FUNGHI[@]}"
do
  for j in  "${FUNGHI[@]}"
  do
    ./compare "${SAMPLING}${i}.txt" "${SAMPLING}${j}.txt" > "${RESULT}${i}_${j}.txt"
  done
done
