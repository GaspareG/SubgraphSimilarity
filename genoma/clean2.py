#!/bin/python

import re
import sys

'''
@chrI_91407_91550_0_1_0_0_0:0:0_1:0:0_0/1
GATCCTCTTGAACCTCCTGGAAATATCACATATAGTAGTTCGAATCTTTCGCTAAATTCAGATGAATTAGACTACTATCAGCGTCATATCGGATTGCAGTTACAGCAGACAGAAGCTTTACTAAAGCACAGTTTGAAAGATGAGGTTCTG
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
'''

def read_in():
    lines = sys.stdin.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].replace('\n','')
    return lines

lines = read_in()

for i in range(0, len(lines)):
  info = lines[i]
  if info[0] == 'A' or info[0] == 'C' or info[0] == 'G' or info[0] == 'T':
    print(info)
