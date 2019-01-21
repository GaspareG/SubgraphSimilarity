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
    #print lines
    return lines

lines = read_in()

print( str(len(lines)>>2) )

for i in range(0, len(lines), 4):
  info = lines[i]
  dig = re.findall(r'\d+', info)
  base = lines[i+1]
  print(dig[0] + " " + dig[1])
  print(base)
