#!/usr/bin/env python
# (c) 2017 Alexander Kanitz, alexander.kanitz@alumni.ethz.ch

import sys

in_file = sys.argv[1]
replacement_file = sys.argv[2]
placeholder = sys.argv[3]

f = open(replacement_file)
replacement_list = f.readlines()
f.close()

with open(in_file) as f:
    counter = 0
    for line in f:
        line = line.rstrip()
        if placeholder in line:
            line = line.replace(placeholder, replacement_list[counter])
            counter += 1
        print(line)
