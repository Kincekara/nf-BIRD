#!/usr/bin/env python3
"""
This script trims a fasta file containing scaffolds by removing scaffolds shorter than a given threshold.
@author: Kutluhan Incekara
email: kutluhan.incekara@ct.gov
"""

import re
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
threshold = int(sys.argv[3])

# Copy lines from the input file to the output file until a scaffold shorter than the threshold is found
with open(input_file, "r") as input:
    with open(output_file, "w") as output:
        for line in input:
            if line.startswith('>'):
                num = re.findall('[0-9]+', line)[1]
                if int(num) < threshold:
                    break
                else:
                    output.write(line)
            else:
                output.write(line)

# Count the number of scaffolds before and after trimming
pre_trim_num = 0
post_trim_num = 0
with open(input_file, "r") as input:
    for line in input:
        if line.startswith('>'):
            pre_trim_num += 1
with open(output_file, "r") as input:
    for line in input:
        if line.startswith('>'):
            post_trim_num += 1
print("Number of scaffolds before trimming:", pre_trim_num)
print("Number of scaffolds after trimming:", post_trim_num)
print("Number of scaffolds removed:", pre_trim_num - post_trim_num)