import sys
import argparse
import os
from operator import itemgetter

#script to add up all the coverage values for either positive or negative strand (unique) within a directory

#create argument parser
parser = argparse.ArgumentParser(description="parser for directory and pos/neg flag")
parser.add_argument('PATH', help="path to the folder containing positive/negative strand coverage data")
parser.add_argument('STRAND', help="set if you want the negative strand. default is positive strand")

args = vars(parser.parse_args())

path = args['PATH']

strand = "negative" if args['STRAND'] == "negative" else "positive"
print "strand: " + strand

print "path: " + path

#total will be a dictionary - total[position] = score
total = {}

#go through each file, adding to the total each time
name_ending = "Unique_" + strand + ".bed"
print "name_ending: " + name_ending
for root, dirs, files in os.walk(path):
    for name in files:
        if name.endswith(name_ending):
            print "found a file of the right name"
            #open file
            full_name = path + name
            COVERAGE_FILE = open(full_name, 'r')
            coverage_lines = COVERAGE_FILE.read().splitlines()
            for line in coverage_lines:
                values = line.split('\t')
                base = int(values[1])
                score = int(values[3])
                try:
                    total[base] += int(score)
                except KeyError:
                    total[base] = int(score)

#open output file
out_file_path = path + "Unique_" + strand + "_combined.bed"
OUT_FILE = open(out_file_path, 'w')


#write total dict to file (sorted)
for base in sorted(total.iterkeys()):
    score = total[base]
    line_values = ['chr', str(base), str(base + 1), str(score)]
    line = "\t".join(line_values) + "\n"
    OUT_FILE.write(line)


