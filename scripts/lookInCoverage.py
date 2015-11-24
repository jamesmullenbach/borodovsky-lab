import argparse


parser = argparse.ArgumentParser(description="parser for coverage file using GFF or PTT files")
parser.add_argument('COVERAGE_FILE', help="path to the coverage file you want to parse")
args = vars(parser.parse_args())

coverage_name = args['COVERAGE_FILE']
BIG_COVERAGE_FILE = open(coverage_name, 'r')
# don't read first line into the sequence
description = BIG_COVERAGE_FILE.readline()
coverage_list = BIG_COVERAGE_FILE.read().splitlines()
for line in coverage_list:
    vals = line.split('\t')
    if (int(vals[1]) != int(vals[2]) - 1):
        print "found one"
