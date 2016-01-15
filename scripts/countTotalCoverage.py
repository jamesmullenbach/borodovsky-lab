import argparse


#count the total coverage for every entry in a coverage file
parser = argparse.ArgumentParser(description="parser for coverage file")
parser.add_argument('COVERAGE_FILE', help="path to the coverage file you want to count")

args = vars(parser.parse_args())

coverage_filename = args['COVERAGE_FILE']
COVERAGE_FILE = open(coverage_filename, 'r')
lines = COVERAGE_FILE.read().splitlines()

total = 0
for line in lines:
    total += int(line.split("\t")[3])

print "total: " + str(total)
