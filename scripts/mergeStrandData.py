import sys
import argparse

#script to combine coverage data files, one for each strand, into one file, keeping strand information

#create argument parser
parser = argparse.ArgumentParser(description="parser for two coverage file inputs")
parser.add_argument('POSITIVE_COVERAGE_FILE', help="path to the positive coverage file you want to parse")
parser.add_argument('NEGATIVE_COVERAGE_FILE', help="path to the negative coverage file you want to parse")
parser.add_argument('OUT_FILE', help="path to the output file you want to create")

args = vars(parser.parse_args())

positive_coverage_name = args['POSITIVE_COVERAGE_FILE']
POS_FILE = open(positive_coverage_name, 'r')
pos_list = POS_FILE.read().splitlines()

negative_coverage_name = args['NEGATIVE_COVERAGE_FILE']
NEG_FILE = open(negative_coverage_name, 'r')
neg_list = NEG_FILE.read().splitlines()

out_name = args['OUT_FILE']
OUT_FILE = open(out_name, 'w')

x_pos = 0
x_neg = 0

while (x_pos < len(pos_list) and x_neg < len(neg_list)):
    pos_line = pos_list[x_pos]
    neg_line = neg_list[x_neg]
    
    pos_base = pos_line[1]
    pos_score = pos_line[3]
    neg_base = neg_line[1]
    neg_score = neg_line[3]

    if (pos_base <= neg_base):
        #write pos_line to file with strand information
        out_line = pos_line + "\t" + "+"
        x_pos += 1
    else:
        #write neg_line to file with strand information
        out_line = neg_line + "\t" + "-"
        x_neg += 1
        
    OUT_FILE.write(out_line)


if (x_pos < len(pos_list)):
    #add the rest of the positive coverage
    while x_pos < len(pos_list):
        out_line = pos_line + "\t" + "+"
        x_pos += 1
        OUT_FILE.write(out_line)
elif (x_neg < len(neg_list)):
    #add the rest of the negative coverage
    while x_neg < len(neg_list):
        out_line = pos_line + "\t" + "-"
        x_neg += 1
        OUT_FILE.write(out_line)
