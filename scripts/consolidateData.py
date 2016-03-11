import argparse
import numpy

parser = argparse.ArgumentParser(description="Parser for sam data that will be converted to profile data")
parser.add_argument("SAM_FILE", help=".sam file you want to convert to profile format")

args = vars(parser.parse_args())

#10^-(this/10) is the probability of incorrect mapping position
MIN_QUALITY_ACCEPTABLE = 10

POSITIVE_OUT_NAME = args['SAM_FILE'].split('.')[0] + "_positive.ribo"
NEGATIVE_OUT_NAME = args['SAM_FILE'].split('.')[0] + "_negative.ribo"

POS_FILE = open(POSITIVE_OUT_NAME, 'w')
NEG_FILE = open(NEGATIVE_OUT_NAME, 'w')

with open(args['SAM_FILE'], 'r') as f:
    header = f.readline()
    sequence_line = f.readline()
    program_line = f.readline()

    #parse reference sequence length from .sam file
    SEQUENCE_LENGTH = int(sequence_line.split('\t')[2].split(':')[1])
    print "sequence length: " + str(SEQUENCE_LENGTH)

    pos_ribo_data = {}
    neg_ribo_data = {}
    for i in range(SEQUENCE_LENGTH):
        pos_ribo_data[i] = 0
        neg_ribo_data[i] = 0
    print "done creating zeros dictionary"

    i = 0
    for line in f:
        vals = line.split('\t')
        flag = int(vals[1])
        start = int(vals[3])
        quality = int(vals[4])
        if quality > MIN_QUALITY_ACCEPTABLE:
            #if opposite strand flag not set, it's positive strand
            strand = "+" if flag & 0x10 == 0 else "-"
            if strand == "+":
                pos_ribo_data[start] += 1
            elif strand == "-":
                neg_ribo_data[start] += 1
        if i % 100000 == 0 and i > 0:
            print "iterations: " + str(i)
        i += 1
    
for start, score in pos_ribo_data.iteritems():
    if score > 0:
        vals = ["chr", str(start), str(start + 1), str(score)]
        POS_FILE.write('\t'.join(vals) + '\n')
for start, score in neg_ribo_data.iteritems():
    if score > 0:
        vals = ["chr", str(start), str(start + 1), str(score)]
        NEG_FILE.write('\t'.join(vals) + '\n')
