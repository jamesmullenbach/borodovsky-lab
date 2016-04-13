import argparse
import numpy

#10^-(this/10) is the probability of incorrect mapping position
MIN_QUALITY_ACCEPTABLE = 10
    
SEQUENCE_LENGTH = 4641652
    
pos_ribo_data = {}
neg_ribo_data = {}
for i in range(SEQUENCE_LENGTH):
    pos_ribo_data[i] = 0
    neg_ribo_data[i] = 0
print "done creating zeros dictionary"

for i in range(1734430, 1734445):
    print "#" * 10 + " SRR " + str(i) + " " + "#" * 10

    POSITIVE_OUT_NAME = "positive.ribo"
    NEGATIVE_OUT_NAME = "negative.ribo"
    
    POS_FILE = open(POSITIVE_OUT_NAME, 'w')
    NEG_FILE = open(NEGATIVE_OUT_NAME, 'w')
   
    SAM_FILE = "SRR" + str(i) + ".sam" 
    with open(SAM_FILE, 'r') as f:
        if i != 1734430:
            f.readline()
            f.readline()
            f.readline()
    
        j = 0
        for line in f:
            vals = line.split('\t')
            flag = int(vals[1])
            start = int(vals[3])
            quality = int(vals[4])
            if quality > MIN_QUALITY_ACCEPTABLE:
                #if opposite strand flag not set, it's positive strand
                strand = "+" if flag & 0x10 == 0 else "-"
                if strand == "+":
                    #end is the start position + the length of the sequence for direct strand
                    end = min(start + len(vals[9]), SEQUENCE_LENGTH - 1)
                    pos_ribo_data[end] += 1
                elif strand == "-":
                    #end is the start position - the length of the sequence for complementary strand
                    end = max(0, start - len(vals[9]))
                    neg_ribo_data[end] += 1
            if j % 1000000 == 0 and j > 0:
                print "iterations: " + str(j)
            j += 1
        
for start, score in pos_ribo_data.iteritems():
    if score > 0:
        vals = ["chr", str(start), str(start + 1), str(score)]
        POS_FILE.write('\t'.join(vals) + '\n')
for start, score in neg_ribo_data.iteritems():
    if score > 0:
        vals = ["chr", str(start), str(start + 1), str(score)]
        NEG_FILE.write('\t'.join(vals) + '\n')
