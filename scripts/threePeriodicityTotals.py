#script to look for three periodicity of ribosome profiling data within genes
#within a gene, examine the 1->4->7 bases, 2->5->8 bases, 3->6->9 bases
import sys
import argparse
import numpy
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="parser for ribosome profiling data file and annotated genes file")
parser.add_argument("GENES_FILE", help="path to the annotated genes file you want to parse")
parser.add_argument("POSITIVE_COVERAGE_FILE", help="path to the direct strand ribosome profiling data file you want to parse")

args = vars(parser.parse_args())

#how many genes to consider
NUM_GENES = 100

#open files
negative_filename = args['POSITIVE_COVERAGE_FILE'].replace("positive", "negative")

POSITIVE_COVERAGE_FILE = open(args['POSITIVE_COVERAGE_FILE'], 'r')
pos_lines = POSITIVE_COVERAGE_FILE.read().splitlines()

NEGATIVE_COVERAGE_FILE = open(negative_filename, 'r')
neg_lines = NEGATIVE_COVERAGE_FILE.read().splitlines()

gene_filename = args['GENES_FILE']
gene_filetype = gene_filename.split('.')[1]
if gene_filetype != "gff" and gene_filetype != "ptt":
    print "only GFF or PTT files allowed for GENES_FILE"
    sys.exit(0)
GENES_FILE = open(args['GENES_FILE'], 'r')

#read off leading three lines in ptt format
if gene_filetype == "ptt":
    GENES_FILE.readline()
    GENES_FILE.readline()
    GENES_FILE.readline()

gene_lines = []
for i in range(NUM_GENES):
    gene_lines.append(GENES_FILE.readline())

ind_pos = 0
ind_neg = 0

x_pos = int(pos_lines[ind_pos].split('\t')[1])
x_neg = int(neg_lines[ind_neg].split('\t')[1])


ones = []
twos = []
threes = []
for i in range(NUM_GENES):
    #parse information from the line
    line = gene_lines[i]
    values = line.split('\t')
    one = 0
    two = 0
    three = 0
    if gene_filetype == "gff":
        start = int(values[3])
        end = int(values[4])
        strand = values[6]
    if gene_filetype == "ptt":
        start = int(values[0].split('..')[0])
        end = int(values[0].split('..')[1])
        strand = values[1]
    
    if strand == "+":
        while (x_pos < start):
            ind_pos += 1
            x_pos = int(pos_lines[ind_pos].split('\t')[1])
        while (x_pos < end):
            vals = pos_lines[ind_pos].split('\t')
            score = int(vals[3])
            x_pos = int(vals[1])
            if x_pos % 3 == 0:
                three += score
            elif x_pos % 3 == 1:
                one += score
            elif x_pos % 3 == 2:
                two += score
            ind_pos += 1
    elif strand == "-":
        #this actually goes backward through the gene, but it's okay in this case because it just adds the coverage scores
        while (x_neg < start):
            ind_neg += 1
            x_neg = int(neg_lines[ind_neg].split('\t')[1])
        while (x_neg < end):
            vals = neg_lines[ind_neg].split('\t')
            score = int(vals[3])
            x_neg = int(vals[1])
            if x_neg % 3 == 0:
                three += score
            elif x_neg % 3 == 1:
                one += score
            elif x_neg % 3 == 2:
                two += score
            ind_neg += 1
    
    ones.append(one)
    twos.append(two)
    threes.append(three)

x_axis = range(NUM_GENES)
#PLOTTING
plt.figure()
one_handle, = plt.plot(x_axis, ones, 'ro', label="ones")
two_handle, = plt.plot(x_axis, twos, 'bo', label="twos")
three_handle, = plt.plot(x_axis, threes, 'go', label="threes")
plt.xlabel("Gene index")
plt.ylabel("Coverage total")
plt.title("Coverage for 1st, 2nd, 3rd bases")
#plt.legend([one_handle, two_handle, three_handle])
plt.show()




 
