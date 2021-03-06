import geneTools
import argparse
import matplotlib.pyplot as plt
import random
import sys

#look for the amount of noise by checking coverage values at OFFSETnt upstream from gene stops
parser = argparse.ArgumentParser(description="parser for a gene and a coverage file")
parser.add_argument("GENES_FILE", help="path to gene file you want to use")
parser.add_argument("POSITIVE_COVERAGE_FILE", help="path to coverage file you want to use")

args = vars(parser.parse_args())

gene_lines, gene_filetype = geneTools.readORFLines(args['GENES_FILE'])

#select a certain number of random genes
NUM_GENES = 500
gene_lines_rand = []
for _ in range(NUM_GENES):
    line_data = geneTools.getLineData(gene_lines[random.randint(0, len(gene_lines) - 1)], gene_filetype)
    if line_data not in gene_lines_rand:
        gene_lines_rand.append(line_data)
gene_lines_rand.sort()

pos_lines, neg_lines = geneTools.readCoverageLines(args['POSITIVE_COVERAGE_FILE'])

ind_pos = 0
ind_neg = 0

x_pos = int(pos_lines[ind_pos].split('\t')[1])
x_neg = int(neg_lines[ind_neg].split('\t')[1])

#ribosome profiling reads are performed by recording OFFSET nt upstream of the 3' end
OFFSET = 30

#array to hold values at k nt upstream of stop codon, OFFSET < k  < 0
nt_data = [0 for _ in range(OFFSET + 1)]
for line in gene_lines_rand:
    start, end, strand = line[0], line[1], line[2]
    coverage = 0
    if strand == "+":
        while x_pos <= end - OFFSET - 1 and ind_pos < len(pos_lines) - 1:
            ind_pos += 1
            x_pos = int(pos_lines[ind_pos].split('\t')[1])
        while x_pos < end and ind_pos < len(pos_lines) - 1:
            score = int(pos_lines[ind_pos].split('\t')[3])
            x_pos = int(pos_lines[ind_pos].split('\t')[1])
            ind_pos += 1
            dist_from_stop = end - x_pos
            if dist_from_stop <= OFFSET and dist_from_stop >= 0:
                nt_data[dist_from_stop] += score
    elif strand == "-":
        #looking for OFFSET nt upstream from the 3' end, which is actually OFFSET after the start as defined by GFF, PTT
        while x_neg < start and ind_neg < len(neg_lines) - 1: 
            ind_neg += 1
            x_neg = int(neg_lines[ind_neg].split('\t')[1])
        while x_neg < start + OFFSET + 1 and ind_neg < len(neg_lines) - 1:
            score = int(neg_lines[ind_neg].split('\t')[3])
            x_neg = int(neg_lines[ind_neg].split('\t')[1]) 
            ind_neg += 1
            dist_from_stop = x_neg - start
            if (dist_from_stop <= OFFSET and dist_from_stop >= 0):
                nt_data[dist_from_stop] += score

#reverse data to get more readable order
nt_data.reverse()
print "nt data: " + str(nt_data)

#PLOTTING
plt.figure()
#plt.plot(noiseData, 'ro')
plt.plot(range(-1 * OFFSET, 1), nt_data, 'ro')
plt.grid(True)
#plt.hist(noiseData, 100)
plt.xlabel('nt upstream from stop codon')
plt.ylabel('score across ' + str(NUM_GENES) + ' randomly selected genes')
#plt.xlabel('gene index')
#plt.ylabel('coverage value')
plt.title("Coverage at positions within " + str(OFFSET) + "nt upstream from 3' ends in " + str(NUM_GENES) + " randomly selected genes")
plt.show() 
