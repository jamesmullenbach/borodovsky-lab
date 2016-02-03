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
NUM_GENES = 5
gene_lines_rand = []
while len(gene_lines_rand) < NUM_GENES:
    index = random.randint(0, len(gene_lines) - 1)
    if index == len(gene_lines) - 1:
        start2 = sys.maxint
    else:
        line1 = geneTools.getLineData(gene_lines[index], gene_filetype)
        line2 = geneTools.getLineData(gene_lines[index + 1], gene_filetype)
        start1, start2 = line1[0], line2[0]
    if start1 < start2 - 50 and line1 not in gene_lines_rand:
        gene_lines_rand.append(line1)
gene_lines_rand.sort()

pos_lines, neg_lines = geneTools.readCoverageLines(args['POSITIVE_COVERAGE_FILE'])

ind_pos = 0
ind_neg = 0

x_pos = int(pos_lines[ind_pos].split('\t')[1])
x_neg = int(neg_lines[ind_neg].split('\t')[1])

#ribosome profiling reads are performed by recording OFFSET nt upstream of the 3' end
OFFSET = 30
bothDirections = True

#array to hold values at k nt upstream of stop codon, OFFSET < k  < 0
nt_data = [0 for _ in range(OFFSET * 2)] if bothDirections else [0 for _ in range(OFFSET)]
for line in gene_lines_rand:
    start, end, strand = line[0], line[1], line[2]
    print "start, end, strand: " + str(start) + ", " + str(end) + ", " + strand
    coverage = 0
    if strand == "+":
        while x_pos <= end - OFFSET and ind_pos < len(pos_lines) - 1:
            ind_pos += 1
            x_pos = int(pos_lines[ind_pos].split('\t')[1])
        stop = end + OFFSET if bothDirections else end
        while x_pos < stop and ind_pos < len(pos_lines) - 1:
            score = int(pos_lines[ind_pos].split('\t')[3])
            x_pos = int(pos_lines[ind_pos].split('\t')[1])
            ind_pos += 1
            dist_from_stop = end - x_pos
            #print "dist from stop: " + str(dist_from_stop)
            if abs(dist_from_stop) < OFFSET:
                if bothDirections:
                    index = -1 * dist_from_stop + OFFSET
                else:
                    index = dist_from_stop
                print "location: " + str(x_pos) + ", score: " + str(score) + ", index: " + str(index)
                nt_data[index] += score
    elif strand == "-":
        #looking for OFFSET nt upstream from the 3' end, which is actually OFFSET after the start as defined by GFF, PTT
        begin = start - OFFSET if bothDirections else start
        while x_neg < begin and ind_neg < len(neg_lines) - 1: 
            ind_neg += 1
            x_neg = int(neg_lines[ind_neg].split('\t')[1])
        while x_neg < start + OFFSET and ind_neg < len(neg_lines) - 1:
            score = int(neg_lines[ind_neg].split('\t')[3])
            x_neg = int(neg_lines[ind_neg].split('\t')[1]) 
            ind_neg += 1
            dist_from_stop = x_neg - start
            if abs(dist_from_stop) < OFFSET:
                if bothDirections:
                    index = -1 * dist_from_stop + OFFSET
                else:
                    index = dist_from_stop
                print "location: " + str(x_neg) + ", score: " + str(score) + ", index: " + str(index)
                nt_data[index] += score

if not bothDirections:
    nt_data.reverse()

xaxis = range(-1 * OFFSET, OFFSET) if bothDirections else range(-1 * OFFSET, 0)

#PLOTTING
plt.figure()
#plt.plot(noiseData, 'ro')
plt.plot(xaxis, nt_data, 'ro')
plt.grid(True)
#plt.hist(noiseData, 100)
plt.xlabel('position relative to stop codon')
plt.ylabel('score across ' + str(NUM_GENES) + ' randomly selected genes')
#plt.xlabel('gene index')
#plt.ylabel('coverage value')
plt.title("Coverage at positions within " + str(OFFSET) + "nt from 3' ends in " + str(NUM_GENES) + " randomly selected genes\nwith no gene starting within 50 nt downstream")
plt.show() 
