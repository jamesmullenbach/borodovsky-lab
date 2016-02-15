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
NUM_GENES = 100
DIRECT_ONLY = True
SAVE_GENES = False
SEQUENCE_FILENAME = "/storage/james/data/ecoli/NC_000913.fa"

#enforce that no gene on either strand starts within ZONE_LENGTH nt of the stop codon
gene_lines_rand = []
ZONE_LENGTH = 50
attempts = 0
index = 0
while len(gene_lines_rand) < NUM_GENES:
    attempts += 1
    index = random.randint(0, len(gene_lines) - 1)
    line1 = geneTools.getLineData(gene_lines[index], gene_filetype)
    if line1[2] == "-":
        if DIRECT_ONLY:
            #do nothing if only looking at direct genes
            continue
        if index == len(gene_lines) - 1:
            #always add the last gene
            end2 = -1 * sysmaxint
        else:
            line2 = geneTools.getLineData(gene_lines[index + 1], gene_filetype)
            start2 = line2[0]
        end1 = line1[1]
        if end1 < start2 - ZONE_LENGTH and line1 not in gene_lines_rand:
            gene_lines_rand.append(line1)
    elif line1[2] == "+":
        if index == 0:
            #always add the first gene
            start2 = sys.maxint
        else:
            line2 = geneTools.getLineData(gene_lines[index - 1], gene_filetype)
            end2 = line2[1]
        start1 = line1[0]
        if start1 > end2 + ZONE_LENGTH and line1 not in gene_lines_rand:
            gene_lines_rand.append(line1)

print "attempts: " + str(attempts)
gene_lines_rand.sort()
print "num genes: " + str(len(gene_lines_rand))

pos_lines, neg_lines = geneTools.readCoverageLines(args['POSITIVE_COVERAGE_FILE'])

ind_pos = 0
ind_neg = 0

x_pos = int(pos_lines[ind_pos].split('\t')[1])
x_neg = int(neg_lines[ind_neg].split('\t')[1])

#ribosome profiling reads are performed by recording OFFSET nt upstream of the 3' end
OFFSET = 30
bothDirections = True

dirs = args['POSITIVE_COVERAGE_FILE'].split('/')
experiment = dirs[5]
dataset = dirs[6]
RESULTS_FILE = open('/home/james/borodovsky-lab/results/weirdGenes_starts_' + experiment + '_' + dataset + '.txt', 'w')

#array to hold values at k nt up/downstream of start codon, OFFSET < k  < 0
nt_data = [0 for _ in range(OFFSET * 2 + 1)] if bothDirections else [0 for _ in range(OFFSET + 1)]
#iterate over randomly selected genes
for line in gene_lines_rand:
    #get data
    start, end, strand = line[0], line[1], line[2]
    #print "start, end, strand: " + str(start) + ", " + str(end) + ", " + strand
    in_gene = 0
    out_gene = 0
    if strand == "+":
        #move thru coverage file until within 30nt of gene stop codon
        while x_pos <= start - OFFSET - 1 and ind_pos < len(pos_lines) - 1:
            ind_pos += 1
            x_pos = int(pos_lines[ind_pos].split('\t')[1])
        #stop at stop codon + 30nt if going in both directions
        stop = start + OFFSET if bothDirections else start
        while x_pos < stop and ind_pos < len(pos_lines) - 1:
            score = int(pos_lines[ind_pos].split('\t')[3])
            x_pos = int(pos_lines[ind_pos].split('\t')[1])
            ind_pos += 1
            dist_from_start = start - x_pos
            #print "dist from stop: " + str(dist_from_stop)
            if abs(dist_from_start) <= OFFSET:
                if bothDirections:
                    index = -1 * dist_from_start + OFFSET
                else:
                    index = dist_from_start
                #print "location: " + str(x_pos) + ", score: " + str(score) + ", index: " + str(index)
                nt_data[index] += score
                if index >= OFFSET:
                    in_gene += score
                else:
                    out_gene += score
                
    elif strand == "-":
        #looking for OFFSET nt upstream from the 3' end, which is actually OFFSET after the start as defined by GFF, PTT
        begin = end - OFFSET if bothDirections else end
        while x_neg < begin and ind_neg < len(neg_lines) - 1: 
            ind_neg += 1
            x_neg = int(neg_lines[ind_neg].split('\t')[1])
        while x_neg < end + OFFSET + 1 and ind_neg < len(neg_lines) - 1:
            score = int(neg_lines[ind_neg].split('\t')[3])
            x_neg = int(neg_lines[ind_neg].split('\t')[1]) 
            ind_neg += 1
            dist_from_start = x_neg - end
            if abs(dist_from_start) < OFFSET:
                if bothDirections:
                    index = -1 * dist_from_start + OFFSET
                else:
                    index = dist_from_start
                #print "location: " + str(x_neg) + ", score: " + str(score) + ", index: " + str(index)
                nt_data[index] += score
                #print "adding score from position " + str(x_neg) 
                if index >= OFFSET:
                    in_gene += score
                else:
                    out_gene += score
    if in_gene < out_gene:
        #convention that genes are read left to right, so left side is out of the gene when examining starts
        subseq = geneTools.getSubsequence(SEQUENCE_FILENAME, start-30, start+2, strand)
        RESULTS_FILE.write('\t'.join([str(start), str(end), strand, str(out_gene), str(in_gene), subseq + "..."]) + '\n')

in_gene_avg = 0
out_gene_avg = 0
for i in range(len(nt_data)):
    if i <= OFFSET - 1:
        out_gene_avg += nt_data[i]
    else:
        in_gene_avg += nt_data[i]

if in_gene_avg > out_gene_avg:
    print "gene coverage higher"
else:
    print "gene coverage lower"
    SAVE_GENES = True

if SAVE_GENES: 
    for line in gene_lines_rand:
        start, end, strand = line[0], line[1], line[2]
        #RESULTS_FILE.write('\t'.join([str(start), str(end), strand, str(in_gene), str(out_gene)]) + '\n')

#plt.figure()
#plt.plot([-1, 1], [out_gene_avg, in_gene_avg], 'ro')
#plt.show()

xaxis = range(-1 * OFFSET, OFFSET + 1) if bothDirections else range(-1 * OFFSET, 1)

#PLOTTING
plt.figure()
#plt.plot(noiseData, 'ro')
plt.plot(xaxis, nt_data, 'ro')
plt.grid(True)
#plt.hist(noiseData, 100)
plt.xlabel('position relative to start codon')
plt.ylabel('score across ' + str(NUM_GENES) + ' randomly selected genes')
#plt.xlabel('gene index')
#plt.ylabel('coverage value')
plt.title("Coverage of " + str(NUM_GENES) + " randomly selected genes\n with no gene starting within " + str(ZONE_LENGTH) + " nt upstream")
plt.show() 
