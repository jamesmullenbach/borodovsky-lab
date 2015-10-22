import sys
import numpy
import matplotlib.pyplot as plt
import pylab
import argparse

#take in GFF file and coverage file
#for each gene in GFF file, search for ranges in coverage file
#over each range, get average and variance, accounting for 0's which are omitted
#need to figure out the strand the ribosomes were on

#if (len(sys.argv) < 3):
#    print "Usage: python " + sys.argv[0] + " COVERAGE_FILE GENES_FILE OUTPUT_FILE_NAME"
#    print "Acceptable file types for GENES_FILE are: .gff, .ptt"
#    sys.exit()

#gene_file_name = sys.argv[2]

parser = argparse.ArgumentParser(description="parser for coverage file using GFF or PTT files")
parser.add_argument('COVERAGE_FILE', help="path to the coverage file you want to parse")
parser.add_argument('GENES_FILE', help="path to the gene file (GFF or PTT) you want to use for genes")
parser.add_argument('OUTPUT_FILE_NAME', help="name of output file to create")
parser.add_argument('--noPeaks', help="use this option to set a threshold value for coverage", type=int)
parser.add_argument('--peaksOnly', help="use this option to consider only peaks in genes")
args = vars(parser.parse_args())

THRESHOLD_COVERAGE = sys.maxint
if args['noPeaks'] is not None:
    THRESHOLD_COVERAGE = args['noPeaks']
    print "running for value " + str(THRESHOLD_COVERAGE)
if args['peaksOnly'] is not None:
    peaksOnly = True

coverage_name = args['COVERAGE_FILE']
BIG_COVERAGE_FILE = open(coverage_name, 'r')
# don't read first line into the sequence
description = BIG_COVERAGE_FILE.readline()
coverage_list = BIG_COVERAGE_FILE.read().splitlines()

gene_file_name = args['GENES_FILE']
gene_filetype = gene_file_name.split('.')[1]
if gene_filetype != 'gff' and gene_filetype != 'ptt':
    print "Acceptable file types for GENES_FILE are: .gff, .ptt"
    sys.exit()
print "file type: " + gene_filetype
GENES_FILE = open(gene_file_name, 'rb')
if gene_filetype == 'ptt':
    #read off the three header lines in .ptt file
    GENES_FILE.readline()
    GENES_FILE.readline()
    GENES_FILE.readline()

GENE_LINES = GENES_FILE.read().splitlines()

out_name = args['OUTPUT_FILE_NAME']
OUT_FILE = open(out_name, 'w')
coverage_line = 0
gene_scores = []
gene_lengths = []
gene_variances = []
gene_stds = []

######BEGIN LOOP THROUGH GFF FILE######
for line in GENE_LINES:
    values = line.split('\t')
    if (len(values) < 5 and gene_filetype == 'gff'):
        print "reached non-uniform line in GFF file: " + str(line)
        break
    #parse interval values and gene ids from gff lines
    start = 0
    end = 0
    gene_id = ""
    if gene_filetype == "gff":
        start = int(values[3])
        end = int(values[4])
        details = values[len(values) - 1]
        gene_id = details.split(' ')[0]
    elif gene_filetype == "ptt":
	start = int(values[0].split('..')[0])
        end = int(values[0].split('..')[1])
        #print "start: " + str(start) + " end: " + str(end)
    gene_length = end - start + 1
    #go through coverage file until you find an index >= start
    if (coverage_line >= len(coverage_list)):
        print "end of coverage file reached"
        break

    ##########AVERAGE CALCULATION#################
    line_data = coverage_list[coverage_line].split('\t')
    base = int(line_data[1])
    gene_score = 0
    starting_coverage_line = coverage_line
    if (base >= start):
        #add value to running total
        coverage_value = min(float(line_data[3]), THRESHOLD_COVERAGE)
        gene_score += coverage_value
        #print "adding value to score: " + str(line_data[3])
    while (base <= end):
	coverage_line += 1
        if (coverage_line >= len(coverage_list)):
            #we have reached the end of the coverage file
            #last gene has some locations without scores
            #break while
            #print "blahbalh hax"
            break
        line_data = coverage_list[coverage_line].split('\t')
        #print "line_data: " + str(line_data)
        base = int(line_data[1])
        #print "base: " + str(base) + ", start: " + str(start)
        if (base >= start):
            #add value to running total
            coverage_value = min(float(line_data[3]), THRESHOLD_COVERAGE)
            gene_score += coverage_value
           #print "adding value to score: " + str(line_data[3])
    gene_average = gene_score / gene_length
    #print "gene average: " + str(gene_average)


    #########VARIANCE CALCULATION###############
    #go back through the array to get the variance
    gene_variance = 0
    coverage_line = starting_coverage_line
    line_data = coverage_list[coverage_line].split('\t')
    base = int(line_data[1])
    score = float(line_data[3])
    last_base = base
    if (base > start):
        #handle edge case in which first nonzero-coverage entry is after start
        gene_variance += ((gene_average) ** 2) * (base - start)
        gene_variance += (gene_average - score) ** 2
    while (base <= end):
        coverage_line += 1
        if (coverage_line >= len(coverage_list)):
            #break while
            break
        line_data = coverage_list[coverage_line].split('\t')
        base = int(line_data[1])
        if (base > start):
            #add value to running total
            if (base - last_base > 1):
                #there was a gap between consecutive coverage lines
                #so add the variance from a value of 0 for each spot in the gap
                #only add starting from the gene
                num_zeroes = 0
                if (last_base < start):
                    num_zeroes = base - start
                else:
                    num_zeroes = base - last_base - 1
                gene_variance += ((gene_average) ** 2) * num_zeroes
            score = float(line_data[3])
            gene_variance += (gene_average - score) ** 2
        last_base = base
    #now the last base is the final base that has a nonzero coverage score
    if (last_base < end):
        gene_variance += ((gene_average) **2) * (end - last_base)
    gene_variance = gene_variance / gene_length
  
    if gene_average < 4000:
        #store values
        gene_lengths.append(gene_length)
        gene_scores.append(gene_average)
        gene_variances.append(gene_variance)
        gene_stds.append(gene_variance ** 0.5)

    #write values to file
    out_line = ">" + gene_id + "|start=" + str(start) + "|end=" + str(end)
    out_line += "|strand=POSITIVE\n"
    out_line += "total score: " + str(gene_score) + '\n'
    out_line += "average score: " + str(gene_average) + '\n'
    out_line += "score variance: " + str(gene_variance) + '\n'
    OUT_FILE.write(out_line)



##############PLOTTING#############
#plot average scores with std error bars
#plt.figure()
plt.plot(gene_lengths, gene_scores, 'ro')
#plt.errorbar(gene_lengths, gene_scores, gene_stds, fmt='bo')
plt.xlabel('Gene length')
plt.ylabel('Gene average coverage score')
image_name = "/home/james/plots/ceilingvalues/" + str(THRESHOLD_COVERAGE) + ".png"
plt.savefig(image_name, bbox_inches='tight')
#plt.show()

#plot gene variances versus length
#plt.figure()
#plt.plot(gene_lengths, gene_variances, 'ro')
#plt.xlabel('Gene length')
#plt.ylabel('Gene score variance')
#plt.show()

#plot gene std deviations versus length
#plt.figure()
#plt.plot(gene_lengths, gene_stds, 'ro')
#plt.xlabel('Gene length')
#plt.ylabel('Gene score standard deviation')
#plt.show()

#plot histogram of gene average score
plt.figure()
plt.hist(gene_scores, 50)
plt.xlabel("Gene average coverage score")
plt.ylabel("Frequency")
image_name = "/home/james/plots/ceilingvalues/" + str(THRESHOLD_COVERAGE) + "hist.png"
plt.savefig(image_name, bbox_inches="tight")
#plt.show()

#plot 2d histogram of gene length and average score
#plt.figure()
#plt.hist2d(gene_lengths, gene_scores, bins=10)
#plt.xlabel("Gene length")
#plt.ylabel("Gene average coverage score")
#plt.show()
