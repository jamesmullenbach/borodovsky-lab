import sys
import numpy
import matplotlib.pyplot as plt
import pylab
import argparse

#parse a coverage file and classify each gene as upward or downward trending
#take in GFF file and coverage file
#for each gene in GFF file, search for ranges in coverage file
#TODO: incorporate strand information

#set up argument parser
parser = argparse.ArgumentParser(description="parser for coverage file that classifies genes as upward or downward trending")
parser.add_argument('COVERAGE_FILE', help="path to the coverage file you want to parse")
parser.add_argument('GENES_FILE', help="path to the gene file (GFF or PTT) you want to use for genes")
parser.add_argument('--out', help="name of output file to create")
args = vars(parser.parse_args())

if args['out'] is not None:
    out_name = args['out']
    OUT_FILE = open(out_name, 'w')

OUTLIER_THRESHOLD = 4000

coverage_name = args['COVERAGE_FILE']
BIG_COVERAGE_FILE = open(coverage_name, 'r')
# don't read first line into the sequence
description = BIG_COVERAGE_FILE.readline()
coverage_list = BIG_COVERAGE_FILE.read().splitlines()

#enforce correct file extension
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

coverage_line = 0
gene_scores = []
gene_lengths = []
gene_upwards = []
gene_upwards_lengths = []
gene_downwards = []
gene_downwards_lengths = []

######BEGIN LOOP THROUGH GENE FILE######
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
        strand = values[6] 
        gene_id = details.split(' ')[0]
    elif gene_filetype == "ptt":
        start = int(values[0].split('..')[0])
        end = int(values[0].split('..')[1])
        strand = values[1]

    gene_length = end - start + 1
    #go through coverage file until you find an index >= start
    if (coverage_line >= len(coverage_list)):
        print "end of coverage file reached"
        break

    ##########AVERAGE CALCULATION#################
    line_data = coverage_list[coverage_line].split('\t')
    base = int(line_data[1])
    gene_score = 0
    gene_left_score = 0
    gene_right_score = 0
    starting_coverage_line = coverage_line
    if (base >= start):
        #add value to running total
        coverage_value = float(line_data[3])
        gene_score += coverage_value
        if (base < (start + gene_length / 2)):
            gene_left_score += coverage_value
        else:
            if (base <= end):
                gene_right_score += coverage_value
    while (base <= end):
        coverage_line += 1
        if (coverage_line >= len(coverage_list)):
            #we have reached the end of the coverage file
            #last gene has some locations without scores
            #break while
            break
        line_data = coverage_list[coverage_line].split('\t')
        base = int(line_data[1])
        if (base >= start):
            #add value to running total
            coverage_value = float(line_data[3])
            gene_score += coverage_value
            #add up left and right half scores of each gene
            if (base < (start + gene_length / 2)):
                gene_left_score += coverage_value
            else:
                if (base <= end):
                    gene_right_score += coverage_value
    gene_average = gene_score / gene_length
    trend = ""

    #ignore outliers
    if gene_average < OUTLIER_THRESHOLD:
        if gene_right_score > gene_left_score:
            if strand == "+":
                #upward trending
                gene_upwards.append(gene_average)
                gene_upwards_lengths.append(gene_length)
                trend = "upwards"
            elif strand == "-":
                #downard trending
                gene_downwards.append(gene_average)
                gene_downwards_lengths.append(gene_length)
                trend = "downwards"
        else:
            if strand == "+":
                #downard trending
                gene_downwards.append(gene_average)
                gene_downwards_lengths.append(gene_length)
                trend = "downwards"
            elif strand == "-":
                #upward trending
                gene_upwards.append(gene_average)
                gene_upwards_lengths.append(gene_length)
                trend = "upwards"
    #write values to file
    out_line = ">" + gene_id + "|start=" + str(start) + "|end=" + str(end)
    out_line += "|strand=POSITIVE\n"
    out_line += "average score: " + str(gene_average) + '\n'
    out_line += "gene trend: " + trend + '\n'
    OUT_FILE.write(out_line)

print "Number of upward trending genes: " + str(len(gene_upwards))
print "Number of downward trending genes: " + str(len(gene_downwards))

##############PLOTTING###############
ROOT_PLOT_DIR = "/home/james/borodovsky-lab/plots/trends/"
#plot average scores with std error bars
plt.figure()
plt.plot(gene_upwards_lengths, gene_upwards, 'go')
plt.plot(gene_downwards_lengths, gene_downwards, 'ro')
plt.xlabel('Gene length')
plt.ylabel('Gene average coverage score')
image_name = ROOT_PLOT_DIR + "trend.png"
plt.show()
#plt.savefig(image_name, bbox_inches='tight')

plt.figure()
plt.plot(gene_upwards_lengths, gene_upwards, 'go')
plt.xlabel('Gene length')
plt.ylabel('Gene average coverage score')
plt.title("Upward trending genes")
plt.show()

plt.figure()
plt.plot(gene_downwards_lengths, gene_downwards, 'ro')
plt.xlabel('Gene length')
plt.ylabel('Gene average coverage score')
plt.title("Downward trending genes")
plt.show()

 
