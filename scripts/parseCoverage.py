import sys
import numpy
import matplotlib.pyplot as plt
import pylab
import argparse
from operator import itemgetter

#script takes in GFF file and coverage file
#for each gene/orf in GFF file, search for ranges in coverage file
#over each range, get average and variance, accounting for 0's which are omitted
#TODO: incorporate strand information
#gene file has genes that are on direct or complementary strand
#if we combine the separate files for each strand into one, and keep the strand information
#then we can use the whole GFF file and know the strand information as well, so we won't add coverage from the opposite strand

#create argument parser with lots of options for plotting
parser = argparse.ArgumentParser(description="parser for coverage file using GFF or PTT files")
parser.add_argument('COVERAGE_FILE', help="path to the coverage file you want to parse")
parser.add_argument('GENES_FILE', help="path to the gene file (GFF or PTT) you want to use for genes")
parser.add_argument('--out', help="name of output file to create")
parser.add_argument('--noPeaks', help="use this option to set a threshold value for coverage", type=int)
parser.add_argument('--peaksOnly', help="use this option to consider only peaks in genes", type=int)
parser.add_argument('--orf', help="set flag if dealing with ORFs, not verified genes, for plotting purposes", action="store_true")
parser.add_argument('--plotScore', help="set flag if plot of average coverage score vs. length is desired", action="store_true")
parser.add_argument('--plotSTD', help="set flag if plot of standard deviation vs. length is desired", action="store_true")
parser.add_argument('--hist', help="give number of bins if histogram of average coverage score is desired", type=int)
parser.add_argument('--histRange', help="give histogram lower range", type=int)
parser.add_argument('--top50', help="use this flag to give an output file to write the top 50 gene/ORFs to")

args = vars(parser.parse_args())

ROOT_PLOT_DIR = "/home/james/borodovsky-lab/plots/"
THRESHOLD_COVERAGE = sys.maxint
#determine the gene average score above which to discard a gene because it is an outlier
OUTLIER_THRESHOLD = 1800
PEAK_VALUE = 0
if args['out'] is not None:
    out_name = args['out']
    OUT_FILE = open(out_name, 'w')
if args['noPeaks'] is not None:
    THRESHOLD_COVERAGE = args['noPeaks']
    ROOT_PLOT_DIR += "ceiling-values/"
    print "running for value " + str(THRESHOLD_COVERAGE)
if args['peaksOnly'] is not None:
    PEAK_VALUE = args['peaksOnly']
    print "PEAK_VALUE: " + str(PEAK_VALUE)
    ROOT_PLOT_DIR += "peak-values/"

coverage_name = args['COVERAGE_FILE']
BIG_COVERAGE_FILE = open(coverage_name, 'r')
# don't read first line into the sequence
description = BIG_COVERAGE_FILE.readline()
coverage_list = BIG_COVERAGE_FILE.read().splitlines()

#open gene/orf file and enforce correct file extension
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
gene_variances = []
gene_stds = []
gene_numPeaks = []
gene_data = []
max_score = 0
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
        strand = values[6]
        details = values[len(values) - 1]
        gene_id = details.split(' ')[0]
    elif gene_filetype == "ptt":
        start = int(values[0].split('..')[0])
        end = int(values[0].split('..')[1])
        strand = values[1]
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
            coverage_value = min(float(line_data[3]), THRESHOLD_COVERAGE)
            gene_score += coverage_value
    gene_average = gene_score / gene_length

    if args['peaksOnly'] is not None:
        ##################PEAKS CALCULATION###############
        gene_peaks = 0
        #reset to starting line in coverage file
        coverage_line = starting_coverage_line 
        line_data = coverage_list[coverage_line].split('\t')
        base = int(line_data[1])
        score = float(line_data[3])
        last_base = base
        if (base >= start):
            #add value if it is a peak
            if (score - gene_average) > PEAK_VALUE:
                gene_peaks += 1
        while (base <= end):
            coverage_line += 1
            if (coverage_line >= len(coverage_list)):
                break
            line_data = coverage_list[coverage_line].split('\t')
            base = int(line_data[1])
            if base >= start:
                coverage_value = float(line_data[3])
                if (coverage_value - gene_average) > PEAK_VALUE:
                    gene_peaks += 1

    #########VARIANCE CALCULATION###############
    #go back through the array to get the variance
    gene_variance = 0
    coverage_line = starting_coverage_line
    line_data = coverage_list[coverage_line].split('\t')
    base = int(line_data[1])
    score = min(THRESHOLD_COVERAGE, float(line_data[3]))
    last_base = base
    if (base > start):
        #handle case in which first nonzero-coverage entry is after start
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
            score = min(THRESHOLD_COVERAGE, float(line_data[3]))
            gene_variance += (gene_average - score) ** 2
        last_base = base
    #now the last base is the final base that has a nonzero coverage score
    if (last_base < end):
        gene_variance += ((gene_average) **2) * (end - last_base)
    gene_variance = gene_variance / gene_length
  
    #remove outliers
    if gene_average < OUTLIER_THRESHOLD:
        #store values
        gene_lengths.append(gene_length)
        gene_scores.append(gene_average)
        gene_variances.append(gene_variance)
        gene_stds.append(gene_variance ** 0.5)
        gene_data.append((start, end, strand,gene_length,  gene_average, gene_variance ** 0.5))
        if gene_average > max_score:
            max_score = gene_average
            max_gene = (start, end, strand, gene_average)
        if args['peaksOnly'] is not None:
            gene_numPeaks.append(gene_peaks)

    if args['out'] is not None:
        #write values to file
        strand_str = "DIRECT" if strand == '+' else "COMPLEMENTARY"
        out_line = ">" + gene_id + "|start=" + str(start) + "|end=" + str(end)
        out_line += "|strand=" + strand_str + "\n"
        out_line += "total score: " + str(gene_score) + '\n'
        out_line += "average score: " + str(gene_average) + '\n'
        out_line += "score variance: " + str(gene_variance) + '\n'
        OUT_FILE.write(out_line)

#print top ORF by score
print "max score gene/ORF: " + str((max_gene[0], max_gene[1])) + " score: " + str(max_gene[3]) + " strand: " + str(max_gene[2])

##############PLOTTING#############

whatisit = "Gene"
if args['orf'] is not None:
    whatisit = "ORF"


#####sort and plot top 50 in terms of average coverage score, and save to file
if args['top50'] is not None:
    sorted_data = sorted(gene_data, key=itemgetter(4), reverse=True) #sort by average coverage descending
    print "sorted_data[0]: " + str(sorted_data[0])
    top50_lengths = [sorted_data[i][3] for i in range(50)]
    top50_averages = [sorted_data[i][4] for i in range(50)]
    plt.figure()
    plt.plot(top50_lengths, top50_averages, 'ro')
    plt.xlabel(whatisit + ' length')
    plt.ylabel(whatisit + ' average coverage score')
    plt.show()
    #write values to file
    OUT_FILE = open(args['top50'], 'w')
    for i in range(50):
        start, end, strand, gene_length, gene_average, gene_std = sorted_data[i]
        details = "avg_coverage=" + str(gene_average)
        data = ["NC_000913", "candidate ORF", ".", str(start), str(end), ".", str(strand), ".", details]
        line = "\t".join(data) + "\n"
        OUT_FILE.write(line)
        #strand_str = "POSITIVE" if strand == '+' else "NEGATIVE"
        #out_line = ">" + gene_id + "|start=" + str(start) + "|end=" + str(end)
        #out_line += "|strand=" + strand_str + "\n"
        #out_line += "gene length: " + str(gene_length) + '\n'
        #out_line += "average score: " + str(gene_average) + '\n'
        #out_line += "score standard deviation: " + str(gene_std) + '\n'
        #OUT_FILE.write(out_line)
    sys.exit(0)


#plot average scores with std error bars
if args['plotScore'] is not None:
    print "plotScore arg was not None"
    plt.figure()
    plt.plot(gene_lengths, gene_scores, 'ro')
    #plt.errorbar(gene_lengths, gene_scores, gene_stds, fmt='bo')
    plt.xlabel(whatisit + ' length')
    plt.ylabel(whatisit + ' average coverage score')
    #image_name = "score.png"
    #if args["noPeaks"] is not None:
    #    image_name = ROOT_PLOT_DIR + str(THRESHOLD_COVERAGE) + image_name
    #else:
    #    image_name = ROOT_PLOT_DIR + image_name
    #plt.savefig(image_name, bbox_inches='tight')
    plt.show()

#plot gene std deviations versus length
if args['plotSTD'] is not None:
    plt.figure()
    plt.plot(gene_lengths, gene_stds, 'ro')
    plt.xlabel(whatisit + ' length')
    plt.ylabel(whatisit + ' score standard deviation')
    image_name = "std.png"
    if args["noPeaks"] is not None:
        image_name = ROOT_PLOT_DIR + str(THRESHOLD_COVERAGE) + image_name
    else:
        image_name = ROOT_PLOT_DIR + image_name
    #plt.savefig(image_name, bbox_inches="tight")
    plt.show()

#plot histogram of gene average score
if args['hist'] is not None:
    plt.figure()
    lower = 10
    if args['histRange'] is not None:
        lower = args['histRange']
    plt.hist(gene_scores, args['hist'], range=(lower, max(gene_scores)))
    plt.xlabel(whatisit + " average coverage score")
    plt.ylabel("Frequency")
    image_name = "hist.png"
    if args["noPeaks"] is not None:
        image_name = ROOT_PLOT_DIR + str(THRESHOLD_COVERAGE) + image_name
    else:
        image_name = ROOT_PLOT_DIR + image_name
    #plt.savefig(image_name, bbox_inches="tight")
    titlestr = "Average score histogram with %i bins and minimum score %i" % (args['hist'], lower)
    plt.title(titlestr)
    plt.show()

#plot number of peaks for each gene against length
if len(gene_numPeaks) > 0:
    plt.figure()
    plt.plot(gene_lengths, gene_numPeaks, 'ro')
    plt.xlabel(whatisit + " length")
    ylabelstr = whatisit + " number of peaks for peak value " + str(PEAK_VALUE)
    plt.ylabel(ylabelstr)
    image_name = ROOT_PLOT_DIR + str(PEAK_VALUE) + "peaks.png"
    plt.savefig(image_name, bbox_inches="tight")

    #plot histogram of number of peaks
    plt.figure()
    plt.hist(gene_numPeaks, 50)
    plt.xlabel(ylabelstr)
    plt.ylabel("Frequency")
    image_name = ROOT_PLOT_DIR + str(PEAK_VALUE) + "peakshist.png"
    plt.savefig(image_name, bbox_inches="tight")

