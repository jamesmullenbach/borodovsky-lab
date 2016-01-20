import argparse

#script to get the coverage of a set number of annotated genes as a proportion of total coverage
#for a given experiment. to be run after addStrandData.py is run for an experimental data folder
#and after countTotalCoverages.py is run for each folder as well
#super ad hoc right now, requires both to be run (Unique_positive_combined.bed, negative, and totals.txt must be present)

#set up argument parser
parser = argparse.ArgumentParser(description="parser for genes file and experimental data directory")
parser.add_argument('GENES_FILE', help="path to the file containing annotated genes")
parser.add_argument('DIRECTORY', help="path to the directory containing combined coverage data files")

args = vars(parser.parse_args())

gene_filename = args['GENES_FILE']
gene_filetype = gene_filename.split('.')[1]
if gene_filetype != 'gff' and gene_filetype != 'ptt':
    print "Acceptable file types for GENES_FILE are: .gff, .ptt"
    sys.exit()
GENES_FILE = open(gene_filename, 'r')
if gene_filetype == 'ptt':
    #read off first 3 lines
    GENES_FILE.readline()
    GENES_FILE.readline()
    GENES_FILE.readline() 

#number of genes to consider
NUM_GENES = 200

#open files
POS_STRAND_FILE = open(args['DIRECTORY'] + "Unique_positive_combined.bed", 'r')
NEG_STRAND_FILE = open(args['DIRECTORY'] + "Unique_negative_combined.bed", 'r')
TOTALS_FILE = open(args['DIRECTORY'] + "totals.txt", 'r')
OUT_FILE = open(args['DIRECTORY'] + "geneCoverageProportions.txt", 'w')


#read out totals into variables
pos_total = int(TOTALS_FILE.readline().split(' ')[2])
neg_total = int(TOTALS_FILE.readline().split(' ')[2])

#read coverage lines into variables
pos_lines = POS_STRAND_FILE.read().splitlines()
neg_lines = NEG_STRAND_FILE.read().splitlines()

ind_pos = 0
ind_neg = 0
x_pos = int(pos_lines[ind_pos].split('\t')[1])
x_neg = int(neg_lines[ind_neg].split('\t')[1])

for i in range(NUM_GENES):
    gene = GENES_FILE.readline()
    values = gene.split('\t')
    gene_score = 0
    gene_proportion = 0
    #subtract 1 since python indexes from 0 and GFF/PTT index from 1
    if gene_filetype == "gff":
        start = int(values[3]) - 1
        end = int(values[4])
        strand = values[6]
    elif gene_filetype == "ptt":
        start = int(values[0].split('..')[0]) - 1
        end = int(values[0].split('..')[1])
        strand = values[1]
    
    #get coverage data for this gene from given combined coverage file
    if strand == "+":
        #get data from positive strand file
        #first go to first positive line that is greater than or equal to start
        while (x_pos < start):
            ind_pos += 1
            x_pos = int(pos_lines[ind_pos].split('\t')[1])
        #add the end of this loop, x_pos is at or past the start
        while (x_pos <= end):
            vals = pos_lines[ind_pos].split('\t')
            score = int(vals[3])
            x_pos = int(vals[1])
            gene_score += score
            ind_pos += 1
        gene_proportion = float(gene_score) / pos_total
    elif strand == "-":
        #get data from negative strand file
        #first go to first negative line that is greater than or equal to start
        while (x_neg < start):
            ind_neg += 1
            x_neg = int(neg_lines[ind_neg].split('\t')[1])
        #add the end of this loop, x_pos is at or past the start
        while (x_neg <= end):
            vals = neg_lines[ind_neg].split('\t')
            score = int(vals[3])
            x_neg = int(vals[1])
            gene_score += score
            ind_neg += 1
        gene_proportion = float(gene_score) / neg_total

    #write to an output file in the same directory
    line = '\t'.join([str(start), str(end), strand, str(gene_proportion)]) + "\n"
    OUT_FILE.write(line)
