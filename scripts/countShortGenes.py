import sys
import argparse
#script to count number of genes lower than 90 bp

parser = argparse.ArgumentParser(description="parser for coverage file using GFF or PTT files")
parser.add_argument('GENES_FILE', help="path to the gene file (GFF or PTT) you want to use for genes")
args = vars(parser.parse_args())

gene_file_name = args['GENES_FILE']
gene_filetype = gene_file_name.split('.')[1]
if gene_filetype != 'gff' and gene_filetype != 'ptt':
    print "Acceptable file types for GENES_FILE are: .gff, .ptt"
    sys.exit()
print "file type: " + gene_filetype
GENES_FILE = open(gene_file_name, 'rb')
MIN_LENGTH = 90
numShortGenes = 0
numGenes = 0
if gene_filetype == 'ptt':
    #read off the three header lines in a standard .ptt file
    GENES_FILE.readline()
    GENES_FILE.readline()
    GENES_FILE.readline()
    for line in GENES_FILE.read().splitlines():
        numGenes += 1
        values = line.split('\t')
        start = int(values[0].split('..')[0])
        end = int(values[0].split('..')[1])
        if (end - start + 1 < MIN_LENGTH):
            numShortGenes += 1
elif gene_filetype == "gff":
    for line in GENES_FILE.read().splitlines():
        numGenes += 1
        values = line.split('\t')
        start = int(values[3])
        end = int(values[4])
        if (end - start + 1 < MIN_LENGTH):
            numShortGenes += 1
print "number of short genes: " + str(numShortGenes)
print "number of long genes: " + str(numGenes - numShortGenes)

