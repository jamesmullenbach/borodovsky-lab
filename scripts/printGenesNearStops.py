#script to print last NUM nucleotides of genes prior to stop codon
import argparse
import sys
import random
import geneTools
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna


parser = argparse.ArgumentParser(description="parser for a gene and a coverage file")
parser.add_argument("GENES_FILE", help="path to gene file you want to use")
parser.add_argument("SEQUENCE_FILE", help="path to sequence file you want to use")

args = vars(parser.parse_args())


####################TEMP#############

NUM = 6

SEQUENCE_FILE = open(args['SEQUENCE_FILE'], 'r')
#read off header line
SEQUENCE_FILE.readline()
sequence_lines = SEQUENCE_FILE.read().splitlines()
seq = "".join(sequence_lines).replace('\n', '')

GENES_FILE = open(args['GENES_FILE'], 'r')
gene_lines = GENES_FILE.read().splitlines()

for line in gene_lines:
    data = line.split('\t')
    start, end, strand = int(data[0]), int(data[1]), data[2]
    if strand == "+":
        print "gene: " + str(start) + "-" + str(end) + " DIRECT: ..." + seq[end-NUM:end]
    elif strand == "-":
        #make it a Seq object and reverse complement
        subseq = Seq(seq[start-1:start+NUM-1], unambiguous_dna)
        subseq = str(Seq.reverse_complement(subseq))
        print "gene: " + str(start) + "-" + str(end) + " COMPLEMENTARY: ..." + subseq
############END########
sys.exit(0)

gene_lines, gene_filetype = geneTools.readORFLines(args['GENES_FILE'])

NUM_GENES = 100
gene_lines_rand = []
attempts = 0
while len(gene_lines_rand) < NUM_GENES:
    index = random.randint(0, len(gene_lines) - 1)
    gene = geneTools.getLineData(gene_lines[index], gene_filetype)
    if gene not in gene_lines_rand:
        gene_lines_rand.append(gene)
        attempts += 1
#print "attempts: " + str(attempts)
gene_lines_rand.sort()

NUM = 6

SEQUENCE_FILE = open(args['SEQUENCE_FILE'], 'r')
#read off header line
SEQUENCE_FILE.readline()
sequence_lines = SEQUENCE_FILE.read().splitlines()
seq = "".join(sequence_lines).replace('\n', '')

for line in gene_lines_rand:
    start, end, strand = line[0], line[1], line[2]
    if strand == "+":
        print "gene: " + str(start) + "-" + str(end) + " DIRECT: ..." + seq[end-NUM:end]
    elif strand == "-":
        #make it a Seq object and reverse complement
        subseq = Seq(seq[start-1:start+NUM-1], unambiguous_dna)
        subseq = str(Seq.reverse_complement(subseq))
        print "gene: " + str(start) + "-" + str(end) + " COMPLEMENTARY: ..." + subseq
