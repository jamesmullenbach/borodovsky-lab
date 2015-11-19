import sys
import argparse
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna

parser = argparse.ArgumentParser(description="script takes in a GFF of ORFs and returns only those with in-frame start codons")
parser.add_argument('SEQUENCE_FILE', help="path to the sequence file you want to parse")
parser.add_argument('ORFS_FILE', help="path to the ORF file (GFF or PTT) you want to use for ORFs")
parser.add_argument('--out', help="name of output file to create")
parser.add_argument('--gff', help="name of output file (in gff format) to create")

args = vars(parser.parse_args())

if (len(sys.argv) < 3):
    print "Usage: python " + sys.argv[0] + "SEQUENCE_FILE GENE_FILE OUTPUT_FILE_NAME"
    sys.exit()

fasta_name = args['SEQUENCE_FILE']
BIG_GENOME_FILE = open(fasta_name, 'r')
# don't read first line into the sequence
description = BIG_GENOME_FILE.readline()
sequence_list = BIG_GENOME_FILE.read().splitlines()
sequence = ""
for line in sequence_list:
    sequence += line

gene_file_name = args['ORFS_FILE']
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

out_name = args['out'] if args['out'] is not None else args['gff']
#out_name = os.path.splitext(fasta_name)[0] + "_genes.txt"
OUT_FILE = open(out_name, 'w')
START_CODONS = ("ATG", "TTG", "GTG")
num = len(GENE_LINES)
numWithStart = 0
for line in GENE_LINES:
    values = line.split('\t')
    #parse interval values and gene ids from gff lines
    #subtract one because python indexes from 0 and GFF indexes from 1
    if gene_filetype == "gff":
        start = int(values[3]) - 1
        end = int(values[4])
        details = values[len(values) - 1]
        gene_id = details.split(' ')[0]
        strand = values[6]
        frame = values[7]
    elif gene_filetype == "ptt":
	start = int(values[0].split('..')[0]) - 1
        end = int(values[0].split('..')[1])
        details = ""
        gene_id = ""
        strand = values[1]
    if (strand == '+'):
        orf = sequence[start:end]
    elif (strand == '-'):
        #make it a Seq and complement it
        seq = Seq(sequence[start:end], unambiguous_dna)
        orf = str(seq.reverse_complement())
    #search for a start codon
    found = False
    i = 0
    while (not found and i <= len(orf) - 90):
        triple = orf[i:i+3]
        for start_codon in START_CODONS:
            if triple == start_codon:
                found = True
                this_start = i
                numWithStart += 1
                break
        i += 3
    if found:
        if args['out'] is not None:
            strand_str = "POSITIVE" if strand == "+" else "NEGATIVE"
            out_line = ">" + gene_id + "|start=" + str(start + this_start + 1) + "|end=" + str(end)
            out_line += "|strand=" + strand_str + "\n" + orf[this_start:len(orf)] + '\n'
            OUT_FILE.write(out_line)
        if args['gff'] is not None:
            details = "stop_codon=" + str(orf[-3:]) + "|start_codon=" + str(orf[this_start:this_start+3])
            if frame is None:
                frame = "."
            data = ["NC_000913", "candidate ORF", ".", str(start + this_start + 1), str(end), ".", str(strand), str(frame), details]
            line = "\t".join(data) + "\n"
            OUT_FILE.write(line)


print "number of ORFs: " + str(num)
print "number with start codons: " + str(numWithStart)
