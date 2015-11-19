import sys
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna

if (len(sys.argv) < 3):
    print "Usage: python " + sys.argv[0] + "SEQUENCE_FILE GENE_FILE OUTPUT_FILE_NAME"
    sys.exit()

fasta_name = sys.argv[1]
BIG_GENOME_FILE = open(fasta_name, 'r')
# don't read first line into the sequence
description = BIG_GENOME_FILE.readline()
sequence_list = BIG_GENOME_FILE.read().splitlines()
sequence = ""
for line in sequence_list:
    sequence += line

gene_file_name = sys.argv[2]
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

out_name = sys.argv[3]
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
                break
        i += 3
    numWithStart += 1
    if found:
        strand_str = "POSITIVE" if strand == "+" else "NEGATIVE"
        out_line = ">" + gene_id + "|start=" + str(start + this_start + 1) + "|end=" + str(end)
        out_line += "|strand=" + strand_str + "\n" + orf[this_start:len(orf)] + '\n'
        OUT_FILE.write(out_line)

print "number of ORFs: " + str(num)
print "number with start codons: " + str(numWithStart)
