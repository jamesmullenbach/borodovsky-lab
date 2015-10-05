import sys
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna

if (len(sys.argv) < 3):
    print "Usage: python " + sys.argv[0] + " SEQUENCE_FILE.fa GFF_FILE.gff OUTPUT_FILE_NAME"
    sys.exit()

fasta_name = sys.argv[1]
BIG_GENOME_FILE = open(fasta_name, 'r')
# don't read first line into the sequence
description = BIG_GENOME_FILE.readline()
sequence_list = BIG_GENOME_FILE.read().splitlines()
sequence = ""
for line in sequence_list:
    sequence += line

gff_name = sys.argv[2]
GENES_FILE = open(gff_name, 'rb')
out_name = sys.argv[3]
#out_name = os.path.splitext(fasta_name)[0] + "_genes.txt"
OUT_FILE = open(out_name, 'w')
for line in GENES_FILE:
    values = line.split('\t')
    #parse interval values and gene ids from gff lines
    #subtract one because python indexes from 0 and GFF indexes from 1
    start = int(values[3]) - 1
    end = int(values[4])
    details = values[len(values) - 1]
    gene_id = details.split(' ')[0]
    if (values[6] == '+'):
        out_line = ">" + gene_id + "|start=" + values[3] + "|end=" + values[4]
        out_line += "|strand=POSITIVE\n" + sequence[start:end] + '\n'
        OUT_FILE.write(out_line)
    elif (values[6] == '-'):
        #make it a Seq and complement it
        seq = Seq(sequence[start:end], unambiguous_dna)
        out_line = ">" + gene_id + "|start=" + values[3] + "|end=" + values[4]
        out_line += "|strand=NEGATIVE\n" + str(seq.complement()[::-1]) + '\n'
        OUT_FILE.write(out_line)

