import argparse
import sys
from Bio import SeqIO

#TODO find out why it's only getting ORFs on negative strand

parser = argparse.ArgumentParser(description="parser for coverage file using GFF or PTT files")
parser.add_argument('SEQUENCE_FILE', help="path to the sequence file you want to parse")
parser.add_argument('GENES_FILE', help="path to the gene file (GFF or PTT) you want to use for genes")
args = vars(parser.parse_args())

sequence_file_name = args['SEQUENCE_FILE']
print sequence_file_name
SEQUENCE_FILE = open(sequence_file_name, 'r')
BIG_GENOME_STRING = SEQUENCE_FILE.read()

gene_file_name = args['GENES_FILE']
gene_filetype = gene_file_name.split('.')[1]
if gene_filetype != 'gff' and gene_filetype != 'ptt':
    print "Acceptable file types for GENES_FILE are: .gff, .ptt"
    sys.exit()
print "file type: " + gene_filetype
GENES_FILE = open(gene_file_name, 'rb')
gene_coords = []
if gene_filetype == 'ptt':
    #read off the three header lines in a standard .ptt file
    GENES_FILE.readline()
    GENES_FILE.readline()
    GENES_FILE.readline()
    for line in GENES_FILE.read().splitlines():
        values = line.split('\t')
        start = int(values[0].split('..')[0])
        end = int(values[0].split('..')[0])
        gene_coords.append((start, end))
elif gene_filetype == "gff":
    for line in GENES_FILE.read().splitlines():
        values = line.split('\t')
        start = int(values[3])
        end = int(values[4])
        gene_coords.append((start, end))

NCBI_TABLE = 11
MIN_ORF_LENGTH = 90

answer = []
i = 0
for seq_record in SeqIO.parse(sequence_file_name, "fasta"):
    seq_len = len(seq_record.seq)
    for strand, nuc in [(+1, seq_record.seq), (-1, seq_record.seq.reverse_complement())]:
        i += 1
        if i % 1000 == 0:
            print i
        for frame in range(3):
            trans = str(nuc[frame:].translate(NCBI_TABLE))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            #length = 3 * ((len(seq_record) - frame) // 3) #multiple of three
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end - aa_start >= MIN_ORF_LENGTH / 3:
                    if strand == 1:
                        start = frame + aa_start * 3
                        end = min(seq_len, frame + aa_end * 3 + 3)
                    else:
                        start = seq_len - frame - aa_end * 3 - 3
                        end = seq_len - frame - aa_start * 3
                    answer.append((start, end, strand, trans[aa_start : aa_end]))
                aa_start = aa_end + 1
            #for pro in nuc[frame:frame + length].translate(NCBI_TABLE).split("*"):
                #if len(pro) >= MIN_ORF_LENGTH / 3:
                    #print("%s...%s - length %i, strand %i, frame %i" % (pro[:30], pro[-3:],  len(pro), strand, frame))
    #print seq_record.id
    #print repr(seq_record.seq)
    #print len(seq_record)
answer.sort()
for start, end, strand, prot in answer:
    if strand == "1":
        print("%s...%s - length %i, strand %i, %i:%i" % (prot[:30], prot[-3:],  len(prot), strand, start, end))


















