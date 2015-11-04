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

START_CODONS = ["ATG", "GTG", "TTG"]
STOP_CODONS = ["TAA", "TGA", "TAG"]
answer = []
i = 0
for seq_record in SeqIO.parse(sequence_file_name, "fasta"):
    seq_len = len(seq_record.seq)
    for strand, nuc in [(+1, seq_record.seq), (-1, seq_record.seq.reverse_complement())]:
        #iterate over both strands
        for frame in range(3):
            #iterate over 3 frames for each separate strand
            frame_seq = nuc[frame:]
            frame_len = len(frame_seq)
            #initialize starting amino acid index and ending amino acid index
            orf_start = 0
            orf_end = 0
            #this loop finds all orfs of min length or greater and adds them to the answer list
            while orf_start < frame_len:
                #search for orf_start
                i = 0
                isStart = False
                while (not isStart and i <= frame_len - 3):
                    triple = frame_seq[i:i+2]
                    for start in START_CODONS:
                        if triple == start:
                            isStart = True
                            orf_start = i
                    i += 3
                #search for orf_end
                isStop = False
                while (not isStop and i <= frame_len - 3):
                    triple = frame_seq[i:i+2]
                    for stop in STOP_CODONS:
                        if triple == stop:
                            isStop = True
                            orf_end = i
                    i += 3
                if orf_end == -1:
                    orf_end = trans_len
                if orf_end - orf_start >= MIN_ORF_LENGTH:
                    if strand == 1:
                        start = frame + orf_start
                        end = min(seq_len, frame + aa_end)
                    else:
                        start = max(0, seq_len - frame - aa_end)
                        end = seq_len - frame - aa_start
                    #append found ORF to answer list
                    data = (start, end, strand, trans[aa_start : aa_end], frame)
                    answer.append(data)
                orf_start = orf_end + 3
            #for pro in nuc[frame:frame + length].translate(NCBI_TABLE).split("*"):
                #if len(pro) >= MIN_ORF_LENGTH / 3:
                    #print("%s...%s - length %i, strand %i, frame %i" % (pro[:30], pro[-3:],  len(pro), strand, frame))
    #print seq_record.id
    #print repr(seq_record.seq)
    #print len(seq_record)
answer.sort()
for start, end, strand, prot, frame in answer:
    if frame == 1 and strand == 1 and start < 10000:
        print("%s...%s - length %i, strand %i, %i:%i" % (prot[:30], prot[-3:],  len(prot), strand, start, end))

#use SeqIO.write to write to a gene file, or just do it yourself
















