import argparse
import sys
from Bio import SeqIO

#TODO find out why it's only getting ORFs on negative strand

parser = argparse.ArgumentParser(description="parser for coverage file using GFF or PTT files")
parser.add_argument('SEQUENCE_FILE', help="path to the sequence file you want to parse")
parser.add_argument('GENES_FILE', help="path to the gene file (GFF or PTT) you want to use for genes")
parser.add_argument('--OUT_FILE', help="name of file you would like to save the found ORFs to in GFF format")
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
        end = int(values[0].split('..')[1])
        strand = 1 if values[1] == '+' else -1
        gene_coords.append((start, end, strand))
elif gene_filetype == "gff":
    for line in GENES_FILE.read().splitlines():
        values = line.split('\t')
        start = int(values[3])
        end = int(values[4])
        strand = 1 if values[6] == '+' else -1
        gene_coords.append((start, end, strand))
gene_coords.sort()
print "gene_coords[0][0], gene_coords[0][1]: %i, %i" % (gene_coords[0][0], gene_coords[0][1])
for gene_coord in gene_coords:
    if (gene_coord[0] < 10000):
        print "gene coord: %i, %i" % (gene_coord[0], gene_coord[1])

NCBI_TABLE = 11
MIN_ORF_LENGTH = 90

START_CODONS = ["ATG", "GTG", "TTG"]
START_REGEX = "(A|G|T)TG"
STOP_CODONS = ["TAA", "TGA", "TAG"]
STOP_REGEX = "T(AA|GA|AG)"
answer = []
i = 0
for seq_record in SeqIO.parse(sequence_file_name, "fasta"):
    seq_len = len(seq_record.seq)
    for strand, nuc in [(+1, seq_record.seq), (-1, seq_record.seq.reverse_complement())]:
        #iterate over both strands
        for frame in range(3):
            #iterate over 3 frames for each separate strand
            print "this frame: " + str(nuc[frame:frame + 30]) + "..." + str(nuc[-3:])
            frame_seq = nuc[frame:]
            frame_len = len(frame_seq)
            #initialize first stop codon and second stop codon
            orf_start = 0
            orf_end = 0
            last_stop = 0
            this_stop = 0
            #this loop finds all orfs of min length or greater and adds them to the answer list
            i = 0
            while last_stop < frame_len:
                #search for orf_start
                #search for first available stop codon
                isStop = False
                while (not isStop and i <= frame_len - 3):
                    triple = frame_seq[i:i+3]
                    for stop in STOP_CODONS:
                        if triple == stop:
                            isStop = True
                            this_stop = i + 3
                            break
                    i += 3
                if (i > frame_len - 3):
                    #no start codon found
                    #so break
                    print "no stop codon found"
                    break
                if this_stop - last_stop + 1 >= MIN_ORF_LENGTH:
                    if strand == 1:
                        abs_start = frame + last_stop
                        abs_end = min(seq_len, frame + this_stop)
                    else:
                        abs_start = max(0, seq_len - frame - this_stop)
                        abs_end = seq_len - frame - last_stop
                    #append found ORF to answer list
                    data = (abs_start, abs_end, strand, frame_seq[last_stop : this_stop], frame)
                    if (abs_start < 10000):
                        print data
                        #print "frame: %s, start: %i, end: %i, start_codon %s, stop_codon %s" % (frame, start, end, nuc[start-1:start+2], nuc[end-3:end])
                    #print "adding data to answer: " + str(data)
                    answer.append(data)
                last_stop = this_stop + 3
                            #for pro in nuc[frame:frame + length].translate(NCBI_TABLE).split("*"):
                #if len(pro) >= MIN_ORF_LENGTH / 3:
                    #print("%s...%s - length %i, strand %i, frame %i" % (pro[:30], pro[-3:],  len(pro), strand, frame))
    #print seq_record.id
    #print repr(seq_record.seq)
    #print len(seq_record)
answer.sort()
#remove ORFs that are verified genes as found in GENES_FILE
#ORFs are scurrently stored in answer, sorted by start nucleotide
#classify ORFs as: 
#    1) intergenic (existing wholly between genes)
#    2) overlapping (ending stop codon exists within a gene)
#    3) in a gene shadow (within a gene that lies on the opposite strand)
gene_ind = 0
gene_orfs = []
non_gene_orfs = []
for i in range(len(answer)):
    start, end, strand, prot, frame = answer[i][0], answer[i][1], answer[i][2], answer[i][3], answer[i][4]
    if (gene_ind >= len(gene_coords)):
        break
    #add one to start because gene files have one-base offset
    if end == gene_coords[gene_ind][1] and start < gene_coords[gene_ind][0]:
        #this ORF contains a gene
        print "this ORF contains a gene: start %i, end %i, strand %i, frame %i" % (start, end, strand, frame)
        answer[i] = (start, end, strand, prot, frame, "gene")
        gene_orfs.append((start, end, strand, prot, frame, "gene"))
        gene_ind += 1
    elif start > gene_coords[gene_ind][0]:
        if end > gene_coords[gene_ind][1]:
            #this ORF starts inside a gene and ends outside one
            answer[i] = (start, end, strand, prot, frame, "overlapping")
        elif end < gene_coords[gene_ind][1] and strand == gene_coords[gene_ind][2] * -1:
            #this ORF lies in the shadow of a gene on the opposite strand
            answer[i] = (start, end, strand, prot, frame, "shadow")
        gene_ind += 1
    elif start < gene_coords[gene_ind][0] and end > gene_coords[gene_ind][1]:
        #this ORF starts outside a gene and ends inside one 
        answer[i] = (start, end, strand, prot, frame, "overlapping")
        gene_ind += 1
    elif start > gene_coords[gene_ind][1]:
        #this ORF starts after the end of current gene
        #so move to the next gene
        gene_ind += 1
        if end < gene_coords[gene_ind][0]:
            #this ORF starts after the end of a gene and ends before the start of the next
            #so it is intergenic
            answer[i] = (start, end, strand, prot, frame, "intergenic")

print "number of genes found: " + str(len(gene_orfs))
for gene_orf in gene_orfs:
    answer.remove(gene_orf)

if args['OUT_FILE'] is not None:
    #write file in GFF format
    out_file_name = args['OUT_FILE']
    OUT_FILE = open(out_file_name, 'w')
    for start, end, strand, prot, frame, category in non_gene_orfs:
        #GFF format: NC_000913 \t non-gene ORF \t . \t start \t end \t . \t strand \t frame \t start_codon=
        details = "start_codon=" + str(prot[0:3])
        strand = "+" if strand == 1 else "-"
        #add one to start because gene files have one-base offset
        data = ["NC_000913", "non-gene ORF", ".", str(start+1), str(end), ".", str(strand), str(frame), details]
        line = "\t".join(data) + "\n"
        OUT_FILE.write(line)
#FORWARD STRAND WORKS
#OPPOSITE STRAND NEEDS TO TEST















