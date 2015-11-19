import argparse
import sys
from Bio import SeqIO

#TODO find out why it's only getting ORFs on negative strand

parser = argparse.ArgumentParser(description="parser for coverage file using GFF or PTT files")
parser.add_argument('SEQUENCE_FILE', help="path to the sequence file you want to parse")
parser.add_argument('GENES_FILE', help="path to the gene file (GFF or PTT) you want to use for genes")
parser.add_argument('--out', help="directory in which you would like to save the found ORFs in GFF format")

args = vars(parser.parse_args())

sequence_file_name = args['SEQUENCE_FILE']
print sequence_file_name
SEQUENCE_FILE = open(sequence_file_name, 'r')
BIG_GENOME_STRING = SEQUENCE_FILE.read()

gene_file_name = args['GENES_FILE']
gene_file_origname, gene_filetype = gene_file_name.split('.')
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
                last_stop = this_stop
                            #for pro in nuc[frame:frame + length].translate(NCBI_TABLE).split("*"):
                #if len(pro) >= MIN_ORF_LENGTH / 3:
                    #print("%s...%s - length %i, strand %i, frame %i" % (pro[:30], pro[-3:],  len(pro), strand, frame))
    #print seq_record.id
    #print repr(seq_record.seq)
    #print len(seq_record)
answer.sort()
#remove ORFs that are verified genes as found in GENES_FILE
#ORFs are currently stored in answer, sorted by start nucleotide
#classify ORFs as: 
#    1) intergenic (existing wholly between genes)
#    2) in a gene shadow (within a gene that lies on the opposite strand)
#    3) overlapping (anything else)
gene_ind = 0
gene_orfs = []
non_gene_orfs = []

#for every orf, look through genes to classify until a gene is found that can classify this ORF
for i in range(len(answer)):
    #add one to start because gene files have one-base offset
    start, end, strand, prot, frame = answer[i][0] + 1, answer[i][1], answer[i][2], answer[i][3], answer[i][4]
    if (gene_ind >= len(gene_coords)):
        break
 
    #a relevant gene is any gene that doesn't start and end before the ORF starts
    #go to first gene that doesn't end before the ORF starts
    relevant_gene = gene_coords[gene_ind]
    while relevant_gene[1] < start:
        gene_ind += 1
        relevant_gene = gene_coords[gene_ind]

    #save gene data to separate vars 
    gene_start = relevant_gene[0]
    gene_end = relevant_gene[1]
    gene_strand = relevant_gene[2]
    
    if gene_start < start:
        if gene_end < end:
            #this ORF is overlapping
            answer[i] = (start, end, strand, prot, frame, "overlapping")
        elif gene_end >= end:
            if strand == gene_strand:
                #this ORF is overlapping
                answer[i] = (start, end, strand, prot, frame, "overlapping")
            else:
                #this ORF is shadow
                answer[i] = (start, end, strand, prot, frame, "shadow")
    elif gene_start >= start and gene_start <= end:
        if gene_end < end:
            #this ORF is overlapping
            answer[i] = (start, end, strand, prot, frame, "overlapping")
        elif gene_end == end:
            if strand == gene_strand:
                #this ORF is a gene 
                answer[i] = (start, end, strand, prot, frame, "gene")
                gene_orfs.append((start, end, strand, prot, frame, "gene"))
            else:
                #this ORF is overlapping
                answer[i] = (start, end, strand, prot, frame, "overlapping")
        else:
            #this ORF is overlapping
            answer[i] = (start, end, strand, prot, frame, "overlapping")
    else:
        #gene starts and ends after ORF ends
        #this ORF is intergenic
        answer[i] = (start, end, strand, prot, frame, "intergenic")

print "number of genes found: " + str(len(gene_orfs))
for gene_orf in gene_orfs:
    answer.remove(gene_orf)

if args['out'] is not None:
    #write file in GFF format
    print "about to write to file"
    intergenic_file_name = args['out'] + gene_file_origname + "_intergenic.gff"
    overlapping_file_name = args['out'] + gene_file_origname + "_overlapping.gff"
    shadow_file_name = args['out'] + gene_file_origname + "_shadow.gff"
    INTERGENIC_FILE = open(intergenic_file_name, 'w')
    OVERLAPPING_FILE = open(overlapping_file_name, 'w')
    SHADOW_FILE = open(shadow_file_name, 'w')
    for start, end, strand, prot, frame, category in answer:
        #GFF format: NC_000913 \t category \t . \t start \t end \t . \t strand \t frame \t start_codon=
        details = "stop_codon=" + str(prot[-3:])
        strand = "+" if strand == 1 else "-"
        data = ["NC_000913", category, ".", str(start), str(end), ".", str(strand), str(frame), details]
        line = "\t".join(data) + "\n"
        #print "line: " + str(line)
        if category == 'intergenic':
            INTERGENIC_FILE.write(line)
        elif category == 'overlapping':
            OVERLAPPING_FILE.write(line)
        elif category == 'shadow':
            SHADOW_FILE.write(line)


