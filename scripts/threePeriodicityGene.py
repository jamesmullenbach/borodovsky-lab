#inspect trends of three-periodicity within a single gene
import sys
import argparse
import numpy
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="parser for ribosome profiling data file and annotated genes file")
parser.add_argument("GENES_FILE", help="path to the annotated genes file you want to parse")
parser.add_argument("POSITIVE_COVERAGE_FILE", help="path to the direct strand ribosome profiling data file you want to parse")

args = vars(parser.parse_args())

#open files
POSITIVE_COVERAGE_FILE = open(args['POSITIVE_COVERAGE_FILE'], 'r')
pos_lines = POSITIVE_COVERAGE_FILE.read().splitlines()

gene_filename = args['GENES_FILE']
gene_filetype = gene_filename.split('.')[1]
if gene_filetype != "gff" and gene_filetype != "ptt":
    print "only GFF or PTT files allowed for GENES_FILE"
    sys.exit(0)
GENES_FILE = open(args['GENES_FILE'], 'r')

#read off leading three lines in ptt format
if gene_filetype == "ptt":
    GENES_FILE.readline()
    GENES_FILE.readline()
    GENES_FILE.readline()

GENE_IND = 100
THRESHOLD = 5000

line = GENES_FILE.readline()

for i in range(GENE_IND):
    line = GENES_FILE.readline()
vals = line.split('\t')
strand = vals[1] if gene_filetype == "ptt" else vals[6]

#ensure positive gene so we don't have to open the negative file
while (strand != "+"):
    line = GENES_FILE.readline()
    vals = line.split('\t')
    strand = vals[1] if gene_filetype == "ptt" else vals[6]
start = int(vals[0].split('..')[0]) if gene_filetype == "ptt" else int(vals[3])
end = int(vals[0].split('..')[1]) if gene_filetype == "ptt" else int(vals[4])

ind_pos = 0

x_pos = int(pos_lines[ind_pos].split('\t')[1])

ones = []
twos = []
threes = []
one_inds = []
two_inds = []
three_inds = []
while x_pos < start:
    ind_pos += 1
    x_pos = int(pos_lines[ind_pos].split('\t')[1])
while x_pos < end:
    vals = pos_lines[ind_pos].split('\t')
    score = int(vals[3])
    x_pos = int(vals[1])
    if score < THRESHOLD:
        if x_pos % 3 == 0:
            three_inds.append(x_pos)
            threes.append(score)
        elif x_pos % 3 == 1:
            one_inds.append(x_pos)
            ones.append(score)
        elif x_pos % 3 == 2:
            two_inds.append(x_pos)
            twos.append(score)
    ind_pos += 1

#PLOTTING
plt.figure()
one_handle, = plt.plot(one_inds, ones, 'r', label="ones")
two_handle, = plt.plot(two_inds, twos, 'b', label="twos")
three_handle, = plt.plot(three_inds, threes, 'g', label="threes")
plt.xlabel("Base")
plt.ylabel("Coverage")
plt.title("Coverage over a single gene for 1st, 2nd, 3rd bases")
#plt.legend([one_handle, two_handle, three_handle])
plt.show()




 
