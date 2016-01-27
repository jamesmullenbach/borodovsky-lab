import geneTools
import argparse
import matplotlib.pyplot as plt

#look for the amount of noise by checking coverage values at OFFSETnt upstream from gene stops
parser = argparse.ArgumentParser(description="parser for a gene and a coverage file")
parser.add_argument("GENES_FILE", help="path to gene file you want to use")
parser.add_argument("POSITIVE_COVERAGE_FILE", help="path to coverage file you want to use")

args = vars(parser.parse_args())

gene_lines, gene_filetype = geneTools.readORFLines(args['GENES_FILE'])

pos_lines, neg_lines = geneTools.readCoverageLines(args['POSITIVE_COVERAGE_FILE'])

ind_pos = 0
ind_neg = 0

x_pos = int(pos_lines[ind_pos].split('\t')[1])
x_neg = int(neg_lines[ind_neg].split('\t')[1])

#ribosome profiling reads are performed by recording OFFSET nt upstream of the 3' end
OFFSET = 14

noiseData = []

for line in gene_lines:
    start, end, strand = geneTools.getLineData(line, gene_filetype)
    if strand == "+":
        while x_pos < end - OFFSET and ind_pos < len(pos_lines) - 1:
            ind_pos += 1
            x_pos = int(pos_lines[ind_pos].split('\t')[1])
        if x_pos > end - OFFSET or x_pos < end - OFFSET:
            noiseData.append(0)
        else:
            score = int(pos_lines[ind_pos].split('\t')[3])
            noiseData.append(score)
    elif strand == "-":
        #looking for OFFSETnt upstream from the 3' end, which is actually OFFSET after the start as defined by GFF, PTT
        while x_neg < start + OFFSET and ind_neg < len(neg_lines) - 1:
            ind_neg += 1
            x_neg = int(neg_lines[ind_neg].split('\t')[1])
        if x_neg > start + OFFSET or x_neg < start + OFFSET:
            #first condition means no data for this position, i.e. data is 0
            #second condition means the datapt is beyond the scope of the coverage file, i.e. data is 0
            noiseData.append(0)
        else:
            score = int(neg_lines[ind_pos].split('\t')[3])
            noiseData.append(score)

#PLOTTING
plt.figure()
#plt.plot(noiseData, 'ro')
plt.hist(noiseData, 100)
plt.xlabel('coverage value bin')
plt.ylabel('number of instances')
#plt.xlabel('gene index')
#plt.ylabel('coverage value')
plt.title("Coverage at positions " + str(OFFSET) + "nt upstream from 3' ends")
plt.show() 
