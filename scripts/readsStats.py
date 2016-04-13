import numpy

for i in range(1734430, 1734445):
    IN_FILE_NAME = "SRR" + str(i) + "_trimmed.fasta"
    with open(IN_FILE_NAME, 'r') as IN_FILE:
        lens = []
        for line in IN_FILE:
            lens.append(len(line.rstrip()))
        print "mean read length: " + str(numpy.mean(lens))
        print "median read length: " + str(numpy.median(lens))
