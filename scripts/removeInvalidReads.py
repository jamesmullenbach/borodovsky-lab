import argparse
import sys
import math

parser = argparse.ArgumentParser(description="pass in the .fastq file you wish to filter")
parser.add_argument("FASTQ_FILE", help="FASTQ file you want to filter")
parser.add_argument("OUT_FILE", help="output file name (FASTQ) that you want to write to")
parser.add_argument("NUM_READS", help="number of reads to consider", type=int)

args = vars(parser.parse_args())

OUT_FILE = open(args["OUT_FILE"], 'w')
FASTQ_FILE = open(args["FASTQ_FILE"], 'r')
NUM_READS = args["NUM_READS"]

read_lengths = []
for i in range(0, NUM_READS):
    line = FASTQ_FILE.readline().rstrip()
    subseq = FASTQ_FILE.readline().rstrip()
    line2 = FASTQ_FILE.readline().rstrip()
    quality = FASTQ_FILE.readline().rstrip()
 
    #print subseq
    #include only those sequences that end in 10 or more A's in output
    #and don't have 10 or more sequential A's anywhere else
    if ("A" * 10) in subseq[-10:]:
        if ("A" * 10 not in subseq[:len(subseq)-10]):
            trimmedSubseq = ""
            trimmedQuality = ""
            #now remove all A's sitting at the end
            foundNonA = False
            for j in xrange(len(subseq)-1, -1, -1):
                c = subseq[j]
                q = quality[j]
                if not foundNonA:
                    if c != 'A':
                        foundNonA = True
                        trimmedSubseq += c
                        trimmedQuality += q
                else:
                    trimmedSubseq += c
                    trimmedQuality += q
            trimmedSubseq = trimmedSubseq[::-1]
            trimmedQuality = trimmedQuality[::-1]
            
            vals = line.split(' ')
            length = vals[2]
            vals[2] = "length=" + str(len(trimmedSubseq))

            vals2 = line2.split(' ')
            length = vals2[2]
            vals2[2] = "length=" + str(len(trimmedQuality))
            
            OUT_FILE.write(' '.join(vals) + '\n')
            OUT_FILE.write(trimmedSubseq + '\n')
            OUT_FILE.write(' '.join(vals2) + '\n')
            OUT_FILE.write(trimmedQuality + '\n')

            read_lengths.append(len(trimmedSubseq))
    if i % 1000000 == 0 and i > 0:
        print("processed " + str(i) + " reads")
mean_length = sum(read_lengths) / float(len(read_lengths))
print "mean read length: " + str(mean_length)

var = 0
for rl in read_lengths:
    var += (rl - mean_length) ** 2
var = var / float(len(read_lengths))
std = math.sqrt(var)
print "read length STD: " + str(std)



