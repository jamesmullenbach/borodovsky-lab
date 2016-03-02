import argparse

parser = argparse.ArgumentParser(description="parser for fastq trimming script")
parser.add_argument("FASTQ_FILE", help="path to fastq file you want to trim")
parser.add_argument("START", help="first index of reads to keep", type=int)
parser.add_argument("STOP", help="last index of reads to keep", type=int)
parser.add_argument("NUM_READS", help="number of reads to consider", type=int)

args = vars(parser.parse_args())

FASTQ_FILE = open(args["FASTQ_FILE"], 'r')
OUT_FILE_NAME = args['FASTQ_FILE'].split('.')[0] + "_" + str(args['START']) + "_" + str(args['STOP']) + ".fastq"
OUT_FILE = open(OUT_FILE_NAME, 'w')

for i in xrange(args['NUM_READS']):
    line1 = FASTQ_FILE.readline().rstrip()
    subseq = FASTQ_FILE.readline().rstrip()
    line2 = FASTQ_FILE.readline().rstrip()
    quality = FASTQ_FILE.readline().rstrip()

    subseq = subseq[args['START']:args['STOP']]
    quality = quality[args['START']:args['STOP']]
    
    vals1 = line1.split(' ')
    vals2 = line2.split(' ')

    vals1[2] = "length=" + str(args['STOP'] - args['START'])
    vals2[2] = "length=" + str(args['STOP'] - args['START'])

    OUT_FILE.write(' '.join(vals1) + '\n')
    OUT_FILE.write(subseq + '\n')
    OUT_FILE.write(' '.join(vals2) + '\n')
    OUT_FILE.write(quality + '\n')
