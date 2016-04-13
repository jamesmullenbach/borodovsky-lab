import argparse
import sys

#parser = argparse.ArgumentParser(description="parser for which file you want to trim")
#parser.add_argument("READS_FILE", help="path to the fasta format reads file you wish to trim")

#args = vars(parser.parse_args())
LINKER_SEQ = "CTGTAGGCACCATCAAT"

#read must match at least this length of the linker sequence in order to be trimmed and kept
MIN_MATCH_LENGTH = 10
#create list of subsequences so we don't have to create them every time
LINKER_SUBSEQS = [LINKER_SEQ[:i] for i in range(MIN_MATCH_LENGTH, 17)]

#trimmed sequence must be at least this length to be kept
MIN_READ_LENGTH = 20

for i in range(1734430, 1734445):
    print "########## SRR " + str(i) + " ##########"
    IN_FILE_NAME = "SRR" + str(i) + ".fasta"

    OUT_FILE_NAME = IN_FILE_NAME.replace(".fasta", "_trimmed.fasta")

    num_reads = 0
    num_kept = 0
    with open(IN_FILE_NAME, 'r') as IN_FILE:
        with open(OUT_FILE_NAME, 'w') as OUT_FILE:
            j = 0
            for line in IN_FILE:
                if j < 0:
                    print "original line: " + line.rstrip()
                #trim linker sequence
                matched = False
                loc = line.find(LINKER_SEQ)
                if loc != -1:
                    #trim the sequence
                    matched = True
                    line = line[:loc]
                else:
                    #search for partial match at end
                    for i in range(len(LINKER_SUBSEQS)):
                        l = len(LINKER_SUBSEQS[i])
                        ind = -1*l-1
                        if j < 0:
                            print "read: " + line[ind:].rstrip() + ", subseq: " + LINKER_SUBSEQS[i]
                        if line[ind:].rstrip() == LINKER_SUBSEQS[i]:
                            line = line[:ind].rstrip()
                            matched = True
                            break
                if matched:
                    if j < 0:
                        print "trimmed line: " + line.rstrip()
                    if len(line.rstrip()) > MIN_READ_LENGTH:
                        OUT_FILE.write(line + '\n') 
                        num_kept += 1
                else:
                    if j < 0:
                        print "no match found!"
                j += 1
                num_reads += 1
                if j % 1000000 == 0:
                    print "% kept: " + str(num_kept / float(num_reads))
