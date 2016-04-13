import sys
from trimOneNT import trim
from initialTrimmedAlignments import align_reads

MIN_LENGTH_FOR_ALIGNMENT = 20
base_dir = "/storage/james/data/ecoli/reads/2015study/"

def main():
    num_remaining = sys.maxint
    numTrimmed = 0
    while (num_remaining > 0):
        print "%%%%%%%%%%%%%%%%% " + str(numTrimmed) + " trimmed %%%%%%%%%%%%%%%%%%%%%"
        #first, trim
        trimReads(numTrimmed)
        #then, align the unaligned
        numTrimmed += 1
        align_reads(base_dir + "trim" + str(numTrimmed))
        for i in range(1734430, 1734445):
            FASTA_NAME = base_dir + "trim" + str(numTrimmed) + "/SRR" + str(i) + "_unaligned.fasta"
            p = Popen(["wc", "-l", FASTA_NAME], stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            NUM_READS = int(out.split(' ')[0])/2
            num_remaining = NUM_READS
            print "number of reads remaining for SRR" + str(i) + ": " + str(num_remaining)

def trimReads(numTrimmed):
    for i in range(1734430, 1734445):
        print "$" * 10 + " trimming SRR" + str(i) + " " + "$" * 10
        in_file = base_dir + "trim" + str(numTrimmed) + "/SRR" + str(i) + "_unaligned.fasta"
        out_file = base_dir + "trim" + str(numTrimmed + 1) + "/SRR" + str(i) + ".fasta"
        trim(in_file, out_file, MIN_LENGTH_FOR_ALIGNMENT)
 

if __name__ == "__main__":
    main()    
