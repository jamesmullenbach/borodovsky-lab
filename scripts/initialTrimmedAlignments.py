from subprocess import *

def align_reads(base_dir):

    BOWTIE_HOME = "/home/tool/bowtie2"

    BOWTIE_COMMAND = BOWTIE_HOME + "/bowtie2"
    SEQUENCE_PREFIX = "NC_000913"

    for i in range(1734430, 1734445):
        FASTA_NAME = base_dir + "SRR" + str(i) + ".fasta"
    
        print "######### SRR" + str(i) + " #############"
        print "counting reads..." 
        p = Popen(["wc", "-l", FASTA_NAME], stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        NUM_READS = int(out.split(' ')[0])/2
        print "number of reads: " + str(NUM_READS) 

        #FILE_NAME = FASTQ_NAME.split('.')[0] + "_" + str(start) + "_" + str(stop) + ".fastq"
        FILE_NAME = FASTA_NAME
        OUT_NAME = FILE_NAME.split('.')[0] + ".sam"
        UNALIGNED_NAME = FILE_NAME.split('.')[0] + "_unaligned.fasta"
        print "------aligning------"
        call([BOWTIE_COMMAND, "-x", SEQUENCE_PREFIX, "-f", "-U", FILE_NAME, "-S", OUT_NAME, "--un", UNALIGNED_NAME, "--no-unal", "-p", "4"])

def main():
    base_dir = "/storage/james/data/ecoli/reads/2015study/"
    align_reads(base_dir)

if __name__ == "__main__":
    main()
