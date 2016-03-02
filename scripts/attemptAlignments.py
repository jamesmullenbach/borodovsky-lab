from subprocess import call

BOWTIE_HOME = "/home/tool/bowtie2"

TRIM_SCRIPT = "/storage/james/data/ecoli/reads/SRR1734441/temp/trimFastqReads.py"

BOWTIE_COMMAND = BOWTIE_HOME + "/bowtie2"
SEQUENCE_PREFIX = "NC_000913"
FASTQ_NAME = "/storage/james/data/ecoli/reads/SRR1734441/temp/SRR1734441_sample.fastq"
OUT_NAME = FASTQ_NAME.split('.')[0] + ".sam"

NUM_READS = 5000

#start with 16 nt intervals, move up
for l in range(16, 50):
    #scan with this interval over entire 50nt read range
    start = 0
    stop = l
    #trim
    print "interval: " + str(start) + ", " + str(stop)
    call(["python", TRIM_SCRIPT, FASTQ_NAME, str(start), str(stop), str(NUM_READS)])

    #now we know the name of the output file, use that for alignment
    FILE_NAME = FASTQ_NAME.split('.')[0] + "_" + str(start) + "_" + str(stop) + ".fastq"
    call([BOWTIE_COMMAND, "-x", SEQUENCE_PREFIX, "-U", FILE_NAME, "-S", OUT_NAME]) 
