from subprocess import call

for i in range(1734430, 1734445):
    call(["/home/tool/sratoolkit/bin/fastq-dump", "SRR" + str(i)])
