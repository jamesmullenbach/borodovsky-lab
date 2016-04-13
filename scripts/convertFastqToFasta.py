for i in range(1734435, 1734445):
    print "SRR" + str(i) + "..."
    with open('SRR' + str(i) + '.fastq', 'r') as FASTQ:
        with open('SRR' + str(i) + '.fasta', 'w') as FASTA:
            i = 0
            for line in FASTQ:
                if (i - 1) % 4 == 0:
                    FASTA.write(line)
                i += 1
