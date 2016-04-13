for i in range(1734430, 1734445):
    print "########### SRR" + str(i) + " ############"
    IN_FILE_NAME = "SRR" + str(i) + "_trimmed.fasta"
    OUT_FILE_NAME = "SRR" + str(i) + ".fasta"
    with open(IN_FILE_NAME, 'r') as IN_FILE:
        with open(OUT_FILE_NAME, 'w') as OUT_FILE:
            OUT_FILE.write(IN_FILE.readline())
            OUT_FILE.write(IN_FILE.readline())
            for line in IN_FILE:
                OUT_FILE.write(">\n")
                OUT_FILE.write(line)
