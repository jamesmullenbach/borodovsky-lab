import sys

def trim(FASTA_FILE_NAME, OUT_FILE_NAME, MIN_LENGTH):

    with open(FASTA_FILE_NAME, 'r') as IN_FILE:
        with open(OUT_FILE_NAME, 'w') as OUT_FILE:
            skip = True
            for line in IN_FILE:
                if skip:
                    skip = False
                    continue
                if not skip:
                    subseq = line.rstrip()[:-1]
                    if len(subseq) > MIN_LENGTH:
                        OUT_FILE.write(">\n")
                        OUT_FILE.write(subseq + "\n")
                    skip = True

if __name__ == "__main__": 
    OUT_FILE_NAME = FASTA_FILE_NAME.split('.')[0] + "_cutOne.fasta"
    trim(sys.argv[1], OUT_FILE_NAME)
