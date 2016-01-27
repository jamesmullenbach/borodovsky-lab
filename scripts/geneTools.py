import sys

#return lines and FILE
def readORFLines(filename):
    filetype = filename.split('.')[1]
    if filetype != "ptt" and filetype != "gff":
        print "Gene/ORF file must be in .GFF or .PTT format"
        sys.exit(0)
    with open(filename, 'r') as FILE:
        if filetype == "ptt":
            FILE.readline()
            FILE.readline()
            FILE.readline()
        lines = FILE.read().splitlines()
    return lines, filetype

#read both positive and coverage files and return the lines
def readCoverageLines(positive_filename):
    with open(positive_filename, 'r') as POSITIVE_COVERAGE_FILE: 
        pos_lines = POSITIVE_COVERAGE_FILE.read().splitlines()

    negative_filename = positive_filename.replace("positive", "negative")
    with open(negative_filename, 'r') as NEGATIVE_COVERAGE_FILE: 
        neg_lines = NEGATIVE_COVERAGE_FILE.read().splitlines()

    return pos_lines, neg_lines

#get data from a GFF or PTT line
def getLineData(line, filetype):
    vals = line.split('\t')
    if filetype == "gff":
        start = int(vals[3])
        end = int(vals[4])
        strand = vals[6]
    elif filetype == "ptt":
        start = int(vals[0].split('..')[0])
        end = int(vals[0].split('..')[1])
        strand = vals[1]
    else:
        print "Gene/ORF file must be in .GFF or .PTT format"
        sys.exit(0)
    return start, end, strand
