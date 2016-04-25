import sys
import math
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna

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

#handy script to print the subsequence of a given FASTA file over the given interval
def getSubsequence(sequence_file, start, end, strand): 
    SEQUENCE_FILE = open(sequence_file, 'r')
    header = SEQUENCE_FILE.readline()
    sequence_list = SEQUENCE_FILE.readlines()
    sequence = ''.join(sequence_list).replace("\n", "")
    if strand == "+":
        res = sequence[start-1:end]
    elif strand == "-":
        #make it a Seq object and reverse complement
        seq = Seq(sequence[start-1:end], unambiguous_dna)
        res = str(Seq.reverse_complement(seq))
    return res

#get coverage from the specified file over the specified interval
#requires that user give the appropriate coverage file for that strand
#return an list of integers
def getIntervalCoverage(coverage_file, start, end): 
    COVERAGE_FILE = open(coverage_file, 'r')
    coverage_list = COVERAGE_FILE.readlines()
    vals = []
    line = 0
    last = start
    list_ind = 0
    while (list_ind < end and line < len(coverage_list)):
        line_vals = coverage_list[line].split('\t')
        list_ind = int(line_vals[1])
        val = int(line_vals[3])
        if list_ind >= start and list_ind < end:
            vals.extend([0] * (list_ind - last - 1))
            last = list_ind
            vals.append(val)
        line += 1
    vals.extend([0] * (end - last - 1))
    return vals

#takes in start location, the range to look through in either direction, and the strand of the gene we want to investigate
#returns a list of integer locations of in-frame start codons (relative to the beginning of the zone investigated
#e.g. if zone size is 51, "start" position will be at index 50), or an empty list if none were found
def checkForStartCodons(start, zone_size, strand):  
    zone_size = int(math.ceil(float(zone_size) / 3) * 3)
    sequence_file = "/storage/james/data/ecoli/NC_000913.fa"
    seq = getSubsequence(sequence_file, start - zone_size, start + zone_size, strand)
    start_codons = ["ATG", "GTG", "TTG"]
    #go through the seq 3 at a time looking for start codons. guaranteed to be in-frame with start, since we enforce zone_size is a multiple of 3
    i = 0
    locations = []
    while i < len(seq):
        codon = seq[i:i+3]
        if codon in start_codons:
            locations.append(i)
        i += 3
    return locations
