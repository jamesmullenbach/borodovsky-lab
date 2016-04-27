import argparse
import geneTools

parser = argparse.ArgumentParser(description="parser for removing one set of (verified) genes from another set of genes")
parser.add_argument("VERIFIED_GENES_FILE", help="path to file containing verified genes")
parser.add_argument("GENES_FILE", help="path to file containing another set of genes")
parser.add_argument("OUT_FILE", help="path to output file you would like to create")

args = vars(parser.parse_args())

v_lines, v_filetype = geneTools.readORFLines(args["VERIFIED_GENES_FILE"])
lines, filetype = geneTools.readORFLines(args["GENES_FILE"])

toRemove = []

for v_line in v_lines:
    v_start, v_end, v_strand = geneTools.getLineData(v_line, v_filetype)
    for line in lines:
        start, end, strand = geneTools.getLineData(line, filetype)
        if v_start == start and v_end == end and v_strand == strand:
            toRemove.append((start, end, strand))

#have list of (start, end, strand) tuples to remove
#now go through genes file and add all lines to ouput file unless it is in toRemove
with open(args["OUT_FILE"], 'w') as of:
    for line in lines:
        start, end, strand = geneTools.getLineData(line, filetype)
        if (start, end, strand) not in toRemove:
            of.write(line + '\n') 
