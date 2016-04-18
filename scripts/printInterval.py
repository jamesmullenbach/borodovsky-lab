import sys
import argparse
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna

#handy script to print the subsequence of a given FASTA file over the given interval

parser = argparse.ArgumentParser(description="print the given interval of the given sequence")
parser.add_argument('SEQUENCE_FILE', help="path to the sequence file you want to parse")
parser.add_argument('START', help="the start point of the interval", type=int)
parser.add_argument('END', help="the endpoint of the interval", type=int)
parser.add_argument('STRAND', help="the strand of this interval (1 or -1)", type=int)

args = vars(parser.parse_args())

SEQUENCE_FILE = open(args['SEQUENCE_FILE'], 'r')
header = SEQUENCE_FILE.readline()
sequence_list = SEQUENCE_FILE.readlines()
sequence = ''.join(sequence_list).replace("\n", "")
if args['STRAND'] == 1:
    seq_str = sequence[args['START']-1:args['END']]
    for i, c in enumerate(seq_str):
        print '\t'.join([str(args['START'] + i), c])
elif args['STRAND'] == -1:
    #make it a Seq object and reverse complement
    seq = Seq(sequence[args['START']-1:args['END']], unambiguous_dna)
    seq_str = str(Seq.reverse_complement(seq))
    for i, c in enumerate(seq_str):
        print '\t'.join([str(args['END'] - i), c])
