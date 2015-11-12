#handy script to print the sequence of a given FASTA file over the given interval
import sys
import argparse
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna

parser = argparse.ArgumentParser(description="print the given interval of the given sequence")
parser.add_argument('SEQUENCE_FILE', help="path to the sequence file you want to parse")
parser.add_argument('START', help="the start point of the interval", type=int)
parser.add_argument('END', help="the endpoint of the interval", type=int)
parser.add_argument('STRAND', help="the strand of this interval", type=int)

args = vars(parser.parse_args())

SEQUENCE_FILE = open(args['SEQUENCE_FILE'], 'r')
header = SEQUENCE_FILE.readline()
sequence_list = SEQUENCE_FILE.readlines()
sequence = ''.join(sequence_list).replace("\n", "")
if args['STRAND'] == 1:
    print sequence[args['START']-1:args['END']]
elif args['STRAND'] == -1:
    seq = Seq(sequence[args['START']-1:args['END']], unambiguous_dna)
    print str(Seq.reverse_complement(seq))
