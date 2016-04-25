import geneTools
import matplotlib.pyplot as plt

def main():
    d_dict = getCoverageAtDistance()
    #do things, compute probabilities, etc.

def getCoverageAtDistance():
    #first, get the ribosome occupancy data for our 700 verified genes

    VERIFIED_FILE_NAME = "/storage/james/data/ecoli/verified.gff"
    POSITIVE_FILE_NAME = "/storage/james/data/ecoli/reads/2015study/positive.ribo"
    NEGATIVE_FILE_NAME = POSITIVE_FILE_NAME.replace("positive", "negative")

    OFFSET = 50
    NUM_GENES = 769

    verified_coverage = {}

    with open(VERIFIED_FILE_NAME, 'r') as f:
        #parse coverage
        i = 0
        for line in f:
            if i > NUM_GENES:
                break;
            vals = line.split('\t')
            start, stop, strand = int(vals[3]), int(vals[4]), vals[6]
            if i % 10 == 0:
                print str(i) + " genes processed so far"
                print "processing gene: (" + str(start) + ", " + str(stop) + ")"
            if strand == "+":
                verified_coverage[(start, stop)] = geneTools.getIntervalCoverage(POSITIVE_FILE_NAME, start-OFFSET, start+OFFSET)
            elif strand == "-":
                verified_coverage[(start, stop)] = geneTools.getIntervalCoverage(NEGATIVE_FILE_NAME, start-OFFSET, start+OFFSET)
            i += 1
        #dictionary: key is distance from start codon. value is a length-769 array of every y (coverage) value at that distance
        #can then turn the array into a histogram
        d_dict = {}
        for i in range(-50, 51):
            d_dict[i] = []
        for (start, stop) in verified_coverage.keys():
            coverage_vals = verified_coverage[(start, stop)]
            print "coverage vals length: " + str(len(coverage_vals))
            for i in range(len(coverage_vals)):
                d_dict[i - 50].append(coverage_vals[i])
        plt.figure()
        plt.hist(d_dict[15], bins=100)
        plt.show()
        
        return d_dict
            
        #print "verified coverage: " + str(verified_coverage)

if __name__ == "__main__":
    main()
