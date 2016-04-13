import argparse
from geneTools import getIntervalCoverage
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="foo")
parser.add_argument("COVERAGE_FILE", help="link to fully specified coverage file you want to use")
parser.add_argument("START", help="start position in the genome")
parser.add_argument("END", help="stop position in the genome")

args = vars(parser.parse_args())
start = int(args['START'])
end = int(args['END'])


vals = getIntervalCoverage(args['COVERAGE_FILE'], start, end)

plt.figure()
x = range(start, end)

print "x length: " + str(len(x)) + " vals length: " + str(len(vals))

if len(vals) < len(x):
    vals.extend([0] * (len(x) - len(vals)))

plt.plot(x, vals, 'ro')
plt.plot(start, 0, 'b*')
plt.xlabel("genome position")
plt.ylabel("coverage score")
plt.title("Coverage over the range specified from the given file")
plt.show()
