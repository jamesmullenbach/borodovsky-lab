import numpy
import matplotlib.pyplot as plt
import argparse

#plot genes proportions on one graph
#no arguments needed

FILE_1 = open('/home/james/borodovsky-lab/Ecoli_strand_data/Li_Gene-Wei_Ecoli_2012/geneCoverageProportions.txt', 'r')
FILE_2 = open('/home/james/borodovsky-lab/Ecoli_strand_data/Liu_Xiaoqiu_Ecoli_lambda_2013/geneCoverageProportions.txt', 'r')
FILE_3 = open('/home/james/borodovsky-lab/Ecoli_strand_data/Oh_Eugene_Ecoli_2011/geneCoverageProportions.txt', 'r')

file_1_lines = FILE_1.read().splitlines()
file_2_lines = FILE_2.read().splitlines()
file_3_lines = FILE_3.read().splitlines()

file_1_data = []
file_2_data = []
file_3_data = []

for line in file_1_lines:
    file_1_data.append(float(line.split('\t')[3]))
for line in file_2_lines:
    file_2_data.append(float(line.split('\t')[3]))
for line in file_3_lines:
    file_3_data.append(float(line.split('\t')[3]))

x_axis = range(len(file_1_lines))

plt.figure()
plt.plot(x_axis, file_1_data, 'r')
plt.plot(x_axis, file_2_data, 'g')
plt.plot(x_axis, file_3_data, 'b')
plt.xlabel("Gene index")
plt.ylabel("Coverage proportion")
plt.title("Gene Coverage proportions - 3 experiments")
plt.show()
