from subprocess import call

root = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR173/"
for i in range(1734445, 1734445):
    fullpath = root + "SRR" + str(i) + "/SRR" + str(i) + ".sra"
    call(["wget", fullpath])
