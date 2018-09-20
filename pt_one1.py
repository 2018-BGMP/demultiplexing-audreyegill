#!/usr/bin/env python

def phred(letter):
    '''This takes an ascii character and returns the corresponding phred score'''
    return ord(letter) - 33

# This is an array of file titles
titles = ('read1', 'index1', 'read2', 'index2')
# This is an array of file paths
files= ('/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz', '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz', '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz', '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz')


# Import packages
import numpy as np
import gzip

def bp_means(file):
    # create np array of 0's with name data_title
    data= np.zeros(101)
    n = 0
    # open file using file_path
    with gzip.open(file, 'rt') as f:
        while True:
            line1 = f.readline()
            line2 = f.readline()
            line3 = f.readline()
            line4 = f.readline().strip()
            if not line1:
                break
            for r in range(101):
                score = line4[r]
                q = phred(score)
                data[r]+=q
            n += 1
        return data/n       

dir='/home/agill/'
file=files[0]
print(file)
array = bp_means(file)
print('array created')
np.savetxt("R1.csv", array, delimiter=",")


