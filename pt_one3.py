#!/usr/bin/env python

def phred(letter):
    '''This takes an ascii character and returns the corresponding phred score'''
    return ord(letter) - 33

# This is an array of file titles
titles = ('read1', 'index1', 'read2', 'index2')
# This is an array of file paths
#files = ('/mnt/c/Users/effyg/gdrive2/Bi622/De-multiplex/test_R1.fastq','/mnt/c/Users/effyg/gdrive2/Bi622/De-multiplex/test_R2.fastq', '/mnt/c/Users/effyg/gdrive2/Bi622/De-multiplex/test_R3.fastq', '/mnt/c/Users/effyg/gdrive2/Bi622/De-multiplex/test_R4.fastq')
files= ('/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz', '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz', '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz', '/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz')


# Import packages
import numpy as np

def bp_means(file, n):
    # create np array of 0's with name data_title
    data= np.zeros(101)
    # open file using file_path
    with open(file, 'r') as f:
        while True:
            line1 = f.readline()
            line2 = f.readline()
            line3 = f.readline()
            line4 = f.readline().strip()
            if not line1:
                break
            # rows correspond to nucleotides in data_title
            # and are denoted by r
            for r in range(101):
                q = phred(line4[r])
                data[r]+=q
        data = data/n
        return data       

n = 4000000
#dir='/mnt/c/Users/effyg/gdrive2/Bi622/De-multiplex/'
dir='/home/agill/'
file=files[2]
print(file)
array = bp_means(file, n)
print('array created')
numpy.savetxt("R3.csv", array, delimiter=",")


