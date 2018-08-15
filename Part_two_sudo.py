# Part Two Sudo
# But also a lot of actual code



#########################
###                   ###
###    Importing      ###
###                   ###
#########################
import numpy as np


###################################
###                             ###
###        Variable Key         ###
###    ----------------------   ###
###      r1 = read 1            ###
###      i1 = index 1           ###
###      r2 = read 2            ###
###      i2 = index 2           ###
###      b = barcode length     ###
###     minqual = minimum       ###
###               acceptable    ###
###               Q score       ###
###   index_list = list of      ###
###                indexes      ###
###                             ###
###################################

b = 8
minqual = 30
index_list = #blah blah

# initialize no_hop = {barcode1_barcode2 : freq} where freq.type = int
no_hop = dict()

# initialize hop = {barcode1_barcode2 : freq} where freq.type = int
hop = dict()

# store the 24 indexes
indexes = []

fr1 = "a string with the file path for reads1 file"
fi1 = "a string with the file path for index1 file"
fr2 = "a string with the file path for reads2 file"
fi2 = "a string with the file path for index2 file"



#########################
###                   ###
###    Functions      ###
###                   ###
#########################

def phred(letter):
    '''This takes an ascii character and returns the corresponding phred score'''
    return ord(letter) - 33

def qual(string):
    '''This turns a quality seqn into a quality average integer (or float if I feel like it)'''
    
    
def barmatch(barcode1, barcode2):
    '''This takes two barcodes of equal length and determines if they match;
    i.e. if they are reverse complement sequences'''
       


#########################
###                   ###
###    Functions      ###
###                   ###
#########################


with open(fr1, 'r') as r1, open(fi1, 'r') as i1, open(fr2, 'r') as r2, open(fi2, 'r') as i2:
    # line count denoted by the following 
    c = 0
    # count of good indexes
    good_ind = 0
    # count of bad indexes
    bad_ind = 0
    while True:
    
        # record of read 1 fastq file
        r1_head = r1.readline().strip()
        r1_read = r1.readline().strip()
        r1_plus = r1.readline().strip()
        r1_qual = r1.readline().strip()

        # record of index 1 fastq 
        i1_head = i1.readline().strip()
        i1_read = i1.readline().strip()
        i1_plus = i1.readline().strip()
        i1_qual = i1.readline().strip()

        # record of read 2 fastq file
        r2_head = r2.readline().strip()
        r2_read = r2.readline().strip()
        r2_plus = r2.readline().strip()
        r2_qual = r2.readline().strip()

        # record of index 2 fastq 
        i2_head = i2.readline().strip()
        i2_read = i2.readline().strip()
        i2_plus = i2.readline().strip()
        i2_qual = i2.readline().strip()   
        
         
        ########
        ########      N FREE START      ########
        ########
        # start with barcodes that are N free
        if 'N' not in i1_read:
            # sets meanq variable using quality function
            meanq = quality(i1_read)
            
            ########
            ########      MIN QUAL START      ########
            ########
            # next check that minimum quality score (minqual) criteria is met for barcodes
            if meanq > minqual:
            
                ########
                ########      MATCH START     ########
                ########
                #next check for matching scores using barmatch fxn
                if len(i1_read)== b and barmatch(i1_read,i2_read)==TRUE and i1_read in index_list:
                    # write out record to good_1_n and good_2_n where n = index number 
                    # something like good.write() for n in index_list so that there are 48 "good" files
                    # add to good_ind counter
                else:
                    # write record to bad_1_n and bad_2_n where n = index number
                    # same model for 48 "bad" files
                    # add to bad_ind counter
                ########
                ########      MATCH END       ########    
                ########
            else:
                # write record to bad_1_n and bad_2_n where n = index number
            ########
            ########      MIN QUAL END      ########
            ########
        else:
            # write record to bad_1_n and bad_2_n where n = index number
        ########
        ########      N FREE START      ########
        ########    
            
            
        # this prevents an infinite loop    
        if not line1:
            break

        c+=1

        
        
        