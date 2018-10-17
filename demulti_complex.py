#!/usr/bin/env python


#########################
###                   ###
###    Importing      ###
###                   ###
#########################

import argparse
import gzip

###################################################################
###                                                             ###
###    Key                                                      ###
###    ----------------------                                   ###
###    r1 = read 1                                              ###
###    i1 = index 1                                             ###
###    r2 = read 2                                              ###
###    i2 = index 2                                             ###
###    ind = index-sample text file                             ###
###    minqual = minimum acceptable Q score                     ###
###    sample_list = {'sample_id': 'barcode'}                   ###
###    barcode_count = {'barcode' : frequency (integer)}        ###
###    counts = {'count statistic' : frequency (integer)}       ###
###    read1_files = {'barcode' : read 1 open file object       ###
###    read2_files = {'barcode' : read 2 open file object       ###
###                                                             ###
###################################################################

#########################
###                   ###
###    Functions      ###
###                   ###
#########################

def phred(letter):
    '''This takes an ascii character and returns the corresponding phred score'''
    return ord(letter) - 33

def qual(string):
    '''This turns a quality seqn into a quality average'''
    quality = 0
    length = len(string)
    for i in string:
        quality += phred(i)
    return quality/length
    
def rev_comp(barcode):
    ''' This takes a barcode and returns its reverse complement'''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    c = ''
    for i in barcode:
        c = complement[i] + c 
    return c
   
def demulti(r1,r2,i1,i2,ind):
    ''' This takes open file objects for two reads, two indexes, and an index text file
    and returns a dictionary of three dictionaries: 
    counts[tracked_values] = frequency
    sample_list[sample ID] = barcode
    barcode_count[barcode] = frequency
    '''
    # This creates a dictionary of the sample IDs on their corresponding barcodes
    sample_list = dict()
    # And a dictionary of barcodes on counts
    barcode_count = dict()

    # Populate sample_list and barcode_count
    for line in ind:
        thing=line.strip().split("\t")
        barcode = thing[4]
        ID=thing[3]
        barcode_count[barcode]=0
        sample_list[ID] = barcode
    del barcode_count['index sequence']        
    del sample_list['index']

    # Open/create undetermined files
    with open("undetermined_R1.fq", "a") as bad_1, open("undetermined_R2.fq", "a") as bad_2:
        record_count = 0
        good_ind = 0
        hopped_ind = 0
        number_N = 0
        number_badqual = 0
        number_invalid = 0
        
        # Open 48 empty append files and create file dictionaries
        read1_files = dict()
        read2_files = dict()
        for ID in sample_list.keys():
            open1 = open("R1_"+ID+"_"+sample_list[ID]+".fq", "a")
            read1_files[sample_list[ID]] = open1
            open2 = open("R2_"+ID+"_"+sample_list[ID]+".fq", "a")
            read2_files[sample_list[ID]] = open2

            
        while True:

            # record of read 1 fastq file
            r1_head = r1.readline().strip()
            r1_read = r1.readline().strip()
            r1.readline().strip() # plus sign
            r1_qual = r1.readline().strip()

            # record of index 1 fastq 
            i1_head = i1.readline().strip()
            i1_read = i1.readline().strip()
            i1.readline().strip() # plus sign
            i1_qual = i1.readline().strip()

            # record of read 2 fastq file
            r2_head = r2.readline().strip()
            r2_read = r2.readline().strip()
            r2.readline().strip() # plus sign
            r2_qual = r2.readline().strip()

            # record of index 2 fastq 
            i2_head = i2.readline().strip()
            i2_read = i2.readline().strip()
            i2.readline().strip() # plus sign
            i2_qual = i2.readline().strip()   
            
            if len(r1_head)==0:
                break
                
            record_count += 1

            
            ########
            ########      N FREE START      ########
            ########
            # start with barcodes that are N free
            if 'N' in i1_read or 'N' in i2_read:
                bad_1.write(r1_head+":"+i1_read+"\treason_N\n"+r1_read+"\n+\n"+r1_qual+"\n")
                bad_2.write(r2_head+":"+i1_read+"\treason_N\n"+r2_read+"\n+\n"+r2_qual+"\n")            
                number_N += 1
                            
            else:
                
                # set meanq variable using quality function
                meanq1 = qual(i1_qual)
                meanq2 = qual(i2_qual)
                
                
                ########
                ########      MIN QUAL START      ########
                ########
                # next check that minimum quality score (minqual) criteria is met for barcodes
                
                if meanq1 > minqual and meanq2 > minqual:
                
                    ########
                    ########      VALID INDEX START     ########
                    ########
                    #next check for valid barcodes
                    
                    if i1_read in barcode_count.keys() and rev_comp(i2_read) in barcode_count.keys():
                        
                        ########
                        ########      MATCH START     ########
                        ########
                        #next check for matching scores using barmatch fxn
                    
                        if i1_read == rev_comp(i2_read):
                            # write out record to good files
                            read1_files[i1_read].write(r1_head+":"+i1_read+"\n"+r1_read+"\n+\n"+r1_qual+"\n")
                            read2_files[i1_read].write(r2_head+":"+i2_read+"\n"+r2_read+"\n+\n"+r2_qual+"\n")
                            good_ind += 1

                            # Add to barcode dictionary
                            barcode_count[i1_read] += 1
                        else:
                            # write record to bad_1 and bad_2
                            bad_1.write(r1_head+":"+i2_read+"\treason_hopped\n"+r1_read+"\n+\n"+r1_qual+"\n")
                            bad_2.write(r2_head+":"+i2_read+"\treason_hopped\n"+r2_read+"\n+\n"+r2_qual+"\n")
                            # Add to hopped count
                            hopped_ind += 1
                            
                        ########
                        ########      MATCH END       ########    
                        ########
                        
                    else:
                        # write record to bad_1 and bad_2
                        bad_1.write(r1_head+":"+i1_read+"\treason_invalid\n"+r1_read+"\n+\n"+r1_qual+"\n")
                        bad_2.write(r2_head+":"+i2_read+"\treason_invalid\n"+r2_read+"\n+\n"+r2_qual+"\n")
                        # Add to invalid count
                        number_invalid += 1
                    
                    ########
                    ########      VALID INDEX END       ########    
                    ########
                
                else:
                    bad_1.write(r1_head+":"+i1_read+"\treason_bad_qual\n"+r1_read+"\n+\n"+r1_qual+"\n")
                    bad_2.write(r2_head+":"+i2_read+"\treason_bad_qual\n"+r2_read+"\n+\n"+r2_qual+"\n")
                    number_badqual += 1
                
                ########
                ########      MIN QUAL END      ########
                ########
            
            ########
            ########      N FREE END      ########
            ########    
                
            
        
        for barcode in barcode_count.keys():
            read1_files[barcode].close()
            read2_files[barcode].close()
    
    counts = {
        'record_count' : record_count,
        'good_ind' : good_ind,
        'hopped_ind' : hopped_ind,
        'number_N' : number_N,
        'number_badqual' : number_badqual,
        'number_invalid' : number_invalid
        }
    all_out = {
        'counts' : counts,
        'sample_list' : sample_list,
        'barcode_count' : barcode_count
        }
    return all_out     
   
def file_type_operations(fr1, fr2, fi1, fi2):
    '''Checks that input files should be either all zipped or all unzipped, 
    and returns True if zipped, False if unzipped'''
    r1_byte = open(fr1, "rb").read(1)
    r2_byte = open(fr2, "rb").read(1)
    i1_byte = open(fi1, "rb").read(1)
    i2_byte = open(fi2, "rb").read(1)
    same_type = r1_byte == r2_byte == i1_byte == i2_byte
    all_zipped = b'\x1f' == r1_byte == r2_byte == i1_byte == i2_byte
    if same_type==False:
        raise ValueError('Files are not in same format. Input fastq files must either be all zipped or all unzipped.') 
    return all_zipped
           
#########################
###                   ###
###     Argparse      ###
###    and whatnot    ###
###                   ###
#########################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r1","--read1", help="Absolute path to read1 fastq.", type=str, required=True)
    parser.add_argument("-r2","--read2", help="Absolute path to read2 fastq.", type=str, required=True)
    parser.add_argument("-i1", "--index1", help="Absolute path to index1 fastq.", type=str, required=True)
    parser.add_argument("-i2", "--index2", help="Absolute path to index2 fastq.", type=str, required=True)
    parser.add_argument("-id", "--indexID", help="Absolute path to index ID file.", dest="indexes", type=str, required=True)
    parser.add_argument("-q", "--minqual", help="Minimum acceptable index quality score (integer). The default is 30.", dest="minqual", type=int, required=False, default=30)
    return parser.parse_args()

args = main()

fr1 = args.read1
fi1 = args.index1
fr2 = args.read2
fi2 = args.index2
indexes = args.indexes
minqual = args.minqual


all_zipped = file_type_operations(fr1, fr2, fi1, fi2)

if all_zipped:    
    with gzip.open(fr1, 'rt') as r1, gzip.open(fi1, 'rt') as i1, gzip.open(fr2, 'rt') as r2, gzip.open(fi2, 'rt') as i2, open(indexes, 'r') as ind:
        output = demulti(r1,r2,i1,i2, ind)
else:
    with open(fr1, 'r') as r1, open(fi1, 'r') as i1, open(fr2, 'r') as r2, open(fi2, 'r') as i2, open(indexes, 'r') as ind:
        output = demulti(r1,r2,i1,i2, ind)
        


counts = output['counts']
barcode_count = output['barcode_count']
rate = counts['hopped_ind']/counts['record_count']
print("Records read: "+str(counts['record_count']))
print("Number of reads with hopped indexes: "+str(counts['hopped_ind']))
print("Number of reads with good indexes: "+str(counts['good_ind']))
print("Number of indexes with Ns: "+str(counts['number_N']))
print("Number of indexes with low quality: "+str(counts['number_badqual']))
print("Rate of index hopping: "+str(rate))

print("Barcode"+"\t"+"Count"+"\t"+"Percent")
for key in barcode_count:
    barcode_percent = barcode_count[key]/counts['record_count']
    print(str(key)+"\t"+str(barcode_count[key])+"\t"+str(barcode_percent))

        

