# Demultiplexing

### Description
This program:

* evaluates index records for Ns 
  * if the barcode for either read contains one or more Ns, the read records are written out to "undetermined" fastq files
  * these records  have "reason_N" added to the header
* index records without Ns are then evaluated for quality of barcode reads
  * if either read has a quality score lower than the minimum specified by the user (using the -q option), the read records are written out to "undetermined" fastq files
  * these records have "reason_bad_qual" added to the header
* index records without Ns and that have minimum quality are then evaluated for valid barcodes
  * if either the barcode from index 1 or the reverse complement of the barcode from index 2 are not in the indexes input file, the read records are written out to "undetermined" fastq file
  *  these read records have "reason_invlaid" added to the header
* index reads are finally evaluated for index hopping
  *  barcodes at this stage that are not reverse complements are indicative of index hopping
  *  these read records have reason_hopped added to their header and are written out to the "undetermined" fastq
*  all other reads are written to a fastq for their respective sample

Barcodes are added to all header lines. 

### Tracked statistics and output
This program tracks a number of statistics (see below for unit test output example). Demultiplexed fastq files follow the naming convention R1_SampleID_barcode.fq and undetermined_R1.fq

```
Records read: 10
Number of reads with hopped indexes: 1
Number of reads with good indexes: 1
Number of indexes with Ns: 3
Number of indexes with low quality: 2
Rate of index hopping: 0.1
Barcode Count   Percent
GTAGCGTA        1       0.1
CTAGCTCA        0       0.0
```

### Input
Input fastq files must either all be zipped or all unzipped. The user need not specify which.
```
usage: demulti_complex.py [-h] -r1 READ1 -r2 READ2 -i1 INDEX1 -i2 INDEX2 -id
                          INDEXES [-q MINQUAL]

optional arguments:
  -h, --help            show this help message and exit
  -r1 READ1, --read1 READ1
                        Absolute path to read1 fastq.
  -r2 READ2, --read2 READ2
                        Absolute path to read2 fastq.
  -i1 INDEX1, --index1 INDEX1
                        Absolute path to index1 fastq.
  -i2 INDEX2, --index2 INDEX2
                        Absolute path to index2 fastq.
  -id INDEXES, --indexID INDEXES
                        Absolute path to index ID file.
  -q MINQUAL, --minqual MINQUAL
                        Minimum acceptable index quality score (integer). The
                        default is 30.
```                        
