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
  *  these records have "reason_invlaid" added to the header
* index reads are finally evaluated for index hopping

### Tracked statistics and output

### Input
