# Extract the 20-nt target recognition sequences of the sgRNAs
This is a Python script to extract the 20-nt target recognition sequences of the sgRNAs by decomposing the fastq files of amplicon-seq results with the combination of 20 6-nt barcodes.The amplicon lengths are supposed to be 133-135 bp (suppose that there are small insertions and/or deletions in sgRNAs of some transgenic lines).In the "combination.txt" file, the 6-nt barcodes in forward and reverse primers are structured respectively in the first and second columns; additionally, the third column indicates the wells of each 96-well plate.
