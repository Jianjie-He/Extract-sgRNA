# Extract the 20-bp target recognition sequences of the sgRNAs
This is a Python script to extract the 20-bp target recognition sequences of the sgRNAs by decomposing the fastq files of amplicon-seq results with the 6-bp barcodes.The amplicon length is supposed to be 133-135 bp (suppose that there are small insertions and/or deletions in sgRNAs of some transgenics).
In the "combination.txt" file, the 6-bp barcodes in forward and reverse primers are structured respectively in the first and second column; additionally, the third column indicates the well ID of each 96-well plates.
