import glob
import re
from Bio import SeqIO

#Python version -- 3.8.6
#BioPython version -- 1.79
#re version -- 2.2.1

#This script file should be in the same folder as the fq files.
path = glob.glob('*.fq')

#Barcodes in the forward and reverse primers.
f = open('combination.txt')
line = [i.strip() for i in f.readlines()]
f.close()

barcode = {}
for i in line:
	barcode.update({i.split('\t')[0]:i.split('\t')[1]})


#Get all of the sequences from the fq files using BioPython.
def unpack(path=None):

	fq_dict = SeqIO.index(path,'fastq')
	dock = []

	for j in fq_dict.keys():
		seq = str(fq_dict[j].seq).replace('\n','')
		dock.append(seq)

	return dock


#Extract the 20-bp target recognition sequences of the sgRNAs by regular expression.
def extract(seq=None):

	table = str.maketrans('ATCG','TAGC')
	regex = re.compile("ATTG\w{18,22}GTTTT")
	#Small insertions and deletions happened in the 20-bp target recognition sequences.
	position,result = '',''
	seq = seq[:135]#Suppose that the number of "135" indicates the longest amplicon.

	for key in barcode.keys():
		barcode_f = key.split(' ')[0]
		barcode_r = key.split(' ')[1]
		if (barcode_f in seq[4:14]) and (barcode_r in seq[-14:-4]):
			position = barcode[key]
			if len(regex.findall(seq)) != 0:
				result = regex.findall(seq)[0]
				if 'N' not in result:
					result = result[4:-5]
				else:
					result = ''
	else:
		seq = str.translate(seq,table)[::-1]
		for key in barcode.keys():
			barcode_f = key.split(' ')[0]
			barcode_r = key.split(' ')[1]
			if (barcode_f in seq[4:14]) and (barcode_r in seq[-14:-4]):
				position = barcode[key]
				if len(regex.findall(seq)) != 0:
					result = regex.findall(seq)[0]
					if 'N' not in result:
						result = result[4:-5]
					else:
						result = ''

	return position+'\t'+result



#Calculate the number of read ocunts of each transgenic line.
def calculate():

	for i in path:
		data = unpack(i)
		container = []
		for j in data:
			sgRNA = extract(j)
			if len(sgRNA) > 3:
				container.append(sgRNA)

		name = i.split('.')[0]+'.extract.txt'
		f = open(name,'w')

		dock = set(container)
		for k in dock:
			num = container.count(k)
			f.write(str(num)+'\t'+k+'\n')
			f.flush()

		f.close()

if __name__ == '__main__':
	calculate()