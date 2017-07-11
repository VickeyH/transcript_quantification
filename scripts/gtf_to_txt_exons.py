import sys
import os

# infile
infile = sys.argv[1]

# outfile
outfile = open(sys.argv[2], 'w')

# Variable
transcript = sys.argv[3]

# Iterate throuhg infile and grab fields
for line in open(infile):
	line = line.rstrip()
	col = line.split('\t')
	CHROM = col[0]
	FEATURE = 'exon_'+col[8].split('exon_number "')[1].split('"')[0]+'_'+transcript
	START_1 = col[3]
	END_1 = col[4]
	START_0 = str(int(START_1)-1)
	END_0 = str(int(END_1)-1)
	outfile.write('\t'.join([CHROM, START_1, END_1, FEATURE]) + '\n')
outfile.close()
