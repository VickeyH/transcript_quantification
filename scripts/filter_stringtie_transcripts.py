import sys
import os
import re

# infile
infile = sys.argv[1]

# outfile
outfile = open(sys.argv[2], 'w')

# outfile header
outfile.write('\t'.join(['#CHROM', 'CALLER', 'FEATURE', 'START', 'END', 'COV', 'FPKM', 'TPM', 'STATUS']) + '\n')

# Iterate throuhg infile and grab fields
for line in open(infile):
	line = line.rstrip()
	col = line.split('\t')
	# Assign vars
	CHROM = col[0]
	CALLER = col[1]
	FEATURE = col[2]
	START = col[3]
	END = col[4]
	ATTR = col[8]
	# Grab vars for transcripts
	if FEATURE == 'transcript':
		#print(ATTR)
		COV = ATTR.split('cov "')[1].split('"')[0]
		FPKM = ATTR.split('FPKM "')[1].split('"')[0]
		TPM = ATTR.split('TPM "')[1].split('"')[0]
		if re.search('reference', ATTR):
			STATUS = 'Known'
		else:
			STATUS = 'Novel'
		outfile.write('\t'.join([CHROM, CALLER, FEATURE, START, END, COV, FPKM, TPM, STATUS]) + '\n')
outfile.close()
