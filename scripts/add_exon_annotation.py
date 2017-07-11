import os
import sys

# infile
infile = sys.argv[1]

# outfile
outfile = open(sys.argv[2], 'w')

# reference files
trans_201 = './data/IKZF1_transcripts_IKZF1-201_exon_overlap.txt'
trans_202 = './data/IKZF1_transcripts_IKZF1-202_exon_overlap.txt'

# Write outfile header
outfile.write('\t'.join(['#CHROM', 'CALLER', 'FEATURE', 'START', 'END', 'COV', 'FPKM', 'TPM', 'STATUS', 'EXONS-201', 'EXONS-202']) + '\n')

# Iterate through file, grab fields, annotate with exon info
for line in open(infile):
	if not line.startswith('#'):
		line = line.rstrip()
		col = line.split('\t')
		CHROM = col[0]
		CALLER = col[1]
		FEATURE = col[2]
		START = col[3]
		END = col[4]
		COV = col[5]
		FPKM = col[6]
		TPM = col[7]
		STATUS = col[8]
		VAR = CHROM + START + END + STATUS
		# Iterate through 201 exons file
		for l_201 in open(trans_201):
			l_201 = l_201.rstrip()
			col_201 = l_201.split('\t')
			CHROM_201 = col_201[0]
			START_201 = col_201[1]
			END_201 = col_201[2]
			STATUS_201 = col_201[4]
			VAR_201 = CHROM_201 + START_201 + END_201 + STATUS_201
			if VAR == VAR_201:
				FEATURE_201 = col_201[3]
				EXON_NUMBER_201 = col_201[5]
				EXONS_201 = col_201[6]
				# Iterate through 202 exons file
				for l_202 in open(trans_202):
					l_202 = l_202.rstrip()
					col_202 = l_202.split('\t')
					CHROM_202 = col_202[0]
					START_202 = col_202[1]
					END_202 = col_202[2]
					FEATURE_202 = col_202[3]
					STATUS_202 = col_202[4]
					VAR_202 = CHROM_202 + START_202 + END_202 + STATUS_202
					if VAR == VAR_202:
						EXON_NUMBER_202 = col_202[5]
						EXONS_202 = col_202[6]
						outfile.write('\t'.join([CHROM, CALLER, FEATURE_201, START, END, COV, FPKM, TPM, STATUS, EXONS_201, EXONS_202]) + '\n')
outfile.close()
