import sys
import os
import re

# transcript file
t_file = sys.argv[1]

# exon file
e_file = sys.argv[2]

# outfile
outfile = open(sys.argv[3], 'w')

# Dictionary of transcripts and exons
e_hits = {}

# Iterate throuhg infile and grab fields
for t_line in open(t_file):
	t_line = t_line.rstrip()
	t_col = t_line.split('\t')
	t_CHROM = t_col[0]
	t_FEATURE = t_col[2]
	t_START = t_col[3]
	t_END = t_col[4]
	t_STATUS = t_col[5]
	for e_line in open(e_file):
		e_line = e_line.rstrip()
		e_col = e_line.split('\t')
		e_CHROM = e_col[0]
		e_START = e_col[1]
		e_END = e_col[2]
		e_FEATURE = e_col[3]
		if re.search('201', e_FEATURE):
			e_ID = e_FEATURE.split('_IKZF1-201')[0].split('exon_')[1]
		elif re.search('202', e_FEATURE):
			e_ID = e_FEATURE.split('_IKZF1-202')[0].split('exon_')[1]
		# Count exons that hit
		if int(e_START) >= int(t_START) and int(e_END) <= int(t_END):
			if t_FEATURE not in e_hits.keys():
				e_hits[t_FEATURE] = set()
				e_hits[t_FEATURE].add(e_ID)
			if t_FEATURE in e_hits.keys():
				e_hits[t_FEATURE].add(e_ID)
		if t_FEATURE in e_hits.keys():
			EXON_HITS = ';'.join(map(str,sorted(e_hits[t_FEATURE])))
			EXON_NUM = str(len(e_hits[t_FEATURE]))
		else:
			EXON_HITS = 'NA'
			EXON_NUM = '0'
	outfile.write('\t'.join([t_CHROM, t_START, t_END, t_FEATURE, t_STATUS, EXON_NUM, EXON_HITS]) + '\n')
outfile.close()
