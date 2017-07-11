#===============================================================================
#
#         FILE: /mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/transcript_quantification/transcript_quantification_main_documentation.sh
#
#        USAGE: for documentation purposes, scripts inside
#
#  DESCRIPTION:  This script serves as a step by step documentation script and development script for quantifying transcripts
#                
# REQUIREMENTS:  ---
#        NOTES:  ---
#       AUTHOR:  Alec Steep, steepale@msu.edu
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  1.0
#      CREATED:  2017.07.06
#     REVISION:  
#===============================================================================

# Based on the results from STRINGTIE, I do not see any significant evidence of dominant negative IKZF1 transcript isoforms



# PROJECT DIRECTORY (TUM Cluster)
proj_dir="/mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/transcript_quantification"
export MDV_DIR="/mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang" # Placed in ~/.bashrc file
cd $proj_dir

# Make proper directories
mkdir -p ./{data,scripts,analysis,jobs}
mkdir ./salmon
mkdir ./data/bam

# Sources to cite/reference:
# Salmon main page: https://combine-lab.github.io/salmon/getting_started/
# Salmon GitHub: https://github.com/COMBINE-lab/salmon
# StringTie main page: https://ccb.jhu.edu/software/stringtie/index.shtml#install
# HISAT2 main page: http://ccb.jhu.edu/software/hisat2/index.shtml

### Install dependencies
# Install Salmon v0.8.2
cd ${HOME}/Apps
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.8.2/Salmon-0.8.2_linux_x86_64.tar.gz
tar xzvf Salmon-0.8.2_linux_x86_64.tar.gz
rm Salmon-0.8.2_linux_x86_64.tar.gz
# Add salmon to PATH (in ${HOME}/.bashrc)
export PATH="/mnt/home/steepale/Apps/Salmon-0.8.2_linux_x86_64/bin:${PATH}"
cd ${proj_dir}

# Install StringTie v1.3.3b
cd ${HOME}/Apps
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.Linux_x86_64.tar.gz
tar xvfz stringtie-1.3.3b.Linux_x86_64.tar.gz
rm stringtie-1.3.3b.Linux_x86_64.tar.gz
cd stringtie-1.3.3b.Linux_x86_64
make release #Although it seems to already be made
export PATH="/mnt/home/steepale/Apps/stringtie-1.3.3b.Linux_x86_64:${PATH}"
cd ${proj_dir}

# Install HISAT2 v2-2.1.0
cd ${HOME}/Apps
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip
rm hisat2-2.1.0-Linux_x86_64.zip
cd hisat2-2.1.0
make
export PATH="/mnt/home/steepale/Apps/hisat2-2.1.0:${PATH}"
cd ${proj_dir}

# Install gffcompare v0.10.1
cd ${HOME}/Apps
wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.10.1.Linux_x86_64.tar.gz
tar xvfz gffcompare-0.10.1.Linux_x86_64.tar.gz
rm gffcompare-0.10.1.Linux_x86_64.tar.gz
export PATH="/mnt/home/steepale/Apps/gffcompare-0.10.1.Linux_x86_64:${PATH}"
cd ${proj_dir}

### We will create our own transcriptome using a reference based approach via StringTie. This will find known and novel transcripts
# Each BAM file needs to be resorted by reference position because they were originally sorted by readname

# Iterate through each germline rnaseq BAM file and create a GTF file
for gbird in `cat ${HOME}/databases/samples/germline_samples_rnaseq_2014_017NNN.txt`
do
echo ${gbird}
stringtie \
${MDV_DIR}/rnaseq_somatic_snvs_indels/data/${gbird}_paired_norRNA_sorted_picard.bam \
-p 1 \
-G ${MDV_DIR}/RNA_DE/data/ref/ensemble/Gallus_gallus.Gallus_gallus-5.0.86.gtf \
-o ./data/${gbird}_stringtie_transcripts.gtf
done

# Iterate through each tumor rnaseq BAM files and create a GTF file
for tbird in `cat /mnt/home/steepale/databases/samples/tumor_samples_rnaseq_2014_017NNN-N_redunant_1s.txt`
do
echo ${tbird}
qsub -v Var=${tbird} -N "stringtie_known_novel_transcripts_"${tbird} ./scripts/stringtie_known_and_novel_transcripts.sh
done

# ./scripts/stringtie_known_and_novel_transcripts.sh
##########################################
#!/bin/bash -login
### Job name
### Resources
#PBS -l nodes=1:ppn=1,walltime=00:01:00:00,mem=4gb
### Send email if the job encounters an error
#PBS –m a
### Output files to where you submitted your batch file
#PBS -e ./jobs/$PBS_JOBNAME.err
#PBS -o ./jobs/$PBS_JOBNAME.log
#PBS -j oe

# Change to working directory
cd $PBS_O_WORKDIR

# Variables
stringtie="${HOME}/Apps/stringtie-1.3.3b.Linux_x86_64/stringtie"

# Perform stringtie analysis
stringtie \
${MDV_DIR}/rnaseq_somatic_snvs_indels/data/${Var}_paired_norRNA_sorted_picard.bam \
-p 1 \
-G ${MDV_DIR}/RNA_DE/data/ref/ensemble/Gallus_gallus.Gallus_gallus-5.0.86.gtf \
-o ./data/${Var}_stringtie_transcripts.gtf

# Collect stats on run
qstat -f ${PBS_JOBID}
##########################################

# Stringtie renames the genes with its own unique ID's. We need to capture those ID's which differ among samples
grep "IKZF1" ./data/*_stringtie_transcripts.gtf | \
grep -w "transcript" | cut -f9 | cut -d';' -f1 | \
sort | uniq | sed 's/"//g' |cut -d ' ' -f2 > \
./data/stringtie_ikzf1_gene_list.txt

# Change the stringtie gene id's given to Ikaros transcripts back to IKZF1 (controls)
for gbird in `cat ${HOME}/databases/samples/germline_samples_rnaseq_2014_017NNN.txt`
do
stringtie_id=`grep "IKZF1" ./data/${gbird}_stringtie_transcripts.gtf | grep -w "transcript" | cut -f9 | cut -d';' -f1 | sort | uniq | sed 's/"//g' |cut -d ' ' -f2`
cat ./data/${gbird}_stringtie_transcripts.gtf | sed "s/${stringtie_id}/IKZF1/g" | grep "IKZF1" > ./data/${gbird}_stringtie_transcripts_IKZF1_renamed.gtf
done

# Change the stringtie gene id's given to Ikaros transcripts back to IKZF1 (unknown mech cohort)
for tbird in `cat ${HOME}/databases/samples/IKZF1_unknown_mech_cohort_rna.txt`
do
echo ${tbird}
stringtie_id=`grep "IKZF1" ./data/${tbird}_stringtie_transcripts.gtf | grep -w "transcript" | cut -f9 | cut -d';' -f1 | sort | uniq | sed 's/"//g' |cut -d ' ' -f2`
cat ./data/${tbird}_stringtie_transcripts.gtf | sed "s/${stringtie_id}/IKZF1/g" | grep "IKZF1" > ./data/${tbird}_stringtie_transcripts_IKZF1_renamed.gtf
done

# Change the stringtie gene id's given to Ikaros transcripts back to IKZF1 (Mutated cohort)
for tbird in `cat ${HOME}/databases/samples/IKZF1_mutated_cohort_rna.txt`
do
echo ${tbird}
stringtie_id=`grep "IKZF1" ./data/${tbird}_stringtie_transcripts.gtf | grep -w "transcript" | cut -f9 | cut -d';' -f1 | sort | uniq | sed 's/"//g' |cut -d ' ' -f2`
cat ./data/${tbird}_stringtie_transcripts.gtf | sed "s/${stringtie_id}/IKZF1/g" | grep "IKZF1" > ./data/${tbird}_stringtie_transcripts_IKZF1_renamed.gtf
done

# Change the stringtie gene id's given to Ikaros transcripts back to IKZF1 (low IKZF1 cohort)
for tbird in `cat ${HOME}/databases/samples/IKZF1_low_expression_cohort_rna.txt`
do
echo ${tbird}
stringtie_id=`grep "IKZF1" ./data/${tbird}_stringtie_transcripts.gtf | grep -w "transcript" | cut -f9 | cut -d';' -f1 | sort | uniq | sed 's/"//g' |cut -d ' ' -f2`
cat ./data/${tbird}_stringtie_transcripts.gtf | sed "s/${stringtie_id}/IKZF1/g" | grep "IKZF1" > ./data/${tbird}_stringtie_transcripts_IKZF1_renamed.gtf
done


# Reformat the data for analysis
# Germline samples cohort
for gbird in `cat ${HOME}/databases/samples/germline_samples_rnaseq_2014_017NNN.txt`
do
python ./scripts/filter_stringtie_transcripts.py \
./data/${gbird}_stringtie_transcripts_IKZF1_renamed.gtf \
./data/${gbird}_stringtie_transcripts_IKZF1_filtered.gtf
echo ${gbird}
cat ./data/${gbird}_stringtie_transcripts_IKZF1_filtered.gtf | grep "Novel"
done

# Unknown Mechanism Cohort
for tbird in `cat ${HOME}/databases/samples/IKZF1_unknown_mech_cohort_rna.txt`
do
python ./scripts/filter_stringtie_transcripts.py \
./data/${tbird}_stringtie_transcripts_IKZF1_renamed.gtf \
./data/${tbird}_stringtie_transcripts_IKZF1_filtered.gtf
echo ${tbird}
cat ./data/${tbird}_stringtie_transcripts_IKZF1_filtered.gtf
done

# Low IKZF1 Expression Cohort
for tbird in `cat ${HOME}/databases/samples/IKZF1_low_expression_cohort_rna.txt`
do
python ./scripts/filter_stringtie_transcripts.py \
./data/${tbird}_stringtie_transcripts_IKZF1_renamed.gtf \
./data/${tbird}_stringtie_transcripts_IKZF1_filtered.gtf
echo ${tbird}
cat ./data/${tbird}_stringtie_transcripts_IKZF1_filtered.gtf | grep "Novel"
done

# Mutated Cohort
for tbird in `cat ${HOME}/databases/samples/IKZF1_mutated_cohort_rna.txt`
do
python ./scripts/filter_stringtie_transcripts.py \
./data/${tbird}_stringtie_transcripts_IKZF1_renamed.gtf \
./data/${tbird}_stringtie_transcripts_IKZF1_filtered.gtf
echo ${tbird}
cat ./data/${tbird}_stringtie_transcripts_IKZF1_filtered.gtf | grep "Novel"
done

# ./scripts/filter_stringtie_transcripts.py
############################
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
#############################

# Determine how many exons per transcript

# Count the number of known exons in galgal5 gtf file and create bed file for it
for transcript in "IKZF1-201" "IKZF1-202"
do
grep "IKZF1" ${MDV_DIR}/RNA_DE/data/ref/ensemble/Gallus_gallus.Gallus_gallus-5.0.86.gtf | \
grep -w "exon" | grep -w ${transcript} > ./data/${transcript}_exons.gtf
done

# Organize the transcripts into one file
grep -v "^CHROM" ./data/*_stringtie_transcripts_IKZF1_filtered.gtf | \
cut -d':' -f2 | cut -f1,2,3,4,5,9 | sort | uniq > ./data/IKZF1_transcripts_raw.txt

# manually adjust transcript names
vi ./data/IKZF1_transcripts_raw.txt

cat ./data/IKZF1_transcripts_raw.txt | sort -k3,3n
#2	StringTie	transcript_01	80915345	80986405	Known
#2	StringTie	transcript_02	80928032	80985307	Known
#2	StringTie	transcript_03	80915293	80974809	Novel
#2	StringTie	transcript_04	80915345	80986405	Novel
#2	StringTie	transcript_05	80915345	80989818	Novel
#2	StringTie	transcript_06	80915622	80986405	Novel
#2	StringTie	transcript_07	80915623	80985307	Novel
#2	StringTie	transcript_08	80915625	80986817	Novel
#2	StringTie	transcript_09	80915638	80985307	Novel
#2	StringTie	transcript_10	80915677	80986405	Novel
#2	StringTie	transcript_11	80915701	80986405	Novel
#2	StringTie	transcript_12	80915715	80986405	Novel
#2	StringTie	transcript_13	80915785	80977751	Novel
#2	StringTie	transcript_14	80915785	80986405	Novel
#2	StringTie	transcript_15	80915789	80986405	Novel
#2	StringTie	transcript_16	80915792	80986405	Novel
#2	StringTie	transcript_17	80915794	80985307	Novel
#2	StringTie	transcript_18	80915796	80986405	Novel
#2	StringTie	transcript_19	80915797	80979112	Novel
#2	StringTie	transcript_20	80915802	80985307	Novel
#2	StringTie	transcript_21	80915802	80986405	Novel
#2	StringTie	transcript_22	80915803	80986405	Novel
#2	StringTie	transcript_23	80915804	80985307	Novel
#2	StringTie	transcript_24	80915805	80985307	Novel
#2	StringTie	transcript_25	80915805	80986405	Novel
#2	StringTie	transcript_26	80915809	80986405	Novel
#2	StringTie	transcript_27	80915848	80950085	Novel
#2	StringTie	transcript_28	80915848	80986405	Novel
#2	StringTie	transcript_29	80916005	80974561	Novel
#2	StringTie	transcript_30	80916005	80986405	Novel
#2	StringTie	transcript_31	80916007	80986405	Novel
#2	StringTie	transcript_32	80916156	80986803	Novel
#2	StringTie	transcript_33	80917179	80917420	Novel
#2	StringTie	transcript_34	80917180	80917416	Novel
#2	StringTie	transcript_35	80927920	80989409	Novel
#2	StringTie	transcript_36	80928032	80954165	Novel
#2	StringTie	transcript_37	80928032	80963686	Novel
#2	StringTie	transcript_38	80928032	80977204	Novel
#2	StringTie	transcript_39	80928032	80979085	Novel
#2	StringTie	transcript_40	80928032	80985307	Novel
#2	StringTie	transcript_41	80928032	80986405	Novel
#2	StringTie	transcript_42	80928032	80986817	Novel
#2	StringTie	transcript_43	80928032	80988889	Novel
#2	StringTie	transcript_44	80934992	80971850	Novel
#2	StringTie	transcript_45	80964752	80985307	Novel
#2	StringTie	transcript_46	80965605	80971838	Novel
#2	StringTie	transcript_47	80965911	80986803	Novel
#2	StringTie	transcript_48	80973999	80986405	Novel
#2	StringTie	transcript_49	80974147	80986817	Novel
#2	StringTie	transcript_50	80974324	80987457	Novel
#2	StringTie	transcript_51	80980254	80986405	Novel
#2	StringTie	transcript_52	80980736	80986405	Novel
#2	StringTie	transcript_53	80982166	80986405	Novel
#2	StringTie	transcript_54	80984560	80989555	Novel
#2	StringTie	transcript_55	80987558	80989761	Novel
#2	StringTie	transcript_56	80987560	80990098	Novel
#2	StringTie	transcript_57	80987561	80988832	Novel
#2	StringTie	transcript_58	80987562	80989797	Novel
#2	StringTie	transcript_59	80987572	80989369	Novel

# Reformat files
for transcript in "IKZF1-201" "IKZF1-202"
do
python ./scripts/gtf_to_txt_exons.py \
./data/${transcript}_exons.gtf \
./data/${transcript}_exons.txt \
${transcript}
done

# ./scripts/gtf_to_txt_exons.py
################################
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
################################

# Determine the overlap manually
for transcript in "IKZF1-201" "IKZF1-202"
do
python ./scripts/count_exons_transcripts.py \
./data/IKZF1_transcripts_raw.txt \
./data/${transcript}_exons.txt \
./data/IKZF1_transcripts_${transcript}_exon_overlap.txt
done

# ./scripts/count_exons_transcripts.py
######################################
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
######################################

# Add the exons to this dataset
# Germline Samples cohort
for bird in `cat ${HOME}/databases/samples/germline_samples_rnaseq_2014_017NNN.txt | head -n1`
do
echo ${bird}
python ./scripts/add_exon_annotation.py \
./data/${bird}_stringtie_transcripts_IKZF1_filtered.gtf \
./data/${bird}_stringtie_transcripts_IKZF1_filtered_exons.txt
sort -k9,9 ./data/${bird}_stringtie_transcripts_IKZF1_filtered_exons.txt | column -t
done

# Unknown Mechanism Cohort
for bird in `cat ${HOME}/databases/samples/IKZF1_unknown_mech_cohort_rna.txt`
do
echo ${bird}
python ./scripts/add_exon_annotation.py \
./data/${bird}_stringtie_transcripts_IKZF1_filtered.gtf \
./data/${bird}_stringtie_transcripts_IKZF1_filtered_exons.txt
sort -k9,9 ./data/${bird}_stringtie_transcripts_IKZF1_filtered_exons.txt | column -t
done

# Low IKZF1 Expression Cohort
for bird in `cat ${HOME}/databases/samples/IKZF1_low_expression_cohort_rna.txt`
do
echo ${bird}
python ./scripts/add_exon_annotation.py \
./data/${bird}_stringtie_transcripts_IKZF1_filtered.gtf \
./data/${bird}_stringtie_transcripts_IKZF1_filtered_exons.txt
sort -k9,9 ./data/${bird}_stringtie_transcripts_IKZF1_filtered_exons.txt | column -t
done

# Mutated Cohort
for bird in `cat ${HOME}/databases/samples/IKZF1_mutated_cohort_rna.txt`
do
echo ${bird}
python ./scripts/add_exon_annotation.py \
./data/${bird}_stringtie_transcripts_IKZF1_filtered.gtf \
./data/${bird}_stringtie_transcripts_IKZF1_filtered_exons.txt
sort -k9,9 ./data/${bird}_stringtie_transcripts_IKZF1_filtered_exons.txt | column -t
done


# ./scripts/add_exon_annotation.py
################################
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
##################################

# Based on the results from STRINGTIE, I do not see any significant evidence of dominant negative IKZF1 transcript isoforms

# Run a preliminary analysis with gffcompare
#for gbird in `cat ${HOME}/databases/samples/germline_samples_rnaseq_2014_017NNN.txt`
#do
#gffcompare -R \
#-r ${MDV_DIR}/RNA_DE/data/ref/ensemble/Gallus_gallus.Gallus_gallus-5.0.86.gtf \
#-o ./data/${gbird}_stringtie_transcripts_IKZF1_renamed \
#./data/${gbird}_stringtie_transcripts_IKZF1_renamed.gtf
#done

# We will begin by using Salmon to quantify transcripts across the transcriptome, but right now our primary focus is IKZF1 isoforms

### Obtain a transcriptome
# We will test with ensembl's default transcriptome first
cd ./data/
wget -r -np -nH ftp://ftp.ensembl.org/pub/release-89/fasta/gallus_gallus/cdna/ 
mv ${proj_dir}/data/pub/release-89/fasta/gallus_gallus/cdna/* ./
rm -r pub
cd ..

# Build an index of the transcriptome with Salmon, although this may not be necessary
salmon index \
-t ./data/Gallus_gallus.Gallus_gallus-5.0.cdna.all.fa.gz \
-i ./data/Gallus_gallus.Gallus_gallus-5.0.cdna.all_index

# Quantify the transcrips (only known reference transcripts; IKZF1-201, -202)
# Germline Samples cohort
for bird in `cat ${HOME}/databases/samples/germline_samples_rnaseq_2014_017NNN.txt`
do
echo ${bird}
# Quantify known transcripts with Salmon
salmon quant \
-i ./data/Gallus_gallus.Gallus_gallus-5.0.cdna.all_index \
-l A \
-1 ${MDV_DIR}/RNA_DE/data/${bird}/${bird}_R1_paired_norRNA.fastq.gz \
-2 ${MDV_DIR}/RNA_DE/data/${bird}/${bird}_R2_paired_norRNA.fastq.gz \
-p 4 \
-o ./data/salmon_quants/${bird}_quant
done

# Unknown Mechanism Cohort
for bird in `cat ${HOME}/databases/samples/IKZF1_unknown_mech_cohort_rna.txt`
do
echo ${bird}
# Quantify known transcripts with Salmon
salmon quant \
-i ./data/Gallus_gallus.Gallus_gallus-5.0.cdna.all_index \
-l A \
-1 ${MDV_DIR}/RNA_DE/data/${bird}/${bird}_R1_paired_norRNA.fastq.gz \
-2 ${MDV_DIR}/RNA_DE/data/${bird}/${bird}_R2_paired_norRNA.fastq.gz \
-p 4 \
-o ./data/salmon_quants/${bird}_quant
done

# Low IKZF1 Expression Cohort
for bird in `cat ${HOME}/databases/samples/IKZF1_low_expression_cohort_rna.txt`
do
echo ${bird}
# Quantify known transcripts with Salmon
salmon quant \
-i ./data/Gallus_gallus.Gallus_gallus-5.0.cdna.all_index \
-l A \
-1 ${MDV_DIR}/RNA_DE/data/${bird}/${bird}_R1_paired_norRNA.fastq.gz \
-2 ${MDV_DIR}/RNA_DE/data/${bird}/${bird}_R2_paired_norRNA.fastq.gz \
-p 4 \
-o ./data/salmon_quants/${bird}_quant
done

# Mutated Cohort
for bird in `cat ${HOME}/databases/samples/IKZF1_mutated_cohort_rna.txt`
do
echo ${bird}
# Quantify known transcripts with Salmon
salmon quant \
-i ./data/Gallus_gallus.Gallus_gallus-5.0.cdna.all_index \
-l A \
-1 ${MDV_DIR}/RNA_DE/data/${bird}/${bird}_R1_paired_norRNA.fastq.gz \
-2 ${MDV_DIR}/RNA_DE/data/${bird}/${bird}_R2_paired_norRNA.fastq.gz \
-p 4 \
-o ./data/salmon_quants/${bird}_quant
done

### Print the counts
# Germline Samples cohort
echo "Germline Samples cohort"
for bird in `cat ${HOME}/databases/samples/germline_samples_rnaseq_2014_017NNN.txt`
do
echo ${bird}
grep -e "ENSGALT00000021364.5" -e "ENSGALT00000081310.1" ./data/salmon_quants/${bird}_quant/quant.sf
done
# Unknown Mechanism Cohort
echo "Unknown Mechanism Cohort"
for bird in `cat ${HOME}/databases/samples/IKZF1_unknown_mech_cohort_rna.txt`
do
echo ${bird}
grep -e "ENSGALT00000021364.5" -e "ENSGALT00000081310.1" ./data/salmon_quants/${bird}_quant/quant.sf
done
# Low IKZF1 Expression Cohort
echo "Low IKZF1 Expression Cohort"
for bird in `cat ${HOME}/databases/samples/IKZF1_low_expression_cohort_rna.txt`
do
echo ${bird}
grep -e "ENSGALT00000021364.5" -e "ENSGALT00000081310.1" ./data/salmon_quants/${bird}_quant/quant.sf
done
# Mutated Cohort
echo "Mutated Cohort"
for bird in `cat ${HOME}/databases/samples/IKZF1_mutated_cohort_rna.txt`
do
echo ${bird}
grep -e "ENSGALT00000021364.5" -e "ENSGALT00000081310.1" ./data/salmon_quants/${bird}_quant/quant.sf
done

# Some of the samples in the Unknown Mechanism Cohort show reduced levels of normal IKZF1 transcripts.
# No clear pattern yet

# Will examine IKZF1 mutations




















