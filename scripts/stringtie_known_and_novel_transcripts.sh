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
