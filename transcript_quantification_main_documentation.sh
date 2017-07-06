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

# PROJECT DIRECTORY (TUM Cluster)
proj_dir="/mnt/research/ADOL/OutsideCollaborations/20160201_Cheng_Steep_Xu_Zhang/transcript_quantification"
cd $proj_dir

# Make proper directories
mkdir -p ./{data,scripts,analysis,jobs}
mkdir ./salmon

# Sources to cite:
# Salmon main page: https://combine-lab.github.io/salmon/getting_started/
# Salmon GitHub: https://github.com/COMBINE-lab/salmon

# We will begin by using Salmon to quantify transcripts across the transcriptome, but right now our primary focus is IKZF1 isoforms

# Install Salmon v0.8.2
cd ${HOME}/Apps
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.8.2/Salmon-0.8.2_linux_x86_64.tar.gz
tar xzvf Salmon-0.8.2_linux_x86_64.tar.gz
rm Salmon-0.8.2_linux_x86_64.tar.gz
# Add salmon to PATH (in ${HOME}/.bashrc)
export PATH="/mnt/home/steepale/Apps/Salmon-0.8.2_linux_x86_64:${PATH}"

# Obtain a transcriptome
cd ./data/
wget -r -np -nH ftp://ftp.ensembl.org/pub/release-89/fasta/gallus_gallus/cdna/ 
mv ${proj_dir}/data/pub/release-89/fasta/gallus_gallus/cdna/* ./
rm -r pub
cd ..

# Build an index of the transcriptome
salmon index \
-t ./data/Gallus_gallus.Gallus_gallus-5.0.cdna.all.fa.gz \
-i ./data/Gallus_gallus.Gallus_gallus-5.0.cdna.all_index








