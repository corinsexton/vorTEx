#!/bin/bash

#SBATCH --job-name=gridss
#SBATCH -A park
#SBATCH --partition short
#SBATCH --mem 24G
#SBATCH -c 8
#SBATCH -t 02:00:00
#SBATCH -o slurm-%x.%j.out


module load gcc/9.2.0 R/4.2.1

# ended up having to hand compile gridsstools with gcc 9.2
# otherwise libcurl version issues


bed=/home/cos689/dev_xTea/HG002_60x/HG002_60x_clip_filtered_candidates.tsv
bam=/home/cos689/dev_xTea/hg002_data/sarek/HG002_60x_hg38.cram
output_vcf=test.vcf

time /home/cos689/software/bin/gridss_extract_overlapping_fragments \
  --targetbed ${bed} \
  --targetmargin 500 \
  -j /home/cos689/software/packages/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
  -o ${bam##*/}.targeted.bam \
  ${bam}


#time echo "RUNNING GRIDSS..."
#gridss \
#	-r /n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/hg38/cgap_matches/Homo_sapiens_assembly38.fa \
#	-j /home/cos689/software/packages/gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \
#	-o ${output_vcf} \
#	-a ${bam}.assembly.bam \
#	${bam##*/}.targeted.bam

#samtools consensus -f fasta in.bam -o cons.fa 

