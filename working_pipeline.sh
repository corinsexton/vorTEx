#!/bin/bash


#SBATCH --job-name=with_chr_hg002_xtea_testing
#SBATCH -A park_contrib
#SBATCH --partition park
#SBATCH --mem 2G
#SBATCH -c 16 
#SBATCH -t 00:30:00
#SBATCH -o slurm-%x.%j.out



#bam_file=SM-N8UQR.cram
#sample_id=SM-N8UQR

bam_file=hg002_data/sarek/HG002_60x_hg38.cram
sample_id=HG002_60x 

#bam_file=small_test.bam
#sample_id=SMALL_TEST

ref=~/data1/SOFTWARE/REFERENCE/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
exclude_regions=~/data1/corinne/ref/encode_exclude_regions_ENCFF356LFX.bed

output_dir=${sample_id}

if [ ! -d ${output_dir} ]; then
	mkdir ${output_dir}
fi
#
## get SV reads https://github.com/arq5x/lumpy-sv/blob/master/src/filter/filter.c
## ~20 minutes on 16G cram
#echo "Running lumpy-filter..."
time get_clip_regions -f ${ref} -e ${exclude_regions} \
	${bam_file}\
	${output_dir}/${sample_id}_clip_positions.tsv \
	16	
#echo "finished!"
#echo
#
### turn clip and disc to fasta
### <1 minute
#echo "Creating fastas..."
##samtools fasta ${output_dir}/${sample_id}_clip.bam > ${output_dir}/${sample_id}_clip.fa
#samtools fasta ${output_dir}/${sample_id}_disc.bam > ${output_dir}/${sample_id}_disc.fa
#echo "finished!"
##echo

# look for TE sequence specific reads
# ~1 minute 16G memory
# memory needed
#echo "Running mmseqs"
#time mmseqs easy-search -s 1 --max-seqs 1 --search-type 3 \
#	--threads 8 \
#	${output_dir}/${sample_id}_softclip.fa \
#	mmseqs/alu_l1_sva/alu_l1_sva_db \
#	${output_dir}/${sample_id}_TE_mmseq.output \
#	${output_dir}/${sample_id}_mmseq_tmp
#echo "finished!"
#echo
