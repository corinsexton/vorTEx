#!/bin/bash


#SBATCH --job-name=xtea_testing
#SBATCH -A park_contrib
#SBATCH --partition park
#SBATCH --mem 16G
#SBATCH -c 8 
#SBATCH -t 01:00:00
#SBATCH -o slurm-%x.%j.out



bam_file=SM-N8UQR.cram
sample_id=SM-N8UQR


ref=~/data1/SOFTWARE/REFERENCE/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

output_dir=${sample_id}

#if [ ! -d ${output_dir} ]; then
#	mkdir ${output_dir}
#fi
#
## get SV reads https://github.com/arq5x/lumpy-sv/blob/master/src/filter/filter.c
## ~20 minutes on 16G cram
#echo "Running lumpy-filter..."
#time lumpy_filter -f ${ref} \
#	${bam_file}\
#	${output_dir}/${sample_id}_clip.bam ${output_dir}/${sample_id}_disc.bam \
#	8
#echo "finished!"
#echo
#
## turn clip and disc to fasta
## <1 minute
#echo "Creating fastas..."
#samtools fasta ${output_dir}/${sample_id}_clip.bam > ${output_dir}/${sample_id}_clip.fa
#samtools fasta ${output_dir}/${sample_id}_disc.bam > ${output_dir}/${sample_id}_disc.fa
#echo "finished!"
#echo

# look for TE sequence specific reads
# ~1 minute 16G memory
# memory needed
echo "Running mmseqs"
time mmseqs easy-search -s 1 --max-seqs 1 --search-type 3 \
	--threads 8 \
	${output_dir}/${sample_id}_clip.fa \
	${output_dir}/${sample_id}_disc.fa \
	mmseqs/mmseq_dbs/l1hs_mmseq_db \
	${output_dir}/${sample_id}_TE_mmseq.output \
	${output_dir}/${sample_id}_mmseq_tmp
echo "finished!"
echo
