#!/bin/bash


#SBATCH --job-name=with_chr_hg002_xtea_testing
#SBATCH -A park_contrib
#SBATCH --partition park
#SBATCH --mem 2G
#SBATCH -c 16 
#SBATCH -t 00:30:00
#SBATCH -o slurm_out/slurm-%x.%j.out



#bam_file=SM-N8UQR.cram
#sample_id=SM-N8UQR

bam_file=hg002_data/sarek/HG002_60x_hg38.cram
sample_id=HG002_60x 

#bam_file=small_test.bam
#sample_id=SMALL_TEST

ref=~/data1/SOFTWARE/REFERENCE/hg38_no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
exclude_regions=~/data1/corinne/ref/encode_exclude_regions_ENCFF356LFX.bed
rmsk=~/data1/corinne/ref/rmsk/all_active_tes.bed

output_dir=${sample_id}

if [ ! -d ${output_dir} ]; then
	mkdir ${output_dir}
fi
#
### get SV reads https://github.com/arq5x/lumpy-sv/blob/master/src/filter/filter.c
### ~10 minutes on 80G cram
#echo "Getting clipped regions..."
#time get_clip_regions -f ${ref} -e ${exclude_regions} -r ${rmsk} \
#	${bam_file}\
#	${output_dir}/${sample_id}_clip_positions.tsv \
#	16	
#echo "finished!"
#echo
##

# took ~6 minutes on 80G cram (60x)
#echo "Calculating local depth..."
#tail -n +2 ${output_dir}/${sample_id}_clip_positions.tsv | \
#	bedtools merge -c 4 -o distinct -d 5000 -i - > ${output_dir}/${sample_id}_clip_positions_5kb.bed
#
#awk '{print $1"\t"$2"\t"$3"\tid_"NR}' ${output_dir}/${sample_id}_clip_positions_5kb.bed > ${output_dir}/${sample_id}_clip_positions_5kb.tsv
#paste <(cut -f4 ${output_dir}/${sample_id}_clip_positions_5kb.tsv) <(cut -f4 ${output_dir}/${sample_id}_clip_positions_5kb.bed) > ${sample_id}/${sample_id}_5kb_regions_key.tsv
#
#./src/2_calculate_local_depth.sh ${output_dir}/${sample_id}_clip_positions_5kb.tsv \
#	${bam_file} \
#	${output_dir}/${sample_id}_local_depth.tsv 16 
#echo "finished!"
#echo


# remove any calls with < 10% clipped read support based on local depth
echo "Filtering read based on depth..."
./src/3_filter_based_on_local_depth.py -l ${output_dir}/${sample_id}_local_depth.tsv \
	-k ${output_dir}/${sample_id}_5kb_regions_key.tsv  \
	-c ${output_dir}/${sample_id}_clip_positions.tsv \
	-o ${output_dir}/${sample_id}_clip_filtered_candidates.tsv
echo "finished!"
echo

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
