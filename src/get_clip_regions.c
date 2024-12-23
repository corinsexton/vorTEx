#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <err.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>
#include <assert.h>


#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <htslib/hts.h>
#include <htslib/khash.h>

#define MAX_LINE_LENGTH 1024

#define MIN_MAPQ 10
#define MIN_CLIP_LEN 9

#define MAX_CHROMOSOMES 24

// Define a hash map to store the intervals blacklist regions
KHASH_MAP_INIT_STR(ranges, int*)

// Define a hashmap for storing the count of clipped reads at each position
KHASH_MAP_INIT_INT64(uint64_t, int)  // A simple hash map with int keys (positions) and int values (counts)

// Example chromosome names (you can expand this list)
const char *chromosome_names[MAX_CHROMOSOMES] = {
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", 
    "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", 
    "chrX", "chrY"
};

// Function to get a chromosome ID based on the name
int get_chromosome_id(const char *chromosome) {
	int i;
    for (i = 0; i < MAX_CHROMOSOMES; i++) {
        if (strcmp(chromosome_names[i], chromosome) == 0) {
            return i;  // Return the index as the unique ID
        }
    }
    return -1;  // Return -1 if chromosome is not found
}

// Function to encode chromosome and position into a unique integer key
uint64_t encode_position(int chromosome_id, int position) {
    // Assuming that chromosome_id can fit in 32 bits and position in 32 bits
    return ((uint64_t)chromosome_id << 32) | (position & 0xFFFFFFFF);
}

// Function to decode the unique integer key back into chromosome ID and position
void decode_position(uint64_t key, int *chromosome_id, int *position) {
    *chromosome_id = key >> 32;
    *position = (int)(key & 0xFFFFFFFF);  // Lower 32 bits are the position
}

//// Function to update the position in the hash map
//void update_position(khash_t(pos) *hash, const char *chromosome, int pos, int start_count, int end_count) {
//    // Get the chromosome ID
//    int chromosome_id = get_chromosome_id(const char *chrom= sam_hdr_tid2name(header, read->core.tid););
//    if (chromosome_id == -1) {
//        fprintf(stderr, "Error: Chromosome not found!\n");
//        return;
//    }
//
//    // Encode the chromosome and position into a unique integer key
//    uint64_t key = encode_position(chromosome_id, pos);
//
//    // Insert or update the hash map
//    khint_t k = kh_get(pos, hash, key);
//    if (k == kh_end(hash)) {
//        // Position not found, create a new entry
//        k = kh_put(pos, hash, key, &found);
//        kh_value(hash, k) = start_count + end_count;  // Initialize with counts
//    } else {
//        // Update existing position
//        kh_value(hash, k) += start_count + end_count;
//    }
//}




// Define a structure for storing RepeatMasker regions (e.g., Alu, SVA, L1)
typedef struct {
    char *chr;
    int start;
    int end;
    char *repeat_class;
} RepeatRegion;

// Structure to store position and its corresponding start and end clipped counts
typedef struct {
	char *chromosome;  // Chromosome name (not used directly in hash map)
    int pos;
    int start_count;
    int end_count;
} PositionClipped;



// Comparator function for sorting RepeatRegions by chromosome and start position
int compare_repeat_regions(const void *a, const void *b) {
    RepeatRegion *r1 = (RepeatRegion *)a;
    RepeatRegion *r2 = (RepeatRegion *)b;

    // First compare chromosome names
    int cmp_chr = strcmp(r1->chr, r2->chr);
    if (cmp_chr != 0) return cmp_chr;

    // If chromosomes are the same, compare start positions
    return r1->start - r2->start;
}

// Function to load RepeatMasker regions from a BED file into a sorted array
int load_repeat_masker_bed(const char *bed_file, RepeatRegion **repeat_regions) {
    FILE *fp = fopen(bed_file, "r");
    if (!fp) {
        perror("Error opening RepeatMasker BED file");
        return -1;
    }

    char line[1024];
    int count = 0;
    RepeatRegion *regions = malloc(sizeof(RepeatRegion) * 1000);  // Initially allocate space for 1000 regions

    while (fgets(line, sizeof(line), fp)) {
        // Skip comment lines
        if (line[0] == '#') continue;

        char chr[100];
        int start, end;
        char repeat_class[100];

        // Parse the BED file format (chr, start, end, repeat_class)
        if (sscanf(line, "%s %d %d %s", chr, &start, &end, repeat_class) != 4) {
            continue;
        }

        // Store the repeat region
        regions[count].chr = strdup(chr);
        regions[count].start = start;
        regions[count].end = end;
        regions[count].repeat_class = strdup(repeat_class);
        count++;

        // Resize memory if needed
        if (count % 1000 == 0) {
            regions = realloc(regions, sizeof(RepeatRegion) * (count + 1000));
        }
    }

    fclose(fp);

    // Sort repeat regions by chromosome and start position (for fast lookup)
    qsort(regions, count, sizeof(RepeatRegion), compare_repeat_regions);

    *repeat_regions = regions;
    return count;
}

// Function to check if a mate is in a repeat region (using binary search)
int is_mate_in_repeat(bam1_t *read, bam_hdr_t *header, RepeatRegion *repeat_regions, int num_regions) {
    if (read->core.mtid == -1) return 0; // No mate

    // Perform a binary search for the correct repeat region
	int i;
    for (i = 0; i < num_regions; i++) {
        if (strcmp(repeat_regions[i].chr, header->target_name[read->core.mtid]) == 0) {
            if (read->core.mpos >= repeat_regions[i].start && 
                read->core.mpos <= repeat_regions[i].end) {
                return 1; // Mate is in a repeat region
            }
        }
    }

    return 0; // Mate is not in a repeat region
}



char int_to_base(int base_int) {
    switch (base_int) {
        case 1: return 'A';
        case 2: return 'C';
        case 4: return 'G';
        case 8: return 'T';
        case 15: return 'T';
        default: return '?';  // 'N' for unknown or invalid base
    }
}



// Function to load regions from the file into a hash map
void load_regions(const char *filename, khash_t(ranges) *region_map) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening region file");
        exit(1);
    }

    char line[MAX_LINE_LENGTH];
    while (fgets(line, sizeof(line), file)) {
        char chrom[256];
        int start, stop;
        if (sscanf(line, "%255s\t%d\t%d", chrom, &start, &stop) != 3) {
            fprintf(stderr, "Error parsing line: %s", line);
            continue;
        }

        // Add the region to the map (chromosome -> list of intervals)
        khiter_t k = kh_get(ranges, region_map, chrom);
		for (k = 0; k < kh_size(region_map); k++) {
		    if (k == kh_end(region_map)) {
            	// Chromosome not found, create a new list
            	int *intervals = malloc(2 * sizeof(int));
            	intervals[0] = start;
            	intervals[1] = stop;
            	kh_put(ranges, region_map, strdup(chrom), &k);
            	kh_value(region_map, k) = intervals;
        	} else {
        	    // Chromosome found, extend the list of intervals
        	    int *intervals = kh_value(region_map, k);
        	    int len = intervals[0] + 1; // Calculate the current length of the interval array
        	    intervals = realloc(intervals, (len + 2) * sizeof(int)); // Resize the array
        	    intervals[len] = start;
        	    intervals[len + 1] = stop;
        	    kh_value(region_map, k) = intervals;
        	}
		}
    }

    fclose(file);
}


void process_cigar(bam1_t *b, char *seq, bam_hdr_t *header, int *clip_side, khash_t(uint64_t) *start_clipped_map, khash_t(uint64_t) *end_clipped_map) {
    uint32_t *cigar = bam_get_cigar(b);  // Get CIGAR operations for this read
    int seq_len = b->core.l_qseq;  // Length of the read sequence
    int b_pos = b->core.pos;
	const char *chrom= sam_hdr_tid2name(header, b->core.tid);
    int chromosome_id = get_chromosome_id(chrom);
	if (chromosome_id == -1) { return; };

    *clip_side = 0;  // Initialize to no clipping

    int pos = 0;  // Position in the sequence
    int i;
    
    // Check if there's a soft clip at the start (left side)
    if ((cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
        *clip_side |= 1;  // Mark left-side soft clip

		uint64_t key = encode_position(chromosome_id, b_pos);

		// Check if start_pos already exists, otherwise add it to the hash map
        int ret;
        khiter_t k = kh_get(uint64_t, start_clipped_map, key);
        if (k == kh_end(start_clipped_map)) {  // Position doesn't exist yet
            k = kh_put(uint64_t, start_clipped_map, key, &ret);
            kh_value(start_clipped_map, k) = 1;
        } else {
            kh_value(start_clipped_map, k)++;
        }

        //int clip_len = cigar[0] >> BAM_CIGAR_SHIFT;  // Get the length of the left soft-clip
		//int j;
        //for (j = 0; j < clip_len; j++) {
        //    seq[pos++] = int_to_base(bam_seqi(bam_get_seq(b), j)); 
        //}
    }

    // Check if there's a soft clip at the end (right side)
    if ((cigar[b->core.n_cigar - 1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
        *clip_side |= 2;  // Mark right-side soft clip

        int last_matched_position = b_pos + b->core.l_qseq;
		uint64_t key = encode_position(chromosome_id, last_matched_position);

		int ret;
        khiter_t k = kh_get(uint64_t, end_clipped_map, key);
        if (k == kh_end(end_clipped_map)) {  // Position doesn't exist yet
            k = kh_put(uint64_t, end_clipped_map, key, &ret);
            kh_value(end_clipped_map, k) = 1;
        } else {
            kh_value(end_clipped_map, k)++;
        }

        //int clip_len = cigar[b->core.n_cigar - 1] >> BAM_CIGAR_SHIFT;  // Get the length of the right soft-clip
        //int start_pos = seq_len - clip_len;  // Start position for right-side soft clip
		//int j;
        //for (j = start_pos; j < seq_len; j++) {
        //    seq[pos++] = int_to_base(bam_seqi(bam_get_seq(b), j)); 
        //}
    }
    //seq[pos] = '\0';  // Null-terminate the sequence string

}

// Function to collect positions and counts from both hash maps into a list
void collect_positions(khash_t(uint64_t) *start_clipped_map, khash_t(uint64_t) *end_clipped_map, PositionClipped **positions, int *num_positions) {
    // Start with the positions from the start_clipped_map
	khiter_t k;
    for (k = kh_begin(start_clipped_map); k != kh_end(start_clipped_map); ++k) {
        if (kh_exist(start_clipped_map, k)) {
            int start_pos = kh_key(start_clipped_map, k);
            int start_count = kh_value(start_clipped_map, k);

            // Check if the position is also in the end_clipped_map
            khiter_t k_end = kh_get(uint64_t, end_clipped_map, start_pos);
            int end_count = (k_end != kh_end(end_clipped_map)) ? kh_value(end_clipped_map, k_end) : 0;

            // Add the position to the list
            (*positions) = realloc(*positions, (*num_positions + 1) * sizeof(PositionClipped));
            (*positions)[*num_positions].pos = start_pos;
            (*positions)[*num_positions].start_count = start_count;
            (*positions)[*num_positions].end_count = end_count;
            (*num_positions)++;
        }
    }

    // Now check for positions that exist only in the end_clipped_map
    for (k = kh_begin(end_clipped_map); k != kh_end(end_clipped_map); ++k) {
        if (kh_exist(end_clipped_map, k)) {
            int end_pos = kh_key(end_clipped_map, k);
            int end_count = kh_value(end_clipped_map, k);

            // Check if the position is already in the start_clipped_map
            khiter_t k_start = kh_get(uint64_t, start_clipped_map, end_pos);
            if (k_start == kh_end(start_clipped_map)) {
                // Position exists only in the end_clipped_map
                (*positions) = realloc(*positions, (*num_positions + 1) * sizeof(PositionClipped));
                (*positions)[*num_positions].pos = end_pos;
                (*positions)[*num_positions].start_count = 0;
                (*positions)[*num_positions].end_count = end_count;
                (*num_positions)++;
            }
        }
    }
}

// Comparison function for sorting the positions by their position values
int compare_positions(const void *a, const void *b) {
    PositionClipped *pos_a = (PositionClipped *)a;
    PositionClipped *pos_b = (PositionClipped *)b;
    return pos_a->pos - pos_b->pos;  // Sort by position
}

void merge_positions(PositionClipped *positions, int *num_positions) {
    if (*num_positions <= 1) return;  // If there's 0 or 1 position, no merge is needed.

    int write_index = 1;  // Start from the second position
	int i;
    for (i = 1; i < *num_positions; i++) {
        // Check if current position is within 50 bp of the previous one
        if (positions[i].pos - positions[write_index - 1].pos <= 50) {
            // Merge counts
            positions[write_index - 1].start_count += positions[i].start_count;
            positions[write_index - 1].end_count += positions[i].end_count;
        } else {
            // Move the current position to the next available spot
            positions[write_index] = positions[i];
            write_index++;
        }
    }

    // Update the number of positions after merging
    *num_positions = write_index;
}



// Function to print the sorted list of positions with their start and end counts
void print_sorted_positions(PositionClipped *positions, int num_positions, const char *filename) {

	FILE *output_file = fopen(filename, "w");
    if (!output_file) {
        perror("Error opening output file");
        return;
    }

    // Print the header line
    fprintf(output_file, "Position\tStart_Clipped\tEnd_Clipped\n");

	int chrom_id, position;

	int i;
    for (i = 0; i < num_positions; i++) {
		if (positions[i].start_count + positions[i].end_count >= 2) {
			decode_position(positions[i].pos,&chrom_id, &position);
            fprintf(output_file, "%s:%d\t%d\t%d\n",chromosome_names[chrom_id], position, positions[i].start_count, positions[i].end_count);
        }
    }
    fclose(output_file);

}

// Function to check if a read overlaps with any region
int overlaps(const bam1_t *read, bam_hdr_t *header,khash_t(ranges) *region_map) {
    // Get the chromosome name of the read
	const char *chrom= sam_hdr_tid2name(header, read->core.tid);

    // Find the corresponding list of regions for the chromosome
    khiter_t k = kh_get(ranges, region_map, chrom);
    if (k == kh_end(region_map)) {
        return 0; // No regions for this chromosome, include the read
    }

    int *regions = kh_value(region_map, k);

    // Check if the read overlaps with any region
	int i;
    for (i = 0; i < regions[0]; i += 2) {
        int start = regions[i];
        int stop = regions[i + 1];
        if (read->core.pos < stop && (read->core.pos + read->core.l_qseq) > start) {
            return 1; // Read overlaps with this region
        }
    }

    return 0; // Read does not overlap with any region
}

void rename_qname(bam1_t *b, bam_hdr_t *header) {
    // Get chromosome name from reference (chromosome is based on reference id)
	const char *chrom= sam_hdr_tid2name(header, b->core.tid);

	const char *mate_ref_name = "none";
	int mate_pos = 0;

	if (b->core.mtid >= 0) {
        // Mate exists (mtid >= 0 indicates the mate's chromosome is set)
        mate_ref_name = sam_hdr_tid2name(header, b->core.mtid);
        mate_pos = b->core.mpos;
    } 
    
    // Get position from the bam1_t record
    int position = b->core.pos + 1; // BAM positions are 0-based, so we add 1 to get 1-based position
    
    // Get CIGAR string (it's an array of operations)
    char cigar_str[1024];
    int cigar_len = 0;
    int i;  // Declare the loop index variable outside the loop
    for (i = 0; i < b->core.n_cigar; i++) {
        uint32_t op = bam_get_cigar(b)[i] & BAM_CIGAR_MASK;  // Corrected: Access b->core.cigar
        uint32_t len = bam_get_cigar(b)[i] >> BAM_CIGAR_SHIFT;
        // Append the CIGAR operation to the string
        cigar_len += snprintf(cigar_str + cigar_len, sizeof(cigar_str) - cigar_len, "%d%c", len, bam_cigar_opchr(op));
    }
    
    // Create a new qname string (CIGAR + chromosome + position)
    char new_qname[1024];
    snprintf(new_qname, sizeof(new_qname), "%s::%s::%d::%s::%d::%s", cigar_str, chrom, position,mate_ref_name,mate_pos,bam_get_qname(b));
    
	bam_set_qname(b,new_qname);
}






int main(int argc, char **argv)
{
    // without these 2 lines, htslib sometimes tries to download a part of the sequence
    // even though the -f reference was provided.
    setenv("REF_CACHE", "", 0);
    setenv("REF_PATH", "fake_value_so_no_download", 0);
    if (argc < 5)
        errx(1,
                "usage\t:%s -f <optional-reference> -e <exclude_bed_file> <bam> <clip_position.tsv> (optional #threads)",
                argv[0]);


    int aligned_reads = 0; // track total aligned reads. useful for QC without 2nd pass over full cram.

	// Initialize hash maps for start and end clipped positions
    khash_t(uint64_t) *start_clipped_map = kh_init(uint64_t);
    khash_t(uint64_t) *end_clipped_map = kh_init(uint64_t);

    int opt;
    char *fasta = NULL;
    char *exclude_regions_file_name;
    while ((opt = getopt(argc, argv, "f:e:")) != -1) {
        switch (opt) {
            case 'f': // Input file
                fasta = optarg;
                break;
            case 'e': // Another option
                exclude_regions_file_name = optarg;
                break;
        }
    }


    char *bam_file_name = argv[optind];
    //char *fasta_file_name= argv[1+optind];
	char *output_file_name= argv[1+optind];

    int threads = 2;
    if (argc == 3+optind) {
        threads = atoi(argv[2+optind]);
    }

 	// Initialize the region map (chromosome -> list of intervals)
    khash_t(ranges) *region_map = kh_init(ranges);

    // Load the regions into the map
    load_regions(exclude_regions_file_name, region_map);



    samFile *in = sam_open(bam_file_name, "rb");
    if(in == NULL) {
        errx(1, "Unable to open BAM/SAM file.");
    }
    // 0x1 | 0x2 | 0x4 | 0x8 | 0x10 | 0x20 | 0x40 | 0x80 | 0x100 | 0x200 | 0x800
    // decode everything but base-qual
    hts_set_opt(in, CRAM_OPT_REQUIRED_FIELDS, 3071);
    hts_set_opt(in, CRAM_OPT_DECODE_MD, 1);
    if (fasta != NULL) {
        hts_set_fai_filename(in, fasta);
    }
    if (threads > 1) {
        hts_set_threads(in, threads);
    }

    bam_hdr_t *hdr = sam_hdr_read(in);


	//// Open the FASTA output file
	//FILE *out = fopen(fasta_file_name, "w");
	//if (out == NULL) {
	//    fprintf(stderr, "Error opening output file %s\n", fasta_file_name);
	//}

    bam1_t *aln = bam_init1();
    int ret;

    while(ret = sam_read1(in, hdr, aln) >= 0) {
		// unmapped, qcfail, duplicate == skip, low mapping qual
        if ((((aln->core.flag) & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP)) != 0) || (aln->core.qual < MIN_MAPQ))
            continue;

		// overlaps with excluded regions
		if (overlaps(aln,hdr, region_map)) {
            continue;
        }

        aligned_reads += 1;


		// CLIPPED READS
    	int n_cigar = aln->core.n_cigar;
    	uint32_t *cigar = bam_get_cigar(aln);

		char seq[1024];  // Buffer to hold soft-clipped sequence (adjust size if necessary)
        int clip_side = 0;  // 1 = left-side soft clip, 2 = right-side soft clip, 3 = both sides

 		// Process the read's CIGAR string and extract soft-clipped bases
        process_cigar(aln, seq, hdr, &clip_side, start_clipped_map,end_clipped_map);

        //// Skip reads with clipping on both sides
        //if (clip_side == 3) {
        //    continue;
        //}

        //// If there's soft clipping (on either side), write it to the FASTA file
        //if (seq[0] != '\0') {  // Only output if there was any soft-clipping
		//	rename_qname(aln,hdr);
        //    fprintf(out, ">%s\n%s\n", bam_get_qname(aln), seq);  // Write to FASTA file
        //    continue;
        //}

		//// DISCORDANTS READS
		////read is not mapped in proper pair (large insert size), the mate IS (not un) mapped and NOT a suppl aln
        //if ((((aln->core.flag) & (BAM_FPROPER_PAIR | BAM_FMUNMAP | BAM_FSUPPLEMENTARY)) == 0) && 
		//		(aln->core.flag & BAM_FPAIRED == 1)) {
		//	//if (is_mate_in_repeat(read, repeat_regions, num_regions)) { }
		//	rename_qname(aln,hdr);
        //    r = sam_write1(disc, hdr, aln);
		//	continue;
        //}


    }

	// Collect positions and counts from both hash maps
    PositionClipped *positions = NULL;
    int num_positions = 0;
	fprintf(stderr, "[vorTEx: get_clip_regions] extracted clips and discordants from %d total aligned reads\n", aligned_reads);
	fprintf(stderr, "[vorTEx: get_clip_regions] collecting positions...");
    collect_positions(start_clipped_map, end_clipped_map, &positions, &num_positions);

    // Sort the positions
	fprintf(stderr, "[vorTEx: get_clip_regions] sorting positions...");
    qsort(positions, num_positions, sizeof(PositionClipped), compare_positions);

    // Merge positions that are within 50 bp of each other
	fprintf(stderr, "[vorTEx: get_clip_regions] merging within 50bp positions...");
    merge_positions(positions, &num_positions);

    // Print the sorted list
    print_sorted_positions(positions, num_positions, output_file_name);


    bam_destroy1(aln);
    bam_hdr_destroy(hdr);
    sam_close(in);

    // Clean up hash maps and positions
    kh_destroy(uint64_t, start_clipped_map);
    kh_destroy(uint64_t, end_clipped_map);
    free(positions);

    fprintf(stderr, "[vorTEx: get_reads] extracted clips and discordants from %d total aligned reads\n", aligned_reads);
    if(ret < -1) {
        errx(1, "get_reads: error reading bam: %s\n", bam_file_name);
    }
    return 0;
}
