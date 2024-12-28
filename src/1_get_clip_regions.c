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
	uint64_t key = ((uint64_t)chromosome_id << 32) | (uint64_t)position;
	uint64_t chr = (uint64_t)chromosome_id << 32;
	//printf("encoded: chr:%d\nshifted chr:%llu\npos:%d\nkey:%llu\n",chromosome_id,chr, position,key);
    return key;
}

// Function to decode the unique integer key back into chromosome ID and position
void decode_position(uint64_t key, int *chromosome_id, int *position) {
    *chromosome_id = (key >> 32);
    *position = key & 0xFFFFFFFF;  // Lower 32 bits are the position
}


// Structure to store position and its corresponding start and end clipped counts
typedef struct {
	int chromosome;  // Chromosome name (not used directly in hash map)
    int pos;
    int start_count;
    int end_count;
	int discordant_count;
} PositionClipped;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Structure to represent a BED entry
typedef struct {
    char chromosome[50];
    int pos1;
    int pos2;
    char repeat_family[100];
} BedEntry;

// Function to read BED entries from a file into an array
BedEntry *read_bed_file(const char *filename, size_t *bed_size) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return NULL;
    }

    BedEntry *bed = NULL;
    *bed_size = 0;
    size_t capacity = 100;
    bed = malloc(capacity * sizeof(BedEntry));
    if (!bed) {
        fprintf(stderr, "Memory allocation failed\n");
        fclose(file);
        return NULL;
    }

    while (fscanf(file, "%49s %d %d %99s", 
                  bed[*bed_size].chromosome, 
                  &bed[*bed_size].pos1, 
                  &bed[*bed_size].pos2, 
                  bed[*bed_size].repeat_family) == 4) {
        (*bed_size)++;
        if (*bed_size == capacity) {
            capacity *= 2;
            bed = realloc(bed, capacity * sizeof(BedEntry));
            if (!bed) {
                fprintf(stderr, "Memory allocation failed\n");
                fclose(file);
                return NULL;
            }
        }
    }
    fclose(file);

    return bed;
}

// Function to compare a query chromosome with a BedEntry
int compare_chromosome(const char *chromosome, BedEntry *entry) {
    return strcmp(chromosome, entry->chromosome);
}

// Optimized binary search to find if a position is within any range
int is_mate_in_repeat(BedEntry *bed, int bed_size, const char *chromosome, int position) {
    int left = 0, right = bed_size - 1;

    // Find the first entry of the target chromosome using binary search
    while (left <= right) {
        int mid = (left + right) / 2;
        int cmp = compare_chromosome(chromosome, &bed[mid]);

        if (cmp > 0) {
            left = mid + 1;
        } else if (cmp < 0) {
            right = mid - 1;
        } else {
            right = mid - 1; // Keep searching for the first occurrence
        }
    }

    // Start searching within the narrowed range
    int start_index = left;
    if (start_index >= bed_size || compare_chromosome(chromosome, &bed[start_index]) != 0) {
        return 0; // Chromosome not found
    }

    // Iterate through the relevant entries and check position
	int i;
    for (i = start_index; i < bed_size && strcmp(bed[i].chromosome, chromosome) == 0; i++) {
        if (position >= bed[i].pos1 && position <= bed[i].pos2) {
            return 1; // Position found
        }
        if (position < bed[i].pos1) {
            break; // Stop early since the BED file is sorted
        }
    }

    return 0; // Position not found
}



//char int_to_base(int base_int) {
//    switch (base_int) {
//        case 1: return 'A';
//        case 2: return 'C';
//        case 4: return 'G';
//        case 8: return 'T';
//        case 15: return 'T';
//        default: return '?';  // 'N' for unknown or invalid base
//    }
//}



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
    } else if ((cigar[b->core.n_cigar - 1] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP) {
		// Check if there's a soft clip at the end (right side)
        *clip_side |= 2;  // Mark right-side soft clip

		int clip_len = cigar[b->core.n_cigar - 1];
		int position = b_pos;

		int i;
		for (i = 0; i < b->core.n_cigar; i++) {
    	    int operation = cigar[i] & BAM_CIGAR_MASK;  // Get the operation type
    	    int length = cigar[i] >> BAM_CIGAR_SHIFT;     // Get the length of the operation

    	    // Handle soft-clipping
    	    if (operation != BAM_CSOFT_CLIP && operation != BAM_CHARD_CLIP && 
					operation != BAM_CINS) {
    	        // If it's a match, update the last matched position
    	        position += length;  // Move the position forward
    	    } 
    	}

		uint64_t key = encode_position(chromosome_id, position);

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
void collect_positions(khash_t(uint64_t) *start_clipped_map, khash_t(uint64_t) *end_clipped_map, khash_t(uint64_t) *discordant_map,
	PositionClipped **positions, int *num_positions) {
    // Start with the positions from the start_clipped_map
	int chrom_id,start_pos,end_pos;
	khiter_t k;

	// Process positions from the start_clipped_map
    for (k = kh_begin(start_clipped_map); k != kh_end(start_clipped_map); ++k) {
        if (kh_exist(start_clipped_map, k)) {
            uint64_t key = kh_key(start_clipped_map, k);
            decode_position(key, &chrom_id, &start_pos);
            int start_count = kh_value(start_clipped_map, k);

            // Check if the position is also in the end_clipped_map
            khiter_t k_end = kh_get(uint64_t, end_clipped_map, key);
            int end_count = (k_end != kh_end(end_clipped_map)) ? kh_value(end_clipped_map, k_end) : 0;

            // Check if the position is also in the discordant_map
            khiter_t k_disc = kh_get(uint64_t, discordant_map, key);
            int discordant_count = (k_disc != kh_end(discordant_map)) ? kh_value(discordant_map, k_disc) : 0;

            // Add or update the position in the list
            (*positions) = realloc(*positions, (*num_positions + 1) * sizeof(PositionClipped));
            (*positions)[*num_positions].chromosome = chrom_id;
            (*positions)[*num_positions].pos = start_pos;
            (*positions)[*num_positions].start_count = start_count;
            (*positions)[*num_positions].end_count = end_count;
            (*positions)[*num_positions].discordant_count = discordant_count;
            (*num_positions)++;
        }
    }

	// Now process positions that exist only in the end_clipped_map
    for (k = kh_begin(end_clipped_map); k != kh_end(end_clipped_map); ++k) {
        if (kh_exist(end_clipped_map, k)) {
            uint64_t key = kh_key(end_clipped_map, k);
            decode_position(key, &chrom_id, &end_pos);
            int end_count = kh_value(end_clipped_map, k);

            // Check if the position is already in the start_clipped_map
            khiter_t k_start = kh_get(uint64_t, start_clipped_map, key);
            if (k_start == kh_end(start_clipped_map)) {
                // Position exists only in the end_clipped_map, no start_count available
                int start_count = 0;

                // Check if the position is also in the discordant_map
                khiter_t k_disc = kh_get(uint64_t, discordant_map, key);
                int discordant_count = (k_disc != kh_end(discordant_map)) ? kh_value(discordant_map, k_disc) : 0;

                // Add the position to the list
                (*positions) = realloc(*positions, (*num_positions + 1) * sizeof(PositionClipped));
                (*positions)[*num_positions].chromosome = chrom_id;
                (*positions)[*num_positions].pos = end_pos;
                (*positions)[*num_positions].start_count = start_count;
                (*positions)[*num_positions].end_count = end_count;
                (*positions)[*num_positions].discordant_count = discordant_count;
                (*num_positions)++;
            }
        }
    }

	// Now process positions from the discordant_map
    for (k = kh_begin(discordant_map); k != kh_end(discordant_map); ++k) {
        if (kh_exist(discordant_map, k)) {
            uint64_t key = kh_key(discordant_map, k);
            decode_position(key, &chrom_id, &start_pos);
            int discordant_count = kh_value(discordant_map, k);

            // Check if the position is already in the start_clipped_map
            khiter_t k_start = kh_get(uint64_t, start_clipped_map, key);
            int start_count = (k_start != kh_end(start_clipped_map)) ? kh_value(start_clipped_map, k_start) : 0;

            // Check if the position is already in the end_clipped_map
            khiter_t k_end = kh_get(uint64_t, end_clipped_map, key);
            int end_count = (k_end != kh_end(end_clipped_map)) ? kh_value(end_clipped_map, k_end) : 0;

            // Add the position to the list
            (*positions) = realloc(*positions, (*num_positions + 1) * sizeof(PositionClipped));
            (*positions)[*num_positions].chromosome = chrom_id;
            (*positions)[*num_positions].pos = start_pos;
            (*positions)[*num_positions].start_count = start_count;
            (*positions)[*num_positions].end_count = end_count;
            (*positions)[*num_positions].discordant_count = discordant_count;
            (*num_positions)++;
        }
    }

}

// Comparison function for sorting the positions by their position values
int compare_position_clipped(const void *a, const void *b) {
    // Cast the pointers to PositionClipped
    PositionClipped *pos1 = (PositionClipped *)a;
    PositionClipped *pos2 = (PositionClipped *)b;
    
    // First compare chromosome names (lexicographical order)
    if (pos1->chromosome != pos2->chromosome) {
        return pos1->chromosome - pos2->chromosome;  // Return the result of chromosome name comparison
    }

    // If chromosomes are the same, compare positions (numerical order)
    return pos1->pos  - pos2->pos;
}

void merge_positions(PositionClipped *positions, int *num_positions) {
    if (*num_positions <= 1) return;  // If there's 0 or 1 position, no merge is needed.

    int write_index = 1;  // Start from the second position
    int i;

    for (i = 1; i < *num_positions; i++) {
        // Case 1: If the current position has discordant_count > 0 and start_count == 0, end_count == 0
        // Merge it with the closest position (before or after)
        if (positions[i].start_count == 0 && positions[i].end_count == 0 && positions[i].discordant_count > 0) {
            // Find the closest position (before or after)
            int closest_index = write_index - 1;
            int closest_distance = abs(positions[i].pos - positions[write_index - 1].pos);  // Initial distance to the previous position

            // Check the next position if it exists
            if (i + 1 < *num_positions) {
                int next_distance = abs(positions[i].pos - positions[i + 1].pos);
                if (next_distance < closest_distance) {
                    closest_index = i + 1;
                    closest_distance = next_distance;
                }
            }

            // Merge the discordant count with the closest position
			if (closest_distance < 2000) {
				positions[closest_index].discordant_count += positions[i].discordant_count;
			} else {
				positions[write_index] = positions[i];
				write_index++;
			}
        }
        // Case 2: Otherwise, if the current position is within 100 bp of the previous one and on the same chromosome
        else if (positions[i].pos - positions[write_index - 1].pos <= 100 &&
                 positions[i].chromosome == positions[write_index - 1].chromosome) {
            // Merge all counts (start_count, end_count, discordant_count)
            positions[write_index - 1].start_count += positions[i].start_count;
            positions[write_index - 1].end_count += positions[i].end_count;
            positions[write_index - 1].discordant_count += positions[i].discordant_count;
        } else {
            // Case 3: If the current position is not within 100 bp of the previous position, move it to the next available spot
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
    fprintf(output_file, "chr\tpos1\tpos2\tid\tstart_clipped\tend_clipped\tdiscordant_in_rpt\n");

	int chrom_id;
	int position;

	int i;
    for (i = 0; i < num_positions; i++) {
		if (positions[i].start_count + positions[i].end_count >= 2) {
			chrom_id = positions[i].chromosome;
            fprintf(output_file, "%s\t%d\t%d\t%s:%d\t%d\t%d\t%d\n",chromosome_names[chrom_id], positions[i].pos,
					positions[i].pos, 
					chromosome_names[chrom_id], positions[i].pos,
					positions[i].start_count, positions[i].end_count, positions[i].discordant_count);
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

//void rename_qname(bam1_t *b, bam_hdr_t *header) {
//    // Get chromosome name from reference (chromosome is based on reference id)
//	const char *chrom= sam_hdr_tid2name(header, b->core.tid);
//
//	const char *mate_ref_name = "none";
//	int mate_pos = 0;
//
//	if (b->core.mtid >= 0) {
//        // Mate exists (mtid >= 0 indicates the mate's chromosome is set)
//        mate_ref_name = sam_hdr_tid2name(header, b->core.mtid);
//        mate_pos = b->core.mpos;
//    } 
//    
//    // Get position from the bam1_t record
//    int position = b->core.pos + 1; // BAM positions are 0-based, so we add 1 to get 1-based position
//    
//    // Get CIGAR string (it's an array of operations)
//    char cigar_str[1024];
//    int cigar_len = 0;
//    int i;  // Declare the loop index variable outside the loop
//    for (i = 0; i < b->core.n_cigar; i++) {
//        uint32_t op = bam_get_cigar(b)[i] & BAM_CIGAR_MASK;  // Corrected: Access b->core.cigar
//        uint32_t len = bam_get_cigar(b)[i] >> BAM_CIGAR_SHIFT;
//        // Append the CIGAR operation to the string
//        cigar_len += snprintf(cigar_str + cigar_len, sizeof(cigar_str) - cigar_len, "%d%c", len, bam_cigar_opchr(op));
//    }
//    
//    // Create a new qname string (CIGAR + chromosome + position)
//    char new_qname[1024];
//    snprintf(new_qname, sizeof(new_qname), "%s::%s::%d::%s::%d::%s", cigar_str, chrom, position,mate_ref_name,mate_pos,bam_get_qname(b));
//    
//	bam_set_qname(b,new_qname);
//}






int main(int argc, char **argv)
{
    // without these 2 lines, htslib sometimes tries to download a part of the sequence
    // even though the -f reference was provided.
    setenv("REF_CACHE", "", 0);
    setenv("REF_PATH", "fake_value_so_no_download", 0);
    if (argc < 5)
        errx(1,
                "usage\t:%s -f <optional-reference> -e <exclude_bed_file> -r <rmsk_file_name> <bam> <clip_position.tsv> (optional #threads)",
                argv[0]);


    int aligned_reads = 0; // track total aligned reads. useful for QC without 2nd pass over full cram.

	// Initialize hash maps for start and end clipped positions
    khash_t(uint64_t) *start_clipped_map = kh_init(uint64_t);
    khash_t(uint64_t) *end_clipped_map = kh_init(uint64_t);
    khash_t(uint64_t) *discordant_map = kh_init(uint64_t);

    int opt;
    char *fasta = NULL;
    char *exclude_regions_file_name;
    char *rmsk_file_name;
    while ((opt = getopt(argc, argv, "f:e:r:")) != -1) {
        switch (opt) {
            case 'f': // Input file
                fasta = optarg;
                break;
            case 'e': // Another option
                exclude_regions_file_name = optarg;
                break;
            case 'r': // Another option
                rmsk_file_name = optarg;
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

    // Read BED file into an array
    size_t bed_size;
    BedEntry *bed = read_bed_file(rmsk_file_name, &bed_size);
    if (!bed) {
        return 1; // Error reading BED file
    }
	fprintf(stderr,"[vorTEx: get_clip_regions] loaded rmsk and excluded regions\n");

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


		// DISCORDANTS READS
		//read is not mapped in proper pair (large insert size), the mate IS (not un) mapped and NOT a suppl aln
        if ((((aln->core.flag) & (BAM_FPROPER_PAIR | BAM_FMUNMAP | BAM_FSUPPLEMENTARY)) == 0) && 
				(aln->core.flag & BAM_FPAIRED == 1) && aln->core.mtid >= 0) {
			const char *mate_chromosome = hdr->target_name[aln->core.mtid];
        	int mate_position = aln->core.mpos;

			if (is_mate_in_repeat(bed, bed_size, mate_chromosome, mate_position)) {
				//fprintf(stderr, "[vorTEx: get_clip_regions] discordant with mate repeat\n", aligned_reads);
				const char *chrom= sam_hdr_tid2name(hdr, aln->core.tid);
    			int chromosome_id = get_chromosome_id(chrom);
				if (chromosome_id == -1) { continue; };
				int b_pos = aln->core.pos;

				uint64_t key = encode_position(chromosome_id, b_pos);

				// Check if start_pos already exists, otherwise add it to the hash map
        		int ret;
        		khiter_t k = kh_get(uint64_t, discordant_map, key);
        		if (k == kh_end(discordant_map)) {  // Position doesn't exist yet
        		    k = kh_put(uint64_t, discordant_map, key, &ret);
        		    kh_value(discordant_map, k) = 1;
        		} else {
        		    kh_value(discordant_map, k)++;
        		}
			}
        }


    }

	// Collect positions and counts from both hash maps
    PositionClipped *positions = NULL;
    int num_positions = 0;
	fprintf(stderr, "[vorTEx: get_clip_regions] extracted clips and discordants from %d total aligned reads\n", aligned_reads);
	fprintf(stderr, "[vorTEx: get_clip_regions] collecting positions...\n");
    collect_positions(start_clipped_map, end_clipped_map, discordant_map, &positions, &num_positions);

    // Sort the positions
	fprintf(stderr, "[vorTEx: get_clip_regions] sorting positions...\n");
    qsort(positions, num_positions, sizeof(PositionClipped), compare_position_clipped);

    // Merge positions that are within 50 bp of each other
	fprintf(stderr, "[vorTEx: get_clip_regions] merging within 50bp positions...\n");
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
	free(bed);

    fprintf(stderr, "[vorTEx: get_reads] extracted clips and discordants from %d total aligned reads\n", aligned_reads);
    if(ret < -1) {
        errx(1, "get_reads: error reading bam: %s\n", bam_file_name);
    }
    return 0;
}
