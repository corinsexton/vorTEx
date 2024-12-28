#!/bin/bash

# Check for correct arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_tsv> <bam_file> <output_file> <number_of_threads>"
    exit 1
fi

input_tsv="$1"
bam_file="$2"
output_file="$3"

# Default to 8 threads if not provided
number_of_threads="${4:-8}"


# Validate that the number of threads is a positive integer
if ! [[ "$number_of_threads" =~ ^[0-9]+$ ]] || [ "$number_of_threads" -le 0 ]; then
    echo "Error: Number of threads must be a positive integer."
    exit 1
fi

# Create a temporary directory in /tmp (or use TMPDIR if defined)
tmp_dir="${TMPDIR:-/tmp}/process_tsv_$(date +%s)"
mkdir -p "$tmp_dir"

# Create a temporary file to store depth information within the temporary directory
temp_depth_file=$(mktemp "$tmp_dir/depth_file.XXXXXX")


# Output header to the output file
echo -e "ids\tchr\tpos1\tpos2\tdepth" > "$output_file"

# Function to calculate average depth over a 1000bp window
calculate_avg_depth() {
    local chr=$1
	local start_pos=$2
    local end_pos=$3
    local bam_file=$4
    local temp_depth_file=$5

	if [ "$start_pos" -eq "$end_pos" ]; then
    	start_pos=$((start_pos- 500))
    	end_pos=$((end_pos+ 500))
	fi
    
    # Run samtools depth on the BAM file to get the coverage for the window
    samtools depth -r "${chr}:${start_pos}-${end_pos}" "$bam_file" > "$temp_depth_file"
    
    # Calculate the total depth in the 1000bp window
	local num_positions=$((end_pos - start_pos))
	local total_depth=$(awk '{s+=$3}END{print s}' "$temp_depth_file")

	avg_depth=$(echo "scale=2; $total_depth / $num_positions" | bc -l)

    # Return the position and average depth
	if [ -z "$total_depth" ] || [ "$total_depth" -eq 0 ]; then
	    echo "Error: total depth is zero or empty"
	    return 1  # or handle this case as needed
	fi

    echo -e "$avg_depth"
}

export -f calculate_avg_depth

# Function to process each line in the TSV file
process_tsv_line() {
    local chr="$1"
    local start_pos="$2"
    local end_pos="$3"
    local id="$4"
    local bam_file="$5"
    local temp_depth_file="$6"
    local output_file="$7"
    
    # Skip header or empty lines
    if [[ "$chr" == "chr" || -z "$chr" ]]; then
        return
    fi


    # Create a unique temporary depth file for each thread
    temp_depth_file=$(mktemp "$TMPDIR/depth_file_XXXXXX")

    # Calculate the average depth for this position
    avg_depth=$(calculate_avg_depth "$chr" "$start_pos" "$end_pos" "$bam_file" "$temp_depth_file")
    
    # Write the results to the output file (use `echo` to append to the output file)
    echo -e "$id\t$chr\t$start_pos\t$end_pos\t$avg_depth" >> "$output_file"
}

export -f process_tsv_line

# Export variables
export bam_file
export temp_depth_file
export output_file

# Read the TSV file and process each line in parallel
cat "$input_tsv" | xargs -I {} -P ${number_of_threads} sh -c "process_tsv_line {} $bam_file $temp_depth_file $output_file"

#cat "$input_tsv" | parallel --no-notice -j 8 "process_tsv_line {} $bam_file $temp_depth_file $output_file"
# Clean up the temporary directory and files
rm -rf "$tmp_dir"


echo "Processing complete. Results written to $output_file"

