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


# Create a temporary file to store depth information
temp_depth_file=$(mktemp)

# Output header to the output file
echo -e "Position\tStart Clipped Count\tEnd Clipped Count\tAvg Depth (1000bp window)" > "$output_file"

# Function to calculate average depth over a 1000bp window
calculate_avg_depth() {
    local pos=$1
    local bam_file=$2
    local temp_depth_file=$3
    
    # Get the range for the 1000bp window around the position
    local start_pos=$((pos - 500))
    local end_pos=$((pos + 500))

    # Run samtools depth on the BAM file to get the coverage for the window
    samtools depth -r "${bam_file}:${start_pos}-${end_pos}" "$bam_file" > "$temp_depth_file"
    
    # Calculate the total depth in the 1000bp window
    local num_positions=1000
	local total_depth=$((awk '{s+=$1}END{print s}' "$temp_depth_file"))

    avg_depth=$(echo "$total_depth / $num_positions" | bc -l)

    # Return the position and average depth
    echo -e "$pos\t$avg_depth"
}

export -f calculate_avg_depth

# Function to process each line in the TSV file
process_tsv_line() {
    local line="$1"
    local bam_file="$2"
    local temp_depth_file="$3"
    local output_file="$4"
    
    # Skip header or empty lines
    if [[ "$line" == "Position" || -z "$line" ]]; then
        return
    fi

    # Parse the line into position, start_count, and end_count
    pos=$(echo "$line" | cut -f1)
    start_count=$(echo "$line" | cut -f2)
    end_count=$(echo "$line" | cut -f3)
    
    # Calculate the average depth for this position
    avg_depth=$(calculate_avg_depth "$pos" "$bam_file" "$temp_depth_file")
    
    # Write the results to the output file (use `echo` to append to the output file)
    echo -e "$pos\t$start_count\t$end_count\t$avg_depth" >> "$output_file"
}

export -f process_tsv_line

# Export variables
export bam_file
export temp_depth_file
export output_file

# Read the TSV file and process each line in parallel
cat "$input_tsv" | parallel --no-notice -j 8 "process_tsv_line {} $bam_file $temp_depth_file $output_file"

# Clean up temporary file
rm "$temp_depth_file"

echo "Processing complete. Results written to $output_file"

