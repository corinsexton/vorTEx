#!/usr/bin/env python
import csv
import argparse

# Read local_depth.tsv and create a dictionary based on id
def read_local_depth(filename):
	id_map = {}  # Mapping from IDs to their local depth regions
	with open(filename, 'r') as f:
		reader = csv.reader(f, delimiter='\t')
		header = next(reader)  # Skip the header row

		for row in reader:
			depth_id, chrom, pos1, pos2, local_depth = row[0], row[1], row[2], row[3], row[4]

			# Convert positions to integers and depth to float
			local_depth = float(local_depth)

			# Store the local depth information for each ID in the id_map
			if depth_id not in id_map:
				id_map[depth_id] = []
			id_map[depth_id] = local_depth
	
	return id_map

def read_local_depth_key(filename,local_depth_dict):
	id_map = {}  # Mapping from IDs to their local depth regions
	with open(filename, 'r') as f:
		for line in f:
			depth_id, clip_id = line.strip().split('\t')

			clip_id_list = clip_id.split(",")

			for id_ in clip_id_list:
				# Store the local depth information for each ID in the id_map
				if id_ not in id_map:
					id_map[id_] = local_depth_dict[depth_id]
				else:
					print(f"Already had depth for this locus: {clip_id}")
					current_depth =  local_depth_dict[depth_id]
					id_map[id_] = min(id_map[id_],local_depth_dict[depth_id],)
	
	return id_map


# Read clip_positions.tsv and assign local depth based on ID
def assign_local_depth_to_clip_positions(depth_dict, clip_positions_file, output_file):
	with open(clip_positions_file, 'r') as clip_file:
		clip_reader = csv.reader(clip_file, delimiter='\t')
		header = next(clip_reader)  # Skip the header row
		
		# Open output file to write results
		with open(output_file, 'w', newline='') as output_f:
			writer = csv.writer(output_f, delimiter='\t')
			writer.writerow(header + ['local_depth'])  # Add a new column for local depth
			
			for clip_row in clip_reader:
				chr_clip, pos1_clip, pos2_clip, clip_id, start_clipped, end_clipped, discordant = clip_row
				start_clipped, end_clipped, discordant = int(start_clipped), int(end_clipped), int(discordant)

				local_depth = depth_dict[clip_id]


				# benchmarking:
				# 9 are not detected ever

				# current xTea cutoff:
				# 5 read support (across any depth) 21 / 1642 truth sets -> 1,009,781 candidates

				# at 20% we lose 129 / 1642 truth sets -> 288,697 candidates
				# at 15% we lose 58 / 1642 truth sets -> 457,381 candidates
				# at 10% we lose 33 / 1642 truth sets -> 743,659 candidates
				# at 8% we lose 25 / 1642 truth sets -> 938,344 candidates
				# at 6.5% we lose 18 / 1642 truth sets -> 1,172,211 candidates
				# at 5% we lose 14 / 1642 truth sets -> 1,539,039 candidates

				# TODO curious about different cutoffs for alus vs l1 vs sva
				#if (start_clipped + end_clipped) / local_depth >= 0.20 or ((start_clipped + end_clipped) / local_depth >= 0.10 and discordant > 1):
				#if (start_clipped + end_clipped) / local_depth >= 0.20:
				if ((start_clipped + end_clipped) / local_depth >= 0.10 and discordant > 1) or (start_clipped + end_clipped) > 0.25:
					writer.writerow(clip_row + [local_depth])

	print(f"Results written to {output_file}")

# Main function to orchestrate the reading and processing
def main():
	# Set up argument parser with flags
	parser = argparse.ArgumentParser(description="Assign local depth from local_depth.tsv to clip_positions.tsv based on genomic intervals.")
	
	# Define command-line arguments with short flags
	parser.add_argument("-l", "--local_depth_file", required=True, help="Path to the local depth TSV file")
	parser.add_argument("-k", "--local_depth_key_file", required=True, help="Path to the local depth key TSV file")
	parser.add_argument("-c", "--clip_positions_file", required=True, help="Path to the clip positions TSV file")
	parser.add_argument("-o", "--output_file", required=True, help="Path to the output TSV file")
	
	# Parse the command-line arguments
	args = parser.parse_args()
	
	# Read the local depth file into a dictionary based on ID
	id_map = read_local_depth(args.local_depth_file)
	depth_dict = read_local_depth_key(args.local_depth_key_file,id_map)
	print("processed depth dict...")
	
	# Process the clip positions and assign local depth
	assign_local_depth_to_clip_positions(depth_dict, args.clip_positions_file, args.output_file)

if __name__ == '__main__':
	main()

