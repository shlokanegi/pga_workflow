import os
import argparse
import pandas as pd
import numpy as np
import json

def parse_arguments():
    parser = argparse.ArgumentParser(description="Translates vg_anchor output files into usable TSVs")
    # parser.add_argument('-p', '--pkl', required=True, help='Anchors pkl file with sentinal and anchor info')    
    parser.add_argument('-j', '--json', required=True, help='anchors.json file, with anchors and aligned reads info')    
    return parser.parse_args()

def get_anchors_to_reads_tsv(json_file):
	"""
	This function converts the json (output of `vg_anchors get-anchors`) into a dictionary and returns it.
	`anchors_read_dict` has anchors as keys. Value is a list with read information [read_name, read_strand, read_start, read_end].
	`read_start` mean position of start of the anchor in the read
	`read_end` mean position of end of the anchor in the read

	Parameters
	----------
	json_file:
		The filepath of the json file

	Returns
	-------
	Dictionary with anchors mapped to all supporting read information
	"""

	records = []
	with open(json_file, 'r') as jf:
		line = jf.read()
		records.append(json.loads(line))
	records = records[0]

	# Print all anchors
	anchor_read_dict = {}
	for i in range(len(records)):
		if records[i][0] not in anchor_read_dict:
			anchor_read_dict[records[i][0]] = []
		anchor_read_dict[records[i][0]]=records[i][1]

	return anchor_read_dict

def main():
	args = parse_arguments()
	# misc tasks
	input_dir = os.path.dirname(args.json)
	output_file = os.path.join(input_dir, "extended_anchor_reads_info.tsv")
	# run function
	anchor_read_dict = get_anchors_to_reads_tsv(args.json)
	with open(output_file, 'w') as f:
		f.write(f"Anchor_path\tAnchor_length\tRead_count\n")
		for anchor, reads in anchor_read_dict.items():
			anchor_length=reads[0][3]-reads[0][2]
			f.write(f"{anchor}\t{anchor_length}\t{len(reads)}\n")
	print("Output anchor to read TSV file saved in {}".format(output_file))

if __name__ == '__main__':
    main()
