#!/usr/bin/env python

import argparse
import pandas as pd
import re
from collections import defaultdict

def main():
    """
    This script takes a read processed TSV file as input and calculates the
    read depth, median MAPQ, and aggregates read IDs/MAPQs for each node.
    """
    parser = argparse.ArgumentParser(
        description="Calculate node read depth and MAPQ stats from a processed TSV file."
    )
    parser.add_argument(
        "processed_tsv",
        type=str,
        help="Path to the input read processed TSV file.",
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="nodes_info.tsv",
        help="Path to the output TSV file. (default: nodes_info.tsv)"
    )
    args = parser.parse_args()

    # --- Step 1: Read the TSV and aggregate data for each node ---
    print("Reading processed TSV and aggregating data...")
    col_names = [
        "read_name", "read_len", "relative_strand", "mapq", "div",
        "path_start", "path_end", "node_list", "orientation_list", "cs_line"
    ]
    try:
        df = pd.read_csv(args.processed_tsv, sep='\t', header=None, names=col_names)
    except FileNotFoundError:
        print(f"Error: Processed TSV not found at {args.processed_tsv}")
        return

    df.dropna(subset=['node_list', 'read_name', 'mapq', 'div'], inplace=True)

    # Use defaultdict for cleaner initialization of the data structure
    node_data = defaultdict(lambda: {'reads': [], 'mapqs': [], 'divs': []})
    
    for _, row in df.iterrows():
        read_name = row['read_name']
        mapq = row['mapq']
        div = row['div']
        
        # Use a set to only count each read once per node
        nodes = set(re.findall(r'\d+', str(row['node_list'])))
        
        if not nodes:
            continue
            
        for node in nodes:
            node_data[node]['reads'].append(read_name)
            node_data[node]['mapqs'].append(mapq)
            node_data[node]['divs'].append(div)

    if not node_data:
        print("No nodes found in the input file. No output generated.")
        return

    # --- Step 2: Process aggregated data into final DataFrame ---
    print(f"Found {len(node_data)} unique nodes. Processing for output...")
    
    output_rows = []
    for node_id, data in node_data.items():
        mapq_list = data['mapqs']
        div_list = data['divs']
        output_rows.append({
            'nodeID': int(node_id),
            'read_depth': len(data['reads']),
            'read_IDs': ",".join(data['reads']),
            'mapqs': ",".join(map(str, mapq_list)),
            'median_mapq': pd.Series(mapq_list).median(),
            'divs': ",".join(map(str, div_list)),
            'median_div': pd.Series(div_list).median()
        })
    
    final_df = pd.DataFrame(output_rows)
    final_df.sort_values('nodeID', inplace=True)

    final_df.to_csv(args.output, sep='\t', index=False)
    print(f"Successfully created node stats TSV file: {args.output}")

if __name__ == "__main__":
    main()
