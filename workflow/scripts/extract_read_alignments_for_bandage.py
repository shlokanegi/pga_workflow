#!/usr/bin/env python

import argparse
import pandas as pd
import re
import os
import shutil
import tempfile

def main():
    """
    This program uses a pre-calculated node stats file and a processed TSV
    of read alignments to generate a zipped archive of Bandage CSVs.
    Each CSV is for a single read and contains only the nodes it traverses.

    USAGE:
    # Example with filtering
    python scripts/extract_read_alignments_for_bandage.py \
        /path/to/read_processed.tsv \
        /path/to/nodes_info.tsv \
        --node-ids-file /path/to/your/node_ids.csv (1, 2, 3, 4, 5, etc...)

    # Example for all reads (no filter)
    python scripts/extract_read_alignments_for_bandage.py \
        /path/to/read_processed.tsv \
        /path/to/nodes_info.tsv
    """
    parser = argparse.ArgumentParser(
        description="Generate a zipped archive of per-read Bandage files."
    )
    parser.add_argument(
        "processed_tsv",
        type=str,
        help="Path to the input read processed TSV file (for read paths).",
    )
    parser.add_argument(
        "nodes_info_tsv",
        type=str,
        help="Path to the TSV file with pre-calculated node stats (e.g., node_read_depths.tsv).",
    )
    parser.add_argument(
        "--node-ids-file",
        type=str,
        default=None,
        help="Optional: Path to a file with a comma-separated list of node IDs to filter reads.",
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="read_traversals.zip",
        help="Path to the output ZIP file. (default: read_traversals.zip)"
    )
    args = parser.parse_args()

    # --- Step 1: Read pre-calculated node depths and read paths ---
    print("Reading node stats and read paths...")
    try:
        depth_df = pd.read_csv(args.nodes_info_tsv, sep='\t', usecols=['nodeID', 'read_depth'])
        depth_map = depth_df.set_index('nodeID')['read_depth']
    except FileNotFoundError:
        print(f"Error: Node info TSV not found at {args.nodes_info_tsv}")
        return
    except (ValueError, KeyError) as e:
        print(f"Error reading {args.nodes_info_tsv}. Ensure it has 'nodeID' and 'read_depth' columns. Details: {e}")
        return

    try:
        col_names = [
            "read_name", "read_len", "relative_strand", "mapq", "div",
            "path_start", "path_end", "node_list", "orientation_list", "cs_line"
        ]
        read_paths_df = pd.read_csv(args.processed_tsv, sep='\t', header=None, names=col_names)
    except FileNotFoundError:
        print(f"Error: Processed TSV not found at {args.processed_tsv}")
        return

    read_paths_df.dropna(subset=['node_list', 'read_name'], inplace=True)

    read_paths = {}
    for _, row in read_paths_df.iterrows():
        nodes = set(re.findall(r'\d+', str(row['node_list'])))
        if nodes:
            read_paths[row['read_name']] = {int(n) for n in nodes}

    if not read_paths:
        print("No valid read paths found in the input file.")
        return

    # --- Step 2: Determine which reads to process ---
    reads_to_process = []
    if args.node_ids_file:
        print("Filtering reads based on node ID file...")
        try:
            with open(args.node_ids_file, 'r') as f:
                node_ids_str = f.read()
            target_nodes = {int(node.strip()) for node in node_ids_str.split(',')}
        except FileNotFoundError:
            print(f"Error: Node ID file not found at {args.node_ids_file}")
            return

        for read_id, visited_nodes in read_paths.items():
            if not visited_nodes.isdisjoint(target_nodes):
                reads_to_process.append(read_id)
        print(f"Found {len(reads_to_process)} reads that traverse the specified nodes.")
    else:
        reads_to_process = list(read_paths.keys())
        print(f"No node list provided. Processing all {len(reads_to_process)} reads.")

    if not reads_to_process:
        print("No reads to process after filtering.")
        return

    # --- Step 3: Generate CSVs in a temporary directory and then zip them ---
    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"Generating {len(reads_to_process)} CSVs in a temporary directory...")
        for i, read_id in enumerate(reads_to_process):
            visited_nodes = sorted(list(read_paths.get(read_id, set())))
            
            if not visited_nodes:
                continue

            # Look up depths for the traversed nodes
            read_depths = depth_map.reindex(visited_nodes).fillna(0).astype(int)

            output_df = pd.DataFrame({
                'nodeID': visited_nodes,
                'read_depth': read_depths,
                'Color': 'red'
            })
            
            output_filename = os.path.join(tmpdir, f"{read_id}.bandage.csv")
            output_df.to_csv(output_filename, index=False)
            
            if (i + 1) % 500 == 0 or (i + 1) == len(reads_to_process):
                print(f"  ...processed {i + 1}/{len(reads_to_process)} reads")

        print("Zipping output files...")
        output_base_name = args.output.removesuffix('.zip')
        shutil.make_archive(output_base_name, 'zip', tmpdir)

    print(f"Successfully created zipped archive: {args.output}")

if __name__ == "__main__":
    main()
