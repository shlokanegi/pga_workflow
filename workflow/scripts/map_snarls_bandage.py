#!/usr/bin/env python

import argparse
import pandas as pd
import re

def main():
    """
    This script identifies reliable and unreliable snarls and
    generates a Bandage CSV to visualize all nodes belonging to them, including
    their read depth. Unreliable snarls are colored blue, and reliable snarls
    are colored red.
    """
    parser = argparse.ArgumentParser(
        description="Generate a Bandage CSV for visualizing reliable and unreliable snarls with read depth."
    )
    parser.add_argument(
        "snarl_dict_csv",
        type=str,
        help="Path to the snarl dictionary CSV file (snarl_start, snarl_end, snarls_inside).",
    )
    parser.add_argument(
        "snarl_compatibility_tsv",
        type=str,
        help="Path to the snarl compatibility fractions TSV file.",
    )
    parser.add_argument(
        "nodes_info_tsv",
        type=str,
        help="Path to the nodes info TSV file (nodeID, read_depth, ...).",
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="snarls.bandage.csv",
        help="Path to the output Bandage CSV file. (default: snarls.bandage.csv)"
    )
    args = parser.parse_args()

    # --- Step 1: Find reliable and unreliable snarl IDs ---
    try:
        compatibility_df = pd.read_csv(args.snarl_compatibility_tsv, sep='\t')
    except FileNotFoundError:
        print(f"Error: Compatibility TSV not found at {args.snarl_compatibility_tsv}")
        return

    unreliable_snarls = compatibility_df[compatibility_df['fraction_compatible'] == 0]
    unreliable_snarl_ids = set(unreliable_snarls['source_snarl_id'])

    reliable_snarls = compatibility_df[compatibility_df['fraction_compatible'] > 0]
    reliable_snarl_ids = set(reliable_snarls['source_snarl_id'])

    if not unreliable_snarl_ids and not reliable_snarl_ids:
        print("No snarls found in the compatibility file. No output file generated.")
        return
    
    print(f"Found {len(unreliable_snarl_ids)} unreliable snarls and {len(reliable_snarl_ids)} reliable snarls.")

    # --- Step 2: Extract nodes from snarls ---
    unreliable_nodes = set()
    reliable_nodes = set()
    print("Reading snarl dictionary and extracting nodes...")
    try:
        with open(args.snarl_dict_csv, 'r') as f:
            for i, line in enumerate(f, 1):
                snarl_id = i
                nodes_found = [int(node) for node in re.findall(r'\d+', line)]
                if snarl_id in unreliable_snarl_ids:
                    unreliable_nodes.update(nodes_found)
                elif snarl_id in reliable_snarl_ids:
                    reliable_nodes.update(nodes_found)
    except FileNotFoundError:
        print(f"Error: Snarl dictionary CSV not found at {args.snarl_dict_csv}")
        return

    all_nodes = unreliable_nodes.union(reliable_nodes)
    print(f"Collected {len(all_nodes)} unique nodes from {len(unreliable_snarl_ids) + len(reliable_snarl_ids)} snarls.")

    if not all_nodes:
        print("No nodes found for the snarls. No output file generated.")
        return
        
    # --- Step 3: Read read depths ---
    try:
        depth_df = pd.read_csv(args.nodes_info_tsv, sep='\t', index_col='nodeID')
    except FileNotFoundError:
        print(f"Error: Node read depths TSV not found at {args.nodes_info_tsv}")
        return
        
    sorted_nodes = sorted(list(all_nodes))
    read_depths = [depth_df.loc[node, 'read_depth'] if node in depth_df.index else 0 for node in sorted_nodes]

    # --- Step 4: Generate and save the Bandage CSV ---
    colors = []
    for node in sorted_nodes:
        if node in unreliable_nodes:
            colors.append('#0000FF')  # Blue for unreliable
        else:
            colors.append('#FF0000')  # Red for reliable

    bandage_df = pd.DataFrame({
        'nodeID': sorted_nodes,
        'read_depth': read_depths,
        'Color': colors
    })
    
    bandage_df.to_csv(args.output, index=False)
    print(f"Successfully created Bandage CSV for snarls: {args.output}")

if __name__ == "__main__":
    main()
