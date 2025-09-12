#!/usr/bin/env python3

import argparse
import configparser
from datetime import datetime

def parse_vg_anchor_log_file(filepath):
    anchor_params = {}
    try:
        with open(filepath, 'r') as f:
            line_number=1
            for line in f:
                stripped_line = line.strip()
                if line_number >= 5 and stripped_line and not stripped_line.startswith('#'):
                    key, value = line.strip().split('=', 1)
                    anchor_params[key.strip()] = value.strip()
                line_number+=1
    except FileNotFoundError:
        print(f"Warning: Could not find file {filepath}. Some parameters will be missing.")
    return anchor_params

def main():
    parser = argparse.ArgumentParser(description="Generate a summary log from various workflow configuration files.")
    parser.add_argument('--params-log', required=True, help="Path to the params_run.log file from anchor generation.")
    parser.add_argument('--shasta-conf', required=True, help="Path to the shasta anchors config file.")
    parser.add_argument('--output-log', required=True, help="Path for the final output log file (pga_run.log).")
    parser.add_argument('--run-mode', required=True, help="Value of RUN_MODE from the main config.")
    parser.add_argument('--region-id', required=True, help="Value of region_id from the main config.")
    parser.add_argument('--asm-preset', required=True, help="Value of MINIMAP:asmPreset from the main config.")
    args = parser.parse_args()

    params_from_run = parse_vg_anchor_log_file(args.params_log)
    shasta_config = configparser.ConfigParser()
    shasta_config.read(args.shasta_conf)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_content = f"""
PARAMETER LOG FOR WORKFLOW RUN
Timestamp: {timestamp}
==================================================

--- Parameters for Anchor Generation ({args.params_log}) ---
MIN_ANCHOR_LENGTH = {params_from_run.get('MIN_ANCHOR_LENGTH', 'N/A')}
EXPECTED_MAP_Q = {params_from_run.get('EXPECTED_MAP_Q', 'N/A')}
MIN_ANCHOR_READS = {params_from_run.get('MIN_ANCHOR_READS', 'N/A')}
HET_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING = {params_from_run.get('HET_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING', 'N/A')}
HOMO_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING = {params_from_run.get('HOMO_FRACTION_READS_RETAINED_THRESHOLD_FOR_MERGING', 'N/A')}
MIN_READS_REQUIRED_FOR_MERGING = {params_from_run.get('MIN_READS_REQUIRED_FOR_MERGING', 'N/A')}
FRACTION_READS_FOR_SNARL_BOUNDARY_EXTENTION = {params_from_run.get('FRACTION_READS_FOR_SNARL_BOUNDARY_EXTENTION', 'N/A')}
MIN_READS_REQUIRED_FOR_BOUNDARY_EXTENSION = {params_from_run.get('MIN_READS_REQUIRED_FOR_BOUNDARY_EXTENSION', 'N/A')}
DROP_FRACTION = {params_from_run.get('DROP_FRACTION', 'N/A')}
MIN_ANCHOR_READCOV = {params_from_run.get('MIN_ANCHOR_READCOV', 'N/A')}

# PHASING CONSISTENCY CHECK ANCHORS/SNARLS CONSTANTS
MIN_SNARL_LINKAGE_THRESHOLD = {params_from_run.get('MIN_SNARL_LINKAGE_THRESHOLD', 'N/A')}
RELIABLE_SNARL_FRACTION_THRESHOLD = {params_from_run.get('RELIABLE_SNARL_FRACTION_THRESHOLD', 'N/A')}
ADD_BACK_HOMO_SNARLS = {params_from_run.get('ADD_BACK_HOMO_SNARLS', 'N/A')}
ERROR_TOLERANCE_IN_COMPATIBILITY_CHECK = {params_from_run.get('ERROR_TOLERANCE_IN_COMPATIBILITY_CHECK', 'N/A')}
ENABLE_UNEQUAL_SET_COMPATIBILITY = {params_from_run.get('ENABLE_UNEQUAL_SET_COMPATIBILITY', 'N/A')}

--- Parameters for Shasta ({args.shasta_conf}) ---
k = {shasta_config.get('Kmers', 'k', fallback='N/A')}
Probability = {shasta_config.get('Kmers', 'probability', fallback='N/A')}
mode3.minAnchorCoverage = {shasta_config.get('Assembly', 'mode3.minAnchorCoverage', fallback='N/A')}
mode3.primaryGraph.maxLoss = {shasta_config.get('Assembly', 'mode3.primaryGraph.maxLoss', fallback='N/A')}
mode3.assemblyGraph.suppressBubbleCleanup = {shasta_config.get('Assembly', 'mode3.assemblyGraph.suppressBubbleCleanup', fallback='N/A')}
mode3.assemblyGraph.bubbleErrorThreshold = {shasta_config.get('Assembly', 'mode3.assemblyGraph.bubbleErrorThreshold', fallback='N/A')}

--- Parameters for Snakemake main workflow ---
RUN_MODE = {args.run_mode}
region_id = {args.region_id}
MINIMAP:asmPreset = {args.asm_preset}

==================================================
"""

    with open(args.output_log, "w") as f:
        f.write(log_content.strip())

if __name__ == '__main__':
    main()
