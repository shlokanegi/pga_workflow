#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import re
import pysam

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Chunk HG002 reference based on Shasta alignments.")
    parser.add_argument("csv_file", help="Input CSV file with alignment data.")
    parser.add_argument("ref_fasta", help="Reference FASTA file to extract sequences from.")
    parser.add_argument("output_dir", help="Output directory to write the TSV and FASTA files.")
    parser.add_argument("output_prefix", help="Prefix for the output TSV and FASTA files.")
    return parser.parse_args()

def calculate_coordinates(df: pd.DataFrame) -> dict[str, dict[str, int]]:
    """
    Calculate the absolute start and end coordinates for each haplotype.
    Egs output: 
    {
        "chr6_PATERNAL": {"start": 100000, "end": 200000},
        "chr6_MATERNAL": {"start": 200000, "end": 300000}
    }
    """
    coords = {}
    
    # Group by the full 'Reference name' to handle each unique reference region
    for ref_name, group in df.groupby('Reference name'):
        # Extract haplotype and offset from the reference name
        match = re.match(r'([^:]+):(\d+)-(\d+)', ref_name)
        if not match:
            print(f"Warning: Could not parse reference name '{ref_name}'. Skipping.")
            continue
            
        haplotype, offset_start, offset_end = match.groups()
        offset_start = int(offset_start)
        offset_end = int(offset_end)
        
        # Find the min start and max end of alignments within this reference chunk
        min_ref_begin = group['Reference begin'].min()
        max_ref_end = group['Reference end'].max()
        
        # Calculate absolute coordinates
        abs_start = offset_start + min_ref_begin
        abs_end = offset_start + max_ref_end
        
        # Store the results
        if haplotype in coords:
             # If a haplotype appears with different regions, take the widest span
            coords[haplotype]['start'] = min(coords[haplotype]['start'], abs_start)
            coords[haplotype]['end'] = max(coords[haplotype]['end'], abs_end)
        else:
            coords[haplotype] = {'start': abs_start, 'end': abs_end}
            
    return coords

def write_coords_tsv(coords: dict[str, dict[str, int]], output_path: str):
    """Write the calculated coordinates to a TSV file."""
    with open(output_path, 'w') as f:
        f.write("haplotype\tstart\tend\n")
        for haplotype, data in coords.items():
            f.write(f"{haplotype}\t{data['start']}\t{data['end']}\n")
    print(f"Wrote coordinates to {output_path}")

def generate_chunked_fasta(coords: dict[str, dict[str, int]], ref_fasta_path: str, output_path: str):
    """
    Generate a chunked reference FASTA file using the calculated coordinates.
    """
    try:
        with pysam.FastaFile(ref_fasta_path) as ref_fasta, open(output_path, 'w') as out_fasta:
            for haplotype, data in coords.items():
                # pysam fetch is 0-based, our coordinates are 1-based.
                # Convert to 0-based start, 0-based end (exclusive) for fetch.
                start = data['start'] - 1
                end = data['end']
                
                # Fetch the sequence
                seq = ref_fasta.fetch(haplotype, start, end)
                
                # Write to the output FASTA file
                header = f">{haplotype}:{data['start']}-{data['end']}"
                out_fasta.write(header + "\n")
                out_fasta.write(seq + "\n")
        print(f"Successfully generated chunked FASTA: {output_path}")
    except FileNotFoundError:
        print(f"Error: Reference FASTA file not found at {ref_fasta_path}")
    except KeyError as e:
        print(f"Error: Contig {e} not found in the reference FASTA. Please check haplotype names.")
    except Exception as e:
        print(f"An error occurred during FASTA generation: {e}")


def main():
    args = parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define output file paths
    output_tsv = os.path.join(args.output_dir, f"{args.output_prefix}.for_hifiasm.coords.tsv")
    output_fasta = os.path.join(args.output_dir, f"{args.output_prefix}.for_hifiasm.fasta")

    try:
        # Load the CSV file
        df = pd.read_csv(args.csv_file)
        
        # Clean up column names by stripping whitespace
        df.columns = df.columns.str.strip()
        
    except FileNotFoundError:
        print(f"Error: Input file not found at {args.csv_file}")
        return
    except Exception as e:
        print(f"Error reading or parsing CSV file: {e}")
        return

    # Calculate coordinates
    coords = calculate_coordinates(df)
    
    if not coords:
        print("Error: No valid coordinates could be calculated. Please check the 'Reference name' column format in the input CSV.")
        return

    # Write coordinates to TSV
    write_coords_tsv(coords, output_tsv)
    
    # Generate the chunked FASTA file
    generate_chunked_fasta(coords, args.ref_fasta, output_fasta)


if __name__ == "__main__":
    main()
