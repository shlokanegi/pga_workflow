import sys
import json
import argparse

def process_gam_stream(gam_json_stream):
    """
    Reads a stream of JSON alignments from vg view -j, extracts refpos info,
    and returns a dictionary keyed by read name.
    """
    read_data = {}
    for line in gam_json_stream:
        try:
            alignment = json.loads(line)
            read_name = alignment.get('name')
            if not read_name:
                continue

            refpos_list = alignment.get('refpos', [])
            
            # --- Initialize ---
            chm13_data = []
            grch38_data = []

            for pos in refpos_list:
                # Assumes reference name format is GENOME#...#CHROMOSOME
                name_parts = pos.get('name', '').split('#')
                if len(name_parts) < 3:
                    continue

                genome = name_parts[0].lower()
                chrom = name_parts[2]
                offset = pos.get('offset', '0')
                strand = '-' if pos.get('is_reverse') else '+'
                
                # --- Sort data into the correct genome bucket ---
                if 'chm13' in genome:
                    chm13_data.append((chrom, offset, strand))
                elif 'grch38' in genome:
                    grch38_data.append((chrom, offset, strand))
            
            read_data[read_name] = {'chm13': chm13_data, 'grch38': grch38_data}

        except (json.JSONDecodeError, KeyError):
            # Skip problematic lines
            sys.stderr.write(f"Warning: Skipping malformed line or record: {line.strip()}\n")
            continue
            
    return read_data

def format_column_data(data_list):
    """
    Takes a list of (chrom, offset, strand) tuples and formats them
    into the separate, list-like column strings.
    """
    if not data_list:
        return "[]", "[]", "[]"
    
    # Unzip the list of tuples into separate lists for each data type
    chroms, offsets, strands = zip(*data_list)
    
    # Format each list into a string like "[item1,item2]"
    chrom_str = f"[{','.join(chroms)}]"
    offset_str = f"[{','.join(map(str, offsets))}]"
    strand_str = f"[{','.join(strands)}]"
    
    return chrom_str, offset_str, strand_str

def add_columns_to_gaf(gaf_file_path, output_file_path, read_data):
    """
    Reads a GAF file, adds the new refpos columns based on the processed
    GAM data, and writes the result to an output file.
    """
    with open(gaf_file_path, 'r') as gaf_in, open(output_file_path, 'w') as gaf_out:
        for line in gaf_in:
            line = line.strip()
            if not line or line.startswith('#'):
                gaf_out.write(line + '\n')
                continue

            columns = line.split('\t')
            read_name = columns[0]
            
            # Retrieve the processed refpos data for this read
            # Default to empty lists if the read name isn't found
            data = read_data.get(read_name, {'chm13': [], 'grch38': []})
            
            # Format the data for each genome into the required column strings
            grch38_chrom_str, grch38_offset_str, grch38_strand_str = format_column_data(data['grch38'])
            chm13_chrom_str, chm13_offset_str, chm13_strand_str = format_column_data(data['chm13'])
            
            # Append the six new columns to the original GAF columns
            new_columns = columns + [
                grch38_chrom_str,
                grch38_offset_str,
                grch38_strand_str,
                chm13_chrom_str,
                chm13_offset_str,
                chm13_strand_str
            ]
            
            gaf_out.write('\t'.join(new_columns) + '\n')

def main():
    """Main function to parse arguments and orchestrate the workflow."""
    parser = argparse.ArgumentParser(
        description="Adds reference position information from a GAM JSON stream to a GAF file. "
                    "This script reads JSON-formatted GAM alignments from standard input.",
        epilog="Example usage: vg view -j alignments.gam | python this_script.py --gaf alignments.gaf --output new.gaf"
    )
    parser.add_argument(
        '--gaf',
        required=True,
        help='Path to the input GAF file.'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Path to the output GAF file with added columns.'
    )
    args = parser.parse_args()

    # 1. Process the GAM JSON stream from stdin
    sys.stderr.write("Reading GAM JSON from stdin...\n")
    read_refpos_data = process_gam_stream(sys.stdin)
    sys.stderr.write(f"Processed data for {len(read_refpos_data)} reads from GAM stream.\n")

    # 2. Read the GAF, add columns, and write to the output file
    sys.stderr.write(f"Reading GAF file '{args.gaf}' and writing to '{args.output}'...\n")
    add_columns_to_gaf(args.gaf, args.output, read_refpos_data)
    sys.stderr.write("Done.\n")

if __name__ == '__main__':
    main()
