from Bio import SeqIO
import pandas as pd
import json
import argparse

# Define complement dictionary
complement_map = str.maketrans("ACGTacgt", "TGCAtgca")

# Function to get complement
def complement(seq: str):
    return seq.translate(complement_map)

# This function remains the same, but it will be called with a pre-loaded sequence.
def get_anchor_sequence(sequence: str, strand: int, start: int, end: int):
    """
    Extracts anchor sequence from a given read sequence based on strand information.
    """
    return (complement(sequence[::-1][start:end]) if strand == 1 else sequence[start:end])

def main():
    parser = argparse.ArgumentParser(
        description="Generates a TSV with anchor sequences",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('-j', '--json', required=True, help="Path to the anchors JSON file")
    parser.add_argument('-m', '--master_tsv', required=True, help="Path to the master table TSV file")
    parser.add_argument('-f', '--fasta', required=True, help="Path to the FASTA file")
    parser.add_argument('-o', '--output_prefix', required=True, help="Output prefix of the TSV file")
    args = parser.parse_args()

    # Step 1: Parse the JSON file to get anchors to read info mapping
    print("Step 1: Parsing the anchors JSON file...")
    with open(args.json) as f:
        anchors_json = json.load(f)
    # Use a dictionary comprehension for a more efficient and concise conversion
    anchor_readinfo_dict = {anchor: read_list[0] for anchor, read_list in anchors_json}

    # -- IMPROVEMENT 1: Read the entire FASTA file into a dictionary ONCE --
    print("Pre-loading FASTA file into memory for fast lookups...")
    fasta_records = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
    print(f"Loaded {len(fasta_records)} records from the FASTA file.")

    # Step 2: Parse the master table TSV
    print("Step 2: Parsing the master table TSV file...")
    master_table = pd.read_csv(args.master_tsv, sep="\t")
    
    # -- BUG FIX: Use a list to select multiple columns with .iloc --
    # Note: It's often safer to select by column name if the order isn't guaranteed.
    # If column names are stable, you could use:
    # master_table = master_table[['col_name_0', 'col_name_1', ...]]
    column_indices_to_keep = [0, 1, 2, 3, 4, 6, 10, 15]
    master_table = master_table.iloc[:, column_indices_to_keep]

    # -- IMPROVEMENT 2: Use .apply() instead of a slow iterrows() loop --
    print("Extracting anchor sequences...")
    
    def extract_sequence_for_row(row):
        anchor_id = row['Anchor_path']
        read_info = anchor_readinfo_dict.get(anchor_id)
        
        if read_info is None:
            print(f"Error: Anchor {anchor_id} not found in the anchors JSON file.")
            # Decide how to handle missing anchors. Here we return None.
            return None 
            
        read_id, read_strand, read_start, read_end = read_info
        
        # Fast lookup from the pre-loaded dictionary
        read_record = fasta_records.get(read_id)
        if read_record is None:
            print(f"Error: Read ID {read_id} for anchor {anchor_id} not found in FASTA file.")
            return None

        # Extract sequence using the helper function
        return get_anchor_sequence(str(read_record.seq), read_strand, read_start, read_end)

    # Create the new column by applying the function to each row
    master_table['Anchor_sequence'] = master_table.apply(extract_sequence_for_row, axis=1)
    
    # Step 3: Save the updated master table to a TSV file
    output_file = f"{args.output_prefix}_anchors_sequences.tsv"
    print(f"Step 3: Saving the updated master table to {output_file}...")
    master_table.to_csv(output_file, sep="\t", index=False)
    print("Processing completed successfully.")

if __name__ == "__main__":
    main()
