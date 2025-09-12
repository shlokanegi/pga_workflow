import sys
import pysam
import pandas as pd

def extract_full_contigs():
    """
    Identifies contigs that align to specific genomic regions and extracts their
    full sequences from a FASTA file.
    """
    if len(sys.argv) != 5:
        print(f"Usage: python {sys.argv[0]} <bam_file> <coords_file> <fasta_file> <output_fasta>")
        sys.exit(1)
    bam_file, coords_file, fasta_file, output_fasta = sys.argv[1:5]

    try:
        coords_df = pd.read_csv(coords_file, sep="\t")
    except FileNotFoundError:
        print(f"Error: Coordinates file not found at {coords_file}")
        sys.exit(1)

    try:
        samfile = pysam.AlignmentFile(bam_file, "rb")
        sequences = pysam.FastaFile(fasta_file)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error opening BAM or FASTA file: {e}")
        sys.exit(1)

    # Step 1: Find all unique contig names that map to the target regions
    mapped_contig_names = set()
    for index, row in coords_df.iterrows():
        contig = row['haplotype']
        start = int(row['start']) - 1000000
        end = int(row['end']) + 1000000

        try:
            fetched_reads = samfile.fetch(contig, start, end)
            for read in fetched_reads:
                mapped_contig_names.add(read.query_name)
        except ValueError:
            print(f"Warning: Contig '{contig}' not found in BAM. Skipping region.")
            continue
    
    print(f"Found {len(mapped_contig_names)} unique contigs aligning to the target regions.")

    # Step 2: Write the full sequences of those contigs to the output file
    with open(output_fasta, "w") as fasta_out:
        for name in mapped_contig_names:
            try:
                sequence = sequences.fetch(name)
                fasta_out.write(f">{name}\n{sequence}\n")
            except KeyError:
                print(f"Warning: Contig '{name}' found in BAM but not in FASTA file. Skipping.")

    samfile.close()
    sequences.close()

if __name__ == "__main__":
    extract_full_contigs()
