import os
import subprocess # Allows external programmes (Clustal Omega and MAFFT) to run
from Bio import AlignIO # Imports AlignIO module from Biopython package (used
# to read, write and manipulate sequence alignments

# Define folder containing fasta files and output folder produced (containing
# MSA results)
fasta_folder = "protein_sequences"
output_folder = "msa_protein_sequences"

# Ensure output folder exists
os.makedirs(output_folder, exist_ok=True) # Creates folder if absent (flags
# error if folder already exists)

# Define path for combined fasta file inside output folder
combined_fasta = os.path.join(output_folder, "combined_sequences.fasta")
# Defines name of file containing FASTA sequences for all species and saves
# it within output folder (msa_protein_sequences)

# Combine all FASTA files into single combined_sequences file
with open(combined_fasta, 'w') as outfile: # Opens combined_sequences file in
    # 'write' mode
    for fasta_file in os.listdir(fasta_folder): # Loops through all files in
        # protein_sequences FASTA folder
        if fasta_file.endswith('.fasta'):
            with open(os.path.join(fasta_folder, fasta_file)) as infile:
                outfile.write(infile.read()) # Reads contents of FASTA files
                # and writes them into combined_sequences file

print(f"Combined all FASTA files into {combined_fasta}.") # Confirms FASTA
# files have been combined into combined_sequences file

# Align combined_sequences FASTA file using Clustal Omega alignment tool
clustal_output = os.path.join(output_folder, "combined_clustal.aln")
subprocess.run([ # Calls Clustal Omega on input file
    "clustalo",
    "-i", combined_fasta, # Input file (combined_sequences)
    "-o", clustal_output, # Output file (combined_clustal)
    "--outfmt", "clu",  # Specifies Clustal format for output
    "--force" # Allows overwriting existing files if necessary
])

# Align combined_sequences FASTA file using MAFFT alignment tool
mafft_output = os.path.join(output_folder, "combined_mafft.aln")
with open(mafft_output, 'w') as mafft_out:
    subprocess.run([ # Calls MAFFT on input file
        "mafft", combined_fasta # Input file (combined_sequences)
    ], stdout=mafft_out) # Output file (combined_mafft)


def summarize_alignment(file_path):
    """
    Summarizes Clustal Omega and MAFFT alignments using Bio.AlignIO.

    Parameters:
        file_path (str): Path to Clustal Omega / MAFFT file.

    Returns:
        Edited and summarised Clustal Omega / MAFFT file.
    """
    file_format = "clustal" if "clustal" in file_path else "fasta" #
    # Determines format of alignment file
    try:
        alignment = AlignIO.read(file_path, file_format) # Reads alignment
        # file and prints summary including number of sequences, alignment
        # length, and the first 50 characters of each sequence
        print(f"\nSummary of alignment from {file_path}:")
        print(f"  Number of sequences: {len(alignment)}")
        print(f"  Alignment length: {alignment.get_alignment_length()}")
        # Calculates alignment score
        alignment_score = calculate_alignment_score(alignment)
        print(
            f"  Alignment score (percent identical positions): {alignment_score:.2f}%")
        print(f"  Sequences:")
        for record in alignment:
            print(
                f"    {record.id}: {record.seq[:50]}"
                f"{'...' if len(record.seq) > 50 else ''}")
    except Exception as e: # Prevents programme crashing and terminating due to
        # error - Alignment file cannot be parsed
        print(f"Failed to parse alignment file {file_path}. Error: {e}")


def calculate_alignment_score(alignment):
    """
    Calculates a basic alignment score as the percentage of identical positions
    across all sequences in the alignment.

    Parameters:
        alignment: Alignment variable created using AlignIO which contains
        aligned sequences.

    Returns:
        Alignment score.
    """
    alignment_length = alignment.get_alignment_length() # Retrieves length of
    # alignment (total number of positions in each sequence of the alignment)
    num_sequences = len(alignment) # Gives number of sequences present in the
    # alignment (number of records in alignment variable)
    identical_positions = 0 # Initialises number of identical positions score
    # to 0

    # Iterates over each position in the alignment
    for i in range(alignment_length):
        residues_at_position = [record.seq[i] for record in alignment]
        # Creates list of all residues at a particular position
        if (len(set(residues_at_position)) == 1 and '-' not in
                residues_at_position): # Converts residues_at_position list
            # to set (removes duplicate residues) and ensures no gaps are
            # present
            identical_positions += 1 # Increments counter if both conditions
            # are true (position is identical)

    # Calculates percentage of identical positions
    score = (identical_positions / alignment_length) * 100
    return score

# Reporting results for the combined sequences
print("\nReporting results for combined sequences:")
summarize_alignment(clustal_output)
summarize_alignment(mafft_output)