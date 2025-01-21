import os

# Define directory and file paths
script_dir = os.path.dirname(os.path.abspath(__file__))
task_2_dir = os.path.join(script_dir, "Task 2")

protein_file = os.path.join(task_2_dir, "Proteins.faa") # Path to RefSeq
# protein file
protein_frame_dir = os.path.join(task_2_dir,
                                     "2d_HFE_protein_reading_frames")
# Path to directory containing translated HFE protein frames from Task 2d
protein_orf_file = os.path.join(protein_frame_dir,
                                        "2d_HFE_ORF_sequence.fasta") # Path
# to HFE ORF sequence file from Task 2d
output_file = os.path.join(task_2_dir, "2e_HFE_protein_comparison.txt")
# Output file path and name

# Ensure the output directory exists
os.makedirs(task_2_dir, exist_ok=True)


def read_protein_sequence(file_path):
    """
    Reads ORF protein sequence from FASTA file.

    Args:
        file_path (str): Path to ORF protein sequence file.

    Returns:
        str: Protein sequence as single string.
    """
    sequence = []
    with open(file_path, "r") as file: # Open "2d_HFE_ORF_sequence.fasta" for
        # reading
        for line in file:
            if not line.startswith(">"):  # Skip header lines
                sequence.append(line.strip())  # Add sequence lines to list
    return "".join(sequence) # Combine list elements into a single string
    # and return it


def extract_refseq_sequence(file_path, ref_id):
    """
    Extracts specific sequence from a RefSeq protein FASTA file ("Protein.faa").

    Args:
        file_path (str): Path to RefSeq protein file.
        ref_id (str): Protein Reference ID of protein.

    Returns:
        str: Protein sequence corresponding to protein reference ID.
    """
    with open(file_path, "r") as file: # Open "Protein.faa" for reading
        capturing = False # Flag to determine if correct sequence is
        # being captured
        refseq_sequence = [] # Empty list to store sequence lines
        for line in file: # Iterate through each line
            if line.startswith(">"): # Check if line is a header
                if ref_id in line: # Check if protein reference ID matches
                    capturing = True # Set the flag to start capturing
                    # sequence
                    refseq_sequence = [] # Clear list for new sequence
                else:
                    capturing = False # Stop capturing if the protein
                    # reference ID does not match
            elif capturing:
                refseq_sequence.append(line.strip()) # Add sequence line
                # without whitespace to list if currently capturing the
                # correct sequence
    if refseq_sequence:
        return "".join(refseq_sequence) # Return joined list as single string
        # if sequence is successfully captured
    raise ValueError(f"Reference ID {ref_id} not found in the file.") # Raise
    # error if no matching protein reference ID is found


def compare_sequences(hfe_sequence, refseq_sequence, output_file):
    """
    Compares protein sequences ("2d_HFE_ORF_sequence.fasta" vs. protein
    sequence extracted from "Protein.faa") and writes differences to an output
    file ("2e_HFE_protein_comparison.txt").

    Args:
        hfe_sequence (str): Protein ORF sequence from Task 2d.
        refseq_sequence (str): Protein sequence from RefSeq.
        output_file (str): Path to file containing comparison results.
    """
    with open(output_file, "w") as output: # Open the output file in write mode
        if hfe_sequence == refseq_sequence: # Check if the sequences are
            # identical
            output.write("The sequences are identical.\n") # Write message
            # indicating identical sequences
        else:
            output.write("The sequences are different.\n") # Write message
            # indicating sequences are not identical
            output.write("Differences:\n")
            # Compare sequences and character by character and write position
            # differences
            for i, (hfe_char, refseq_char) in enumerate(zip(hfe_sequence,
            refseq_sequence), start=1):
                if hfe_char != refseq_char:
                    output.write(f"Position {i}: "
                                 f"HFE_ORF_sequence has '{hfe_char}', "
                                 f"HFE_RefSeq has '{refseq_char}'\n") # Write
                    # position and difference if characters at the same
                    # position are different
            # Handle length differences
            if len(hfe_sequence) > len(refseq_sequence):
                output.write(f"HFE ORF sequence is longer. Extra residues: "
                             f"{hfe_sequence[len(refseq_sequence):]}\n")
                # Write extra residues in the HFE sequence if it is longer
            elif len(refseq_sequence) > len(hfe_sequence):
                output.write(f"HFE RefSeq sequence is longer. Extra residues: "
                             f"{refseq_sequence[len(hfe_sequence):]}\n")
                # Write extra residues in the RefSeq sequence if it is longer


def main():
    try:
        # Read HFE ORF sequence
        if not os.path.exists(protein_orf_file): # Check if
            # "2d_HFE_ORF_sequence.fasta" ORF file exists
            raise FileNotFoundError(f"File not found: {protein_orf_file}")
            # Raise error if absent
        hfe_sequence = read_protein_sequence(protein_orf_file) # Read ORF
        # sequence from file (Call read_protein_sequence function)

        # Extract RefSeq protein sequence
        if not os.path.exists(protein_file): # Check if "Protein.faa" RefSeq
            # file exists
            raise FileNotFoundError(f"File not found: {protein_file}")
            # Raise error if absent
        refseq_sequence = extract_refseq_sequence(protein_file,
                                                  "NP_000401.1") # Extract
        # RefSeq sequence for the given protein reference ID (Call
        # extract_refseq_sequence function)

        # Compare sequences and write results
        compare_sequences(hfe_sequence, refseq_sequence, output_file)
        # Compare sequences and save the results (Call compare_sequences
        # function)
        print(f"Comparison complete. Results successfully saved to"
              f" {output_file}")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()