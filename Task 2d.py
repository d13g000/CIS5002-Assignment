import os

# Define directory and file paths
script_dir = os.path.dirname(os.path.abspath(__file__))
task_2_dir = os.path.join(script_dir, "Task 2")

mrna_file = os.path.join(task_2_dir, "2c_HFE_gene_mRNA.fasta")  # File path
# for mRNA sequence generated in script for Task 2c
output = os.path.join(task_2_dir, "2d_HFE_protein_reading_frames")
# Output folder path and name

# Translation table
translation_table = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}  # Dictionary mapping mRNA codons to their corresponding amino acids


def translate_sequence(sequence):
    """
    Translates mRNA sequence into protein sequence.

    Args:
        sequence (str): mRNA sequence.

    Returns:
        str: Translated protein sequence.
    """
    protein = []  # Empty list to store protein sequence
    for i in range(0, len(sequence) - 2, 3):  # Iterate over mRNA sequence in
        # codons (triplets)
        codon = sequence[i:i + 3]  # Extract codon (3 nucleotides)
        amino_acid = translation_table.get(codon, '')  # Get corresponding
        # amino acid for extracted codon
        if amino_acid:  # Skip invalid codons
            protein.append(amino_acid)  # Add amino acid to protein sequence
            # list
    return ''.join(protein)  # Return protein sequence as string of joined
    # amino acids


def reverse_complement(sequence):
    """
    Computes reverse complement of mRNA sequence.

    Args:
        sequence (str): mRNA sequence.

    Returns:
        str: Reverse complement of mRNA sequence.
    """
    complement = str.maketrans("AUCG", "UAGC")  # Complement table
    return sequence.translate(complement)[::-1]  # Return translated and
    # reversed sequence


def translate_mrna(mrna_sequence):
    """
    Translates mRNA sequence in all six reading frames.

    Args:
        mrna_sequence (str): mRNA sequence.

    Returns:
        dict: Dictionary with frame numbers (1-6) and their corresponding
        protein sequences.
    """
    protein_sequences = {}  # Empty dictionary to store protein sequences from
    # each reading frame

    # Forward frames (1, 2, 3)
    for frame in range(3):
        frame_sequence = mrna_sequence[frame:]
        protein_sequences[frame + 1] = translate_sequence(frame_sequence)
        # Offset sequence for each frame by single nucleotide and translate
        # it (Call translate_sequence function)

    # Reverse complement frames (4, 5, 6)
    reverse_comp_sequence = reverse_complement(mrna_sequence)  # Compute
    # reverse complement (Call reverse complement function)
    for frame in range(3):
        frame_sequence = reverse_comp_sequence[frame:]
        protein_sequences[frame + 4] = translate_sequence(frame_sequence)

    return protein_sequences  # Dictionary of protein sequences generated from
    # each frame


def find_orf_sequence(protein_sequence):
    """
    Extracts open reading frame (ORF) from protein sequence where ORF starts
    with 'M' (methionine) and ends at the first '*' (stop codon).

    Args:
        protein_sequence (str): Protein sequence.

    Returns:
        str: ORF sequence (from 'M' to the first amino acid before '*').
    """
    start_index = protein_sequence.find('M')  # Find first 'M' (start codon)
    if start_index == -1:  # No start codon found
        return ""  # Return empty string

    stop_index = protein_sequence.find('*', start_index)  # Find first '*'
    # (stop codon) after first start codon
    if stop_index == -1:  # No stop codon found
        return protein_sequence[start_index:]  # Return protein sequence from
        # start to the end

    return protein_sequence[start_index:stop_index]  # Return ORF
    # sequence excluding stop codon


def find_orfs(protein_sequences):
    """
    Finds longest ORF among all reading frames.

    Args:
        protein_sequences (dict): Dictionary of protein sequences produced
        using each reading frame.

    Returns:
        tuple: Longest ORF sequence and its respective reading frame.
    """
    longest_orf = ""  # Empty string to store ORF
    longest_frame = -1  # Initialise frame number

    for frame, protein_sequence in protein_sequences.items():  # Iterate
        # through each reading frame
        orf = find_orf_sequence(protein_sequence)  # Extract ORF for each
        # frame (Call find_orf_sequence function)
        if len(orf) > len(longest_orf):  # Check if ORF extracted is longest
            longest_orf = orf  # Update ORF
            longest_frame = frame  # Update frame number

    return longest_orf, longest_frame  # Return longest ORF and its respective
    # frame


def translate_mrna_to_orf(mrna_file, output_folder):
    """
    Translates mRNA sequence for all six frames and finds ORF.

    Args:
        mrna_file (str): Path to FASTA file containing gene mRNA sequence (
        "2c_HFE_gene_mrna.fasta").
        output_folder (str): Path to folder where the protein sequences will
        be saved ("2d_HFE_protein_reading_frames").
    """
    try:
        # Read mRNA sequence from the FASTA file
        with open(mrna_file, "r") as f:
            lines = f.readlines()  # Read all lines from
            # "2c_HFE_gene_mrna.fasta" FASTA file
            header = lines[0].strip()  # Extract FASTA file header
            mrna_sequence = "".join(line.strip() for line in lines[1:])
            # Combine mRNA sequence lines

        mrna_sequence = "".join([base for base in mrna_sequence if
                                 base.isupper()]) # Remove lowercase letters
        # (soft-masked/non-coding regions)

        # Translate mRNA to protein sequences in all six frames
        protein_sequences = translate_mrna(mrna_sequence)  # Call
        # translate_mrna function

        # Find longest ORF
        longest_orf, longest_frame = find_orfs(protein_sequences)  # Call
        # find_orfs function

        print(f"The ORF is found in Frame {longest_frame}")  # Print frame
        # belonging to ORF
        print(f"ORF Protein Sequence ({len(longest_orf)} aa): {longest_orf}")
        # Print ORF sequence and length in amino acids

        # Create "2d_HFE_protein_reading_frames" output folder if not present
        os.makedirs(output_folder, exist_ok=True)

        # Write protein sequences for each frame to separate FASTA files
        for frame, protein_sequence in protein_sequences.items():
            output_file = os.path.join(output_folder,
                f"2d_HFE_protein_reading_frame_{frame}.fasta")  # Name files
            # according to frame where 1, 2, 3 are forward frames and 4, 5,
            # 6 are reverse complement frames
            with open(output_file, "w") as f:
                f.write(f"{header} (Frame {frame})\n")  # Write header
                # including extracted FASTA file header and frame number
                f.write(protein_sequence + "\n")  # Write protein sequence
                # associated to frame
            print(f"Protein sequence for Frame {frame} successfully saved to "
                  f"{output_file}")  # Print "output_file" location

        # Save ORF to a separate file
        orf_output_file = os.path.join(output_folder,
                                       "2d_HFE_ORF_sequence.fasta")
        with open(orf_output_file, "w") as f:
            f.write(f"{header} (ORF: Frame {longest_frame})\n")
            # Write header including extracted FASTA file header and ORF frame
            # number
            f.write(longest_orf + "\n")  # Write ORF sequence
        print(f"ORF sequence successfully saved to {orf_output_file}")

    except Exception as e:
        print(f"Error in translation: {e}")


def main():
    # Step 1: Translate mRNA to protein and find ORF
    translate_mrna_to_orf(mrna_file, output)  # Call translate_mrna_to_orf
    # function


if __name__ == "__main__":
    main()