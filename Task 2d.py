import os

# Genetic Code Table for Translation
GENETIC_CODE = {
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
}

def translate_sequence(sequence):
    """
    Translates an mRNA sequence into a protein sequence.
    Stop codons (*) are included in the translation.

    Args:
        sequence (str): The mRNA sequence.

    Returns:
        str: The translated protein sequence.
    """
    protein = []
    for i in range(0, len(sequence) - 2, 3):  # Process in codons (triplets)
        codon = sequence[i:i+3]
        amino_acid = GENETIC_CODE.get(codon, '')
        if amino_acid:  # Skip invalid codons
            protein.append(amino_acid)
    return ''.join(protein)

def translate_mrna(mrna_sequence):
    """
    Translates the entire mRNA sequence into protein sequences in all three reading frames.

    Args:
        mrna_sequence (str): The mRNA sequence.

    Returns:
        dict: A dictionary with frame numbers (1, 2, 3) and the corresponding protein sequences.
    """
    # Filter to only keep uppercase letters (valid mRNA nucleotides)
    mrna_sequence = ''.join([base for base in mrna_sequence if base.isupper()])

    protein_sequences = {}

    for frame in range(3):
        # Extract the sequence for the current reading frame
        frame_sequence = mrna_sequence[frame:]
        # Translate the sequence
        protein_sequence = translate_sequence(frame_sequence)
        protein_sequences[frame + 1] = protein_sequence

    return protein_sequences

def find_orf_sequence(protein_sequence):
    """
    Extracts the ORF from a given protein sequence.
    The ORF starts with 'M' (methionine) and ends at the first '*' (stop codon).
    This version excludes the stop codon.

    Args:
        protein_sequence (str): The full protein sequence.

    Returns:
        str: The ORF sequence (from 'M' to the first codon before '*').
    """
    start_index = protein_sequence.find('M')
    if start_index == -1:
        return ""  # No start codon found

    stop_index = protein_sequence.find('*', start_index)
    if stop_index == -1:
        return protein_sequence[start_index:]  # No stop codon found, return from start to the end

    return protein_sequence[start_index:stop_index]  # Exclude the stop codon

def find_orfs(protein_sequences):
    """
    Outputs the ORF for each reading frame and finds the correct ORF, based on the longest protein sequence.

    Args:
        protein_sequences (dict): Dictionary of protein sequences for each reading frame.

    Returns:
        tuple: The correct ORF protein sequence and its reading frame.
    """
    longest_orf = ""
    longest_frame = -1

    for frame, protein_sequence in protein_sequences.items():
        # Extract the ORF from the protein sequence
        orf = find_orf_sequence(protein_sequence)
        # Identify the longest ORF
        if len(orf) > len(longest_orf):
            longest_orf = orf
            longest_frame = frame

    return longest_orf, longest_frame

def translate_mrna_to_orf(mrna_file, output_folder):
    """
    Translates the mRNA sequence in all three frames and finds the correct ORF.

    Args:
        mrna_file (str): Path to the input FASTA file containing the mRNA sequence.
        output_folder (str): Path to the folder where the separate protein sequences will be saved.
    """
    try:
        # Read the mRNA sequence from the FASTA file
        with open(mrna_file, "r") as f:
            lines = f.readlines()
            header = lines[0].strip()  # FASTA header
            mrna_sequence = "".join(line.strip() for line in lines[1:])  # Join mRNA sequence lines

        # Translate mRNA to protein sequences in all frames
        protein_sequences = translate_mrna(mrna_sequence)

        # Find the valid ORF
        longest_orf, longest_frame = find_orfs(protein_sequences)

        # Print the longest ORF and its frame
        print(f"The ORF is found in Frame {longest_frame}")
        print(f"ORF Protein Sequence:\n{longest_orf}")

        # Create output folder if it doesn't exist
        os.makedirs(output_folder, exist_ok=True)

        # Write the protein sequences for each frame to separate FASTA files
        for frame, protein_sequence in protein_sequences.items():
            output_file = os.path.join(output_folder, f"HFE_protein_frame_{frame}.fasta")
            with open(output_file, "w") as f:
                f.write(f"{header} (ORF: Frame {longest_frame})\n")
                f.write(protein_sequence + "\n")
            print(f"Protein sequence for Frame {frame} saved to {output_file}")

        # Save the longest ORF to a separate file
        orf_output_file = os.path.join(output_folder, "HFE_ORF_sequence.fasta")
        with open(orf_output_file, "w") as f:
            f.write(f"{header} (ORF Sequence)\n")
            f.write(longest_orf + "\n")
        print(f"ORF sequence saved to {orf_output_file}")

    except Exception as e:
        print(f"Error in translation: {e}")

def main():
    # Define paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    task_2b_dir = os.path.join(script_dir, "Task 2")
    mrna_file = os.path.join(task_2b_dir, "HFE_mrna.fasta")
    output_folder = os.path.join(task_2b_dir, "HFE_protein_frames")

    # Translate mRNA to find the ORF and output the sequences
    translate_mrna_to_orf(mrna_file, output_folder)

if __name__ == "__main__":
    main()