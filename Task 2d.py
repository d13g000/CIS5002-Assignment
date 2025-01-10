import os
from Bio.Seq import Seq


def translate_mrna(mrna_sequence):
    """
    Translates the entire mRNA sequence into protein sequences in all three reading frames.
    This ensures only uppercase letters are considered for translation.

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
        protein_sequence = Seq(frame_sequence).translate()  # Translate the whole sequence including stop codons
        protein_sequences[frame + 1] = str(protein_sequence)

    return protein_sequences


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
        # Identify the longest ORF
        if len(protein_sequence) > len(longest_orf):
            longest_orf = protein_sequence
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
            mrna_sequence = "".join(
                line.strip() for line in lines[1:])  # Join mRNA sequence lines

        # Translate mRNA to protein sequences in all frames
        protein_sequences = translate_mrna(mrna_sequence)

        # Find the valid ORF
        longest_orf, longest_frame = find_orfs(protein_sequences)

        # Create output folder if it doesn't exist
        os.makedirs(output_folder, exist_ok=True)

        # Write the protein sequences for each frame to separate FASTA files
        for frame, protein_sequence in protein_sequences.items():
            output_file = os.path.join(output_folder,
                                       f"HFE_protein_frame_{frame}.fasta")
            with open(output_file, "w") as f:
                f.write(f"{header} (ORF: Frame {longest_frame})\n")
                f.write(protein_sequence + "\n")
            print(f"Protein sequence for Frame {frame} saved to {output_file}")

        # Print the frame that contains the ORF (longest protein sequence)
        print(f"The ORF is found in Frame {longest_frame}")

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
