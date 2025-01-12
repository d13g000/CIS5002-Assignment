import os

# Define file paths
base_dir = os.path.dirname(os.path.abspath(__file__))  # Current directory
task_2_dir = os.path.join(base_dir, "Task 2")
protein_faa_file = os.path.join(task_2_dir, "protein.faa")
hfe_protein_frame_dir = os.path.join(task_2_dir, "HFE_protein_frames")
hfe_protein_frame_1_file = os.path.join(hfe_protein_frame_dir,
                                        "HFE_ORF_sequence.fasta")
output_file = os.path.join(task_2_dir, "HFE_protein_comparison.txt")

# Function to read a protein sequence from a file
def read_protein_sequence(file_path):
    sequence = []
    with open(file_path, "r") as file:
        for line in file:
            if not line.startswith(">"):  # Skip header lines
                sequence.append(line.strip())  # Add sequence lines
    return "".join(sequence)

# Read the HFE protein frame 1 sequence
if os.path.exists(hfe_protein_frame_1_file):
    hfe_sequence = read_protein_sequence(hfe_protein_frame_1_file)
else:
    raise FileNotFoundError(f"File not found: {hfe_protein_frame_1_file}")

# Read the RefSeq protein file and extract the NP_000401.1 sequence
if os.path.exists(protein_faa_file):
    refseq_sequence = None
    with open(protein_faa_file, "r") as file:
        lines = file.readlines()
        capturing = False
        for line in lines:
            if line.startswith(">"):
                if "NP_000401.1" in line:
                    capturing = True
                    refseq_sequence = []
                else:
                    capturing = False
            elif capturing:
                refseq_sequence.append(line.strip())
    if refseq_sequence is not None:
        refseq_sequence = "".join(refseq_sequence)
    else:
        raise ValueError("NP_000401.1 sequence not found in protein.faa file")
else:
    raise FileNotFoundError(f"File not found: {protein_faa_file}")

# Compare sequences and write the results to an output file
with open(output_file, "w") as output:
    if hfe_sequence == refseq_sequence:
        output.write("The sequences are identical.\n")
    else:
        output.write("The sequences are different.\n")
        output.write("Differences:\n")
        # Identify and write differences
        for i, (hfe_char, refseq_char) in enumerate(zip(hfe_sequence, refseq_sequence), start=1):
            if hfe_char != refseq_char:
                output.write(f"Position {i}: HFE_ORF_sequence has "
                             f"'{hfe_char}', HFE_RefSeq has '{refseq_char}'\n")
        # Handle length differences
        if len(hfe_sequence) > len(refseq_sequence):
            output.write(f"HFE_ORF_sequence is longer. Extra residues:"
                         f" {hfe_sequence[len(refseq_sequence):]}\n")
        elif len(refseq_sequence) > len(hfe_sequence):
            output.write(f"HFE_RefSeq sequence is longer. Extra residues:"
                         f" {refseq_sequence[len(hfe_sequence):]}\n")

print(f"Comparison complete. Results saved to {output_file}")