import os
from Bio import SeqIO
import re

# Define paths to the files
folder = "reference_sequence_and_annotation"
annotation_file = os.path.join(folder, "annotations.gtf")
genome_file = os.path.join(folder, "genome.fasta")

# Define the HFE gene details
gene_id = "3077"  # Numeric part of GeneID:3077

def parse_gtf(gtf_file, target_gene_id):
    """
    Parses a GTF file to extract exons for the specified gene ID based on db_xref.
    """
    exons = []
    target_pattern = f"GeneID:{target_gene_id}"  # Match GeneID:3077 in db_xref

    with open(gtf_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue

            # Check if db_xref contains the target GeneID
            attributes = fields[8]
            if re.search(target_pattern, attributes):
                # Debugging: Print matched attributes
                print(f"Matched: {attributes}")

                exons.append({
                    "seqname": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6]
                })

    if not exons:
        print("Debugging: No exons found matching the specified GeneID.")
    return exons

def extract_sequence(genome_file, exons):
    """
    Extracts sequences for the provided exons from the genome file.
    """
    # Load the genome reference sequence into a dictionary
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    # Extract exon sequences
    sequences = []
    for exon in exons:
        seq_record = genome[exon["seqname"]]
        sequence = seq_record.seq[exon["start"] - 1:exon["end"]]
        if exon["strand"] == "-":
            sequence = sequence.reverse_complement()
        sequences.append(sequence)

    # Concatenate all exon sequences
    return "".join(str(seq) for seq in sequences)

def main():
    # Parse the GTF file to extract exons for the HFE gene
    exons = parse_gtf(annotation_file, gene_id)

    if not exons:
        print("No exons found for the HFE gene.")
        return  # Exit the function gracefully

    print(f"Found {len(exons)} exons for the HFE gene.")  # Debugging: Number of exons

    # Extract the sequence from the genome file
    hfe_sequence = extract_sequence(genome_file, exons)

    # Define the output file path
    output_file = os.path.join(folder, "HFE_gene_sequence.fasta")

    # Save the extracted sequence to the file
    with open(output_file, "w") as file:
        file.write(">HFE_gene_sequence\n")
        file.write(hfe_sequence)

    print(f"HFE gene sequence extracted and saved to {output_file}.")

if __name__ == "__main__":
    main()
