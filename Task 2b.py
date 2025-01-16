import os

# Define the directory and file paths
script_dir = os.path.dirname(os.path.abspath(__file__))
task_2b_dir = os.path.join(script_dir, "Task 2")

# File names within the Task_2b folder
REFERENCE_GENOME_FILE = os.path.join(task_2b_dir, "Reference genome.fna")
ANNOTATION_FILE = os.path.join(task_2b_dir, "Annotations.gtf")
OUTPUT_FILE = os.path.join(task_2b_dir, "HFE_gene.fasta")

# GeneID and Transcript ID for the HFE gene
HFE_GENE_ID = "3077"
HFE_TRANSCRIPT_ID = "NM_000410.4"


def extract_genomic_location(gtf_file, gene_id, transcript_id):
    """
    Extract the genomic location of the HFE gene using both Gene ID and Transcript ID from the GTF file.

    Args:
        gtf_file (str): Path to the GTF file.
        gene_id (str): GeneID of the HFE gene.
        transcript_id (str): Transcript ID of the HFE gene.

    Returns:
        dict: A dictionary containing chromosome, start, end, and strand of the gene.
    """
    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            attributes = fields[8]

            # Filter for the specific Gene ID and Transcript ID
            if f"GeneID:{gene_id}" in attributes and f'transcript_id "{transcript_id}"' in attributes:
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                return {"chrom": chrom, "start": start, "end": end, "strand": strand}

    raise ValueError(
        f"GeneID {gene_id} or Transcript ID {transcript_id} not found in the GTF file."
    )


def parse_fasta(file_path):
    """
    Parse a FASTA file into a dictionary mapping sequence IDs to sequences.

    Args:
        file_path (str): Path to the FASTA file.

    Returns:
        dict: A dictionary where keys are sequence IDs and values are sequences.
    """
    sequences = {}
    with open(file_path, "r") as f:
        current_id = None
        current_seq = []
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].strip().split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id is not None:
            sequences[current_id] = "".join(current_seq)
    return sequences


def extract_sequence(reference_file, chrom, start, end, strand):
    """
    Extract the sequence of the specified region from the reference genome.

    Args:
        reference_file (str): Path to the reference genome file.
        chrom (str): Chromosome name.
        start (int): Start position (1-based).
        end (int): End position (1-based).
        strand (str): Strand information ('+' or '-').

    Returns:
        str: The extracted sequence (including lowercase letters).
    """
    sequences = parse_fasta(reference_file)
    if chrom not in sequences:
        raise ValueError(f"Chromosome {chrom} not found in the reference genome.")
    sequence = sequences[chrom][start - 1:end]  # Convert to 0-based indexing
    if strand == "-":
        complement = str.maketrans("ACGTacgt", "TGCAtgca")
        return sequence.translate(complement)[::-1]
    return sequence


def main():
    try:
        # Step 1: Extract genomic location of the HFE gene
        gene_location = extract_genomic_location(
            ANNOTATION_FILE, HFE_GENE_ID, HFE_TRANSCRIPT_ID
        )
        print(f"Genomic location of HFE gene: {gene_location}")

        # Step 2: Extract the sequence from the reference genome
        hfe_sequence = extract_sequence(
            REFERENCE_GENOME_FILE,
            gene_location["chrom"],
            gene_location["start"],
            gene_location["end"],
            gene_location["strand"],
        )

        # Step 3: Write the sequence to a FASTA file
        with open(OUTPUT_FILE, "w") as f:
            f.write(
                f">HFE_gene|GeneID:{HFE_GENE_ID}|TranscriptID:{HFE_TRANSCRIPT_ID}|{gene_location['chrom']}:{gene_location['start']}-{gene_location['end']}({gene_location['strand']})\n"
            )
            f.write(hfe_sequence.strip() + "\n")

        print(f"HFE gene sequence saved to {OUTPUT_FILE}")

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()