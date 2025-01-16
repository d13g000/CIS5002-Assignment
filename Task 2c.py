import os

# Define file paths dynamically
script_dir = os.path.dirname(os.path.abspath(__file__))
task_2b_dir = os.path.join(script_dir, "Task 2")
REFERENCE_GENOME_FILE = os.path.join(task_2b_dir, "Reference genome.fna")
ANNOTATION_FILE = os.path.join(task_2b_dir, "Annotations.gtf")
HFE_GENE_SEQUENCE_FILE = os.path.join(task_2b_dir, "HFE_gene.fasta")
OUTPUT_DIR = os.path.join(task_2b_dir, "HFE_mrna_variants")

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

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

def transcribe_to_mrna(dna_sequence):
    """
    Transcribe a DNA sequence to mRNA by replacing T/t with U/u.

    Args:
        dna_sequence (str): Input DNA sequence.

    Returns:
        str: Transcribed mRNA sequence.
    """
    return dna_sequence.replace("T", "U").replace("t", "u")

def parse_gtf_for_transcripts(gtf_file, gene_id):
    """
    Parse the GTF file to extract exons for all transcript variants of a given gene.

    Args:
        gtf_file (str): Path to the GTF file.
        gene_id (str): GeneID of the target gene.

    Returns:
        dict: A dictionary with transcript IDs as keys and their exon locations as values.
              Each value is a list of tuples (chrom, start, end, strand).
    """
    transcripts = {}
    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            chrom, feature, start, end, strand, attributes = fields[0], fields[2], int(fields[3]), int(fields[4]), fields[6], fields[8]
            if feature == "exon" and f"GeneID:{gene_id}" in attributes:
                transcript_id = None
                for attr in attributes.split(";"):
                    if "transcript_id" in attr:
                        transcript_id = attr.split('"')[1]
                        break
                if transcript_id:
                    if transcript_id not in transcripts:
                        transcripts[transcript_id] = {"chrom": chrom, "strand": strand, "exons": []}
                    transcripts[transcript_id]["exons"].append((start, end))
    return transcripts

def extract_exon_sequences(reference_file, chrom, exons, strand):
    """
    Extract the combined sequence for a list of exons from the reference genome.

    Args:
        reference_file (str): Path to the reference genome file.
        chrom (str): Chromosome name.
        exons (list): List of tuples representing exon start and end positions.
        strand (str): Strand information ('+' or '-').

    Returns:
        str: The combined mRNA sequence.
    """
    sequences = parse_fasta(reference_file)
    if chrom not in sequences:
        raise ValueError(f"Chromosome {chrom} not found in the reference genome.")
    sequence = "".join(sequences[chrom][start - 1:end].upper() for start, end in sorted(exons))
    if strand == "-":
        complement = str.maketrans("ACGTacgt", "UGCAugca")
        sequence = sequence.translate(complement)[::-1]  # Reverse complement for negative strand
    return sequence

def convert_hfe_gene_to_mrna(hfe_gene_file, output_file):
    """
    Convert the HFE gene sequence file into an mRNA file.

    Args:
        hfe_gene_file (str): Path to the HFE gene DNA sequence file.
        output_file (str): Path to save the mRNA sequence.
    """
    sequences = parse_fasta(hfe_gene_file)
    for seq_id, dna_sequence in sequences.items():
        mrna_sequence = transcribe_to_mrna(dna_sequence)
        with open(output_file, "w") as f:
            f.write(f">{seq_id}\n")
            f.write(mrna_sequence + "\n")
    print(f"Converted HFE gene sequence to mRNA and saved to {output_file}")

def main():
    # Define the GeneID for which to extract mRNA variants
    gene_id = "3077"  # HFE GeneID

    try:
        # Step 1: Convert the HFE gene sequence file into mRNA
        hfe_mrna_file = os.path.join(task_2b_dir, "HFE_gene_mRNA.fasta")
        convert_hfe_gene_to_mrna(HFE_GENE_SEQUENCE_FILE, hfe_mrna_file)

        # Step 2: Parse the GTF file to get transcript variants
        transcripts = parse_gtf_for_transcripts(ANNOTATION_FILE, gene_id)
        print(f"Found {len(transcripts)} transcript variants for GeneID {gene_id}.")

        # Step 3: Extract sequences for each transcript and save to files
        for transcript_id, data in transcripts.items():
            chrom, strand, exons = data["chrom"], data["strand"], data["exons"]
            dna_sequence = extract_exon_sequences(REFERENCE_GENOME_FILE, chrom, exons, strand)

            # Transcribe DNA to mRNA
            mrna_sequence = transcribe_to_mrna(dna_sequence)

            # Save the mRNA sequence to a FASTA file
            output_file = os.path.join(OUTPUT_DIR, f"{transcript_id}_mRNA.fasta")
            with open(output_file, "w") as f:
                f.write(f">{transcript_id}|GeneID:{gene_id}|{chrom}|{strand}\n")
                f.write(mrna_sequence.strip() + "\n")

            print(f"Saved mRNA sequence for transcript {transcript_id} to {output_file}")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()