import os
from Bio import SeqIO
from Bio.Seq import Seq

# Define file paths dynamically
script_dir = os.path.dirname(os.path.abspath(__file__))
task_2b_dir = os.path.join(script_dir, "Task 2")
REFERENCE_GENOME_FILE = os.path.join(task_2b_dir, "Reference genome.fna")
ANNOTATION_FILE = os.path.join(task_2b_dir, "Annotations.gtf")
HFE_GENE_SEQUENCE_FILE = os.path.join(task_2b_dir, "HFE_gene.fasta")
OUTPUT_DIR = os.path.join(task_2b_dir, "HFE_mrna_variants")

# Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def transcribe_dna_to_mrna(dna_sequence_file, output_mrna_file):
    """
    Transcribes a DNA sequence into mRNA and saves it to a file.

    Args:
        dna_sequence_file (str): Path to the input FASTA file containing the DNA sequence.
        output_mrna_file (str): Path to the output file to save the mRNA sequence.
    """
    try:
        # Read the DNA sequence from the FASTA file
        with open(dna_sequence_file, "r") as f:
            lines = f.readlines()
            header = lines[0].strip()  # FASTA header
            dna_sequence = "".join(
                line.strip() for line in lines[1:])  # Join DNA sequence lines

        # Transcribe DNA to mRNA
        dna_seq = Seq(dna_sequence)
        mrna_sequence = dna_seq.transcribe()

        # Write the mRNA sequence to the output file
        with open(output_mrna_file, "w") as f:
            f.write(f"{header}\n")  # Reuse the FASTA header
            f.write(str(mrna_sequence) + "\n")

        print(f"mRNA sequence saved to {output_mrna_file}")

    except Exception as e:
        print(f"Error in transcription: {e}")

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
    with open(reference_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            if record.id == chrom:
                # Concatenate exon sequences and filter out lowercase letters
                sequence = "".join([str(record.seq[start - 1:end]).upper() for start, end in sorted(exons)])
                # Return reverse complement if on negative strand
                return str(Seq(sequence).reverse_complement()) if strand == "-" else sequence
    raise ValueError(f"Chromosome {chrom} not found in the reference genome.")

def main():
    # Define the GeneID for which to extract mRNA variants
    gene_id = "3077"  # HFE GeneID

    try:
        # Step 1: Transcribe DNA to mRNA (using HFE gene sequence)
        output_mrna_file = os.path.join(task_2b_dir, "HFE_mrna.fasta")
        transcribe_dna_to_mrna(HFE_GENE_SEQUENCE_FILE, output_mrna_file)

        # Step 2: Parse the GTF file to get transcript variants
        transcripts = parse_gtf_for_transcripts(ANNOTATION_FILE, gene_id)
        print(f"Found {len(transcripts)} transcript variants for GeneID {gene_id}.")

        # Step 3: Extract sequences for each transcript and save to files
        for transcript_id, data in transcripts.items():
            chrom, strand, exons = data["chrom"], data["strand"], data["exons"]
            mrna_sequence = extract_exon_sequences(REFERENCE_GENOME_FILE, chrom, exons, strand)

            # Write each transcript's mRNA to a separate FASTA file
            output_file = os.path.join(OUTPUT_DIR, f"{transcript_id}_mRNA.fasta")
            with open(output_file, "w") as f:
                f.write(f">{transcript_id}|GeneID:{gene_id}|{chrom}|{strand}\n")
                f.write(mrna_sequence + "\n")

            print(f"Saved mRNA sequence for transcript {transcript_id} to {output_file}")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
