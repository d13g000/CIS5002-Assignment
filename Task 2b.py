import os
from Bio import SeqIO
from Bio.Seq import Seq

# Define the directory and file paths
script_dir = os.path.dirname(os.path.abspath(__file__))
task_2b_dir = os.path.join(script_dir, "Task 2")

# File names within the Task_2b folder
REFERENCE_GENOME_FILE = os.path.join(task_2b_dir, "Reference genome.fna")
ANNOTATION_FILE = os.path.join(task_2b_dir, "Annotations.gtf")
OUTPUT_FILE = os.path.join(task_2b_dir, "HFE_gene_sequence.fasta")

# GeneID for the HFE gene
HFE_GENE_ID = "3077"

def extract_hfe_genomic_location(gtf_file, gene_id):
    """
    Extract the genomic location of the HFE gene from the GTF file.

    Args:
        gtf_file (str): Path to the GTF file.
        gene_id (str): GeneID of the HFE gene.

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
            if f"GeneID:{gene_id}" in attributes and fields[2] == "gene":
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                return {"chrom": chrom, "start": start, "end": end, "strand": strand}

    raise ValueError(f"GeneID {gene_id} not found in the GTF file.")

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
        str: The extracted sequence with lowercase letters removed.
    """
    with open(reference_file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            if record.id == chrom:
                sequence = record.seq[start - 1 : end]  # Convert to 0-based indexing
                filtered_sequence = "".join(
                    [char for char in str(sequence) if char.isupper()])
                if strand == "-":
                    return str(Seq(filtered_sequence).reverse_complement())
                else:
                    return filtered_sequence

    raise ValueError(f"Chromosome {chrom} not found in the reference genome.")

def main():
    try:
        # Step 1: Extract genomic location of the HFE gene
        gene_location = extract_hfe_genomic_location(ANNOTATION_FILE, HFE_GENE_ID)
        print(f"Genomic location of HFE gene: {gene_location}")

        # Step 2: Extract the sequence from the reference genome
        hfe_sequence = extract_sequence(
            REFERENCE_GENOME_FILE,
            gene_location["chrom"],
            gene_location["start"],
            gene_location["end"],
            gene_location["strand"]
        )

        # Step 3: Write the sequence to a FASTA file
        with open(OUTPUT_FILE, "w") as f:
            f.write(f">HFE_gene|GeneID:{HFE_GENE_ID}|{gene_location['chrom']}:{gene_location['start']}-{gene_location['end']}({gene_location['strand']})\n")
            f.write(hfe_sequence + "\n")

        print(f"HFE gene sequence saved to {OUTPUT_FILE}")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()