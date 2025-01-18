import os

# Define directory and file paths
script_dir = os.path.dirname(os.path.abspath(__file__))
task_2_dir = os.path.join(script_dir, "Task 2")

reference_genome_file = os.path.join(task_2_dir, "Reference genome.fna")
annotations_file = os.path.join(task_2_dir, "Annotations.gtf")
gene_sequence_file = os.path.join(task_2_dir, "2b_HFE_gene.fasta") # File path
# for gene sequence generated in script for Task 2b
output = os.path.join(task_2_dir, "2c_HFE_mRNA_variants") # Output file path
# and name

# Ensure output directory exists
os.makedirs(output, exist_ok=True) # Create output directory if non-existent

# Define GeneID
gene_id = "3077"  # GeneID for HFE

def parse_fasta(file_path):
    """
    Parses FASTA file into dictionary mapping sequence IDs to sequences.

    Args:
        file_path (str): Path to FASTA file.

    Returns:
        dict: Dictionary containing sequence IDs as keys and sequences as
        values.
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
    Transcribes DNA sequence to mRNA by replacing T/t with U/u.

    Args:
        dna_sequence (str): DNA sequence to transcribe.

    Returns:
        str: Transcribed mRNA sequence.
    """
    return dna_sequence.replace("T", "U").replace("t", "u") # Replace thymine
    # (T/t) with uracil (U/u)


def convert_gene_to_mrna(hfe_gene_file, output_file):
    """
    Converts "2b_HFE_gene.fasta" file into mRNA file ("2c_HFE_mrna_variants").

    Args:
        hfe_gene_file (str): Path to HFE gene DNA sequence file (
        "2b_HFE_gene.fasta").
        output_file (str): Output file path to save mRNA sequence (
        "2c_HFE_mrna_variants").
    """
    sequences = parse_fasta(hfe_gene_file) # Parse gene sequence file (Call
    # parse_fasta function)
    for seq_id, dna_sequence in sequences.items(): # Iterate over each sequence
        mrna_sequence = transcribe_to_mrna(dna_sequence) # Transcribe DNA to
        # mRNA (Call transcribe_to_mrna function)
        with open(output_file, "w") as f: # Open "2c_HFE_mrna_variants" file
            # for writing
            f.write(f">{seq_id}\n") # Write header including sequence ID
            f.write(mrna_sequence + "\n") # Write transcribed mRNA sequence
    print(f"Converted HFE gene sequence to mRNA successfully and saved to"
          f" {output_file}") # Print to notify successful transcription and
    # "output" file location


def parse_gtf_for_transcripts(gtf_file, gene_id):
    """
    Parses "Annotations.gtf" file to extract exons for all transcript
    variants of gene.

    Args:
        gtf_file (str): Path to GTF file.
        gene_id (str): GeneID of target gene.

    Returns:
        dict: Dictionary with transcript IDs as keys and their exon
        locations as values (each value is list of tuples (chromosome number,
        start location, end location and, strand).
    """
    transcripts = {} # Empty dictionary to store transcripts
    with (open(gtf_file, "r") as f): # Open GTF file for reading
        for line in f: # Iterate through each line in the file
            if line.startswith("#"): # Skip comment lines
                continue
            fields = line.strip().split("\t") # Split each line into fields
            # using tabs
            if len(fields) < 9: # Skip lines with less than 9 fields
                # (insufficient fields - invalid format)
                continue

            chrom, feature, start, end, strand, attributes = fields[0
            ], fields[2], int(fields[3]), int(fields[4]), fields[6], fields[
                8] # Extract relevant fields
            if feature == "exon" and f"GeneID:{gene_id}" in attributes: #
                # Check for exons of target gene
                transcript_id = None
                for attr in attributes.split(";"): # Extract transcript ID
                    # from attributes
                    if "transcript_id" in attr:
                        transcript_id = attr.split('"')[1] # Parse transcript ID
                        break
                if transcript_id:
                    if transcript_id not in transcripts: # Initialise
                        # transcript if not present
                        transcripts[transcript_id] = {"chrom": chrom,
                                                      "strand": strand,
                                                      "exons": []}
                    transcripts[transcript_id]["exons"].append((start,
                                                                end)) # Add
                    # exon start and end positions to list within dictionary
    return transcripts # Transcript information dictionary


def extract_exon_sequences(reference_file, chrom, exons, strand):
    """
    Extracts the combined sequence for a list of exons from the reference
    genome.

    Args:
        reference_file (str): Path to "Reference genome".
        chrom (str): Chromosome number.
        exons (list): Tuple representing exon start and end positions.
        strand (str): Strand information (+/-).

    Returns:
        str: mRNA sequence.
    """
    sequences = parse_fasta(reference_file) # Parse "Reference genome" (Calls
    # parse_fasta function)
    if chrom not in sequences:
        raise ValueError(f"Chromosome {chrom} not found in the reference "
                         f"genome.")
    # Combine exon sequences in order
    sequence = "".join(sequences[chrom][start - 1:end].upper() for start,
    end in sorted(exons))
    if strand == "-": # Handle reverse complement for for negative strand
        complement = str.maketrans("ACGTacgt", "UGCAugca")
        sequence = sequence.translate(complement)[::-1]  # Reverse
        # complement for negative strand
    return sequence # Return forward/reverse sequence


def main():
    try:
        # Step 1: Transcribes DNA gene sequence into mRNA
        hfe_mrna_file = os.path.join(task_2_dir, "2c_HFE_gene_mRNA.fasta")
        # File path for transcribed mRNA sequence
        convert_gene_to_mrna(gene_sequence_file, hfe_mrna_file) # Call
        # convert_gene_to_mrna function

        # Step 2: Parses "Annotations" GTF file to determine transcript variants
        transcripts = parse_gtf_for_transcripts(annotations_file, gene_id)
        # Call parse_gtf_for_transcripts function
        print(f"Found {len(transcripts)} transcript variants for GeneID "
              f"{gene_id}.") # Print number of transcripts found

        # Step 3: a) Extracts transcript sequences and saves them into FASTA
        # files
        for transcript_id, data in transcripts.items():
            chrom, strand, exons = data["chrom"], data["strand"], data["exons"]
            dna_sequence = extract_exon_sequences(reference_genome_file,
                                                  chrom, exons, strand)
            # Call extract_exon_sequences function

            # b) Transcribes DNA to mRNA
            mrna_sequence = transcribe_to_mrna(dna_sequence) # Call
            # transcribe_to_mrna function

            # c) Saves mRNA sequences to a FASTA files
            output_file = os.path.join(output, f"2c_{transcript_id}_mRNA.fasta")
            # Output file name and path (within "output" folder)
            with open(output_file, "w") as f: # Open transcript ID output file
                f.write(f">{transcript_id}|GeneID:{gene_id}|{chrom}|"
                        f"{strand}\n") # Write header including transcript ID
                # and gene information
                f.write(mrna_sequence.strip() + "\n") # Write extracted mRNA
                # sequence

            print(f"Saved mRNA sequence for transcript {transcript_id} to "
                  f"{output_file}") # Print "output_file" file location

    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()