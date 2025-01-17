import os # Import os library to be able to work with file paths and conduct
# system-related operations

# Define directory and file paths
script_dir = os.path.dirname(os.path.abspath(__file__)) # Initiates to
# directory where script is
task_2_dir = os.path.join(script_dir, "Task 2") # Path to "Task 2" folder
# within same directory as script

# File names within the Task 2 folder
reference_genome_file = os.path.join(task_2_dir, "Reference genome.fna")
# File path for "Reference genome"
annotations_file = os.path.join(task_2_dir, "Annotations.gtf") # File path
# for "Annotations"
output = os.path.join(task_2_dir, "2b_HFE_gene.fasta") # Output file path
# and name

# Define GeneID and TranscriptID
gene_id = "3077" # GeneID for HFE
transcript_id = "NM_000410.4" # TranscriptID for HFE

def extract_genomic_location(gtf_file, gene_id, transcript_id):
    """
    Extracts the genomic location of the HFE gene using both its GeneID and
    TranscriptID from the "Annotations" GTF file.

    Args:
        gtf_file (str): Path to GTF file.
        gene_id (str): HFE GeneID.
        transcript_id (str): HFE TranscriptID.

    Returns:
        dict: Dictionary containing HFE gene information (chromosome number,
        start location, end location, and strand (+/-)).
    """
    with open(gtf_file, "r") as f: # Open GTF file for reading ("r")
        for line in f: # Iterate through each line in the file
            if line.startswith("#"): # Skip comment lines
                continue
            fields = line.strip().split("\t") # Split each line into fields
            # using tabs
            if len(fields) < 9: # Skip lines with less than 9 fields
                # (insufficient fields - invalid format)
                continue

            attributes = fields[8] # Extract attributes field (9th column)

            # Filter for specific GeneID and TranscriptID
            if f"GeneID:{gene_id}" in attributes and f'transcript_id "{
            transcript_id}"' in attributes: # Check if GeneID and
                # TranscriptID match desired values
                # Extract location details
                chrom = fields[0] # Chromosome number
                start = int(fields[3]) # Start location/position
                end = int(fields[4]) # End location/position
                strand = fields[6] # Strand information (+/-)
                return {"chrom": chrom, "start": start, "end": end,
                        "strand": strand} # HFE gene information dictionary

    raise ValueError(
        f"GeneID {gene_id} or Transcript ID {transcript_id} not "
        f"found in the GTF file.") # Raise error if GeneID and/or
    # TranscriptID are not found

def parse_fasta(file_path):
    """
    Parses FASTA file into dictionary that maps sequence IDs to sequences.

    Args:
        file_path (str): Path to FASTA file.

    Returns:
        dict: Dictionary containing sequence IDs as keys and sequences as
        values.
    """
    sequences = {} # Empty dictionary to store sequences
    with open(file_path, "r") as f: # Open FASTA file for reading
        current_id = None # Store current sequence ID
        current_seq = [] # Empty list to store sequence lines
        for line in f: # Iterate through each line in the file
            if line.startswith(">"): # Detect start of new sequence
                if current_id is not None: # Save previous sequence if
                    # sequence ID exists
                    sequences[current_id] = "".join(current_seq) # Join
                    # sequence lines
                current_id = line[1:].strip().split()[0] # Extract sequence ID
                current_seq = [] # Reset sequence list
            else:
                current_seq.append(line.strip()) # Add sequence line to list
        if current_id is not None: # Add last sequence after loop
            sequences[current_id] = "".join(current_seq)
    return sequences # Return sequence dictionary


def extract_sequence(reference_file, chrom, start, end, strand):
    """
    Extracts sequence from specified region within reference genome.

    Args:
        reference_file (str): File path for "Reference genome".
        chrom (str): Chromosome number.
        start (int): Start location.
        end (int): End location.
        strand (str): Strand information (+/-).

    Returns:
        str: Extracted HFE gene sequence.
    """
    sequences = parse_fasta(reference_file) # Parse "Reference genome" file
    # into dictionary (Call parse_fasta function)
    if chrom not in sequences: # Check if specified chromosome number exists
        # in genome
        raise ValueError(f"Chromosome {chrom} not found in the reference "
                         f"genome.") # Raise error if chromosome number is not
        # found
    sequence = sequences[chrom][start - 1:end]  # Extract sequence
    # (initially convert from 1-based to 0-based indexing)
    if strand == "-": # Compute reverse complement if strand is negative (-)
        complement = str.maketrans("ACGTacgt", "TGCAtgca") # Complement table
        return sequence.translate(complement)[::-1] # Reverse complement
        # sequence and return it
    return sequence # Return forward sequence (+ strand)


def main():
    try:
        # Step 1: Extracts HFE gene location using "Annotations" file and,
        # geneID and transcriptID information
        gene_location = extract_genomic_location(annotations_file, gene_id,
                                                 transcript_id) # Call
        # extract_genomic_location function
        print(f"Genomic location of HFE gene: {gene_location}") # Print gene
        # location details (chromosome number, start and end locations and,
        # strand)

        # Step 2: Extracts sequence from reference genome using "Reference
        # genome" file and location details extracted when calling
        # extract_genomic_location function
        hfe_sequence = extract_sequence(
            reference_genome_file,
            gene_location["chrom"],
            gene_location["start"],
            gene_location["end"],
            gene_location["strand"],
        ) # Call extract_sequence function

        # Step 3: Writes extracted sequence into FASTA file
        with open(output, "w") as f: # Open "output" file
            f.write(f">HFE_gene|GeneID:{gene_id}|TranscriptID:{transcript_id}|"
                f"{gene_location['chrom']}:{gene_location['start']}-"
                f"{gene_location['end']}({gene_location['strand']})\n")
            # Write header including gene information and extraction details
            f.write(hfe_sequence.strip() + "\n") # Write extracted gene
            # sequence

        print(f"HFE gene sequence successfully saved to {output}") # Print
        # "output" file location

    except Exception as e:
        print(f"Error: {e}") # Print error if encountered


if __name__ == "__main__":
    main() # Call main function