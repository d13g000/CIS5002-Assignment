import os  # Library for interacting with the operating system
import requests  # Library for sending HTTP requests to fetch data from web
# APIs
import subprocess  # Library for running external programs
from Bio import AlignIO  # Biopython module for reading, writing, and analysing
# sequence alignments

# Function to fetch UniProt IDs for different species of given protein
def fetch_protein_id(protein_name, species_name):
    """
    Constructs a query to the UniProt API to find the UniProt ID (accession
    number) for specific protein and species.

    Parameters:
        protein_name (str): Name of the protein ("Heme oxygenase 1").
        species_name (str): Name of the species.

    Returns:
        str: The UniProt ID of the protein for the species, or None if not
        found.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"  # Base URL for
    # querying UniProt
    query = (f"query={protein_name} AND (organism_name:"
             f"{species_name})&fields=accession") # Constructs query string
    # to search for the protein in the specific species
    url = f"{base_url}?{query}"  # Completes URL with the query included

    try:
        response = requests.get(url) # Sends GET request to UniProt API
        if response.status_code == 200: # Checks if request was successful
            data = response.json() # Parses JSON response into a Python
            # dictionary
            if "results" in data and data["results"]: # Checks if there are
                # results in the response
                return data["results"][0]["primaryAccession"] # Returns the
                # first result's primary accession (UniProt ID)
        else:
            print(f"Failed to fetch protein ID for {protein_name} in"
                  f" {species_name}.")
            print("Response:", response.text) # Prints error message if the
            # request fails
    except Exception as e: # Handles unexpected errors
        print(f"An error occurred while fetching protein ID: {e}")
    return None  # Returns None if no ID is found

# Function to download protein sequences from UniProt
def download_protein_sequence(protein_id, output_file, format="fasta"):
    """
    Fetches the protein sequence for a given UniProt ID and saves it to a
    file in the specified format.

    Parameters:
        protein_id (str): UniProt ID of the protein to download.
        output_file (str): File path where the sequence should be saved.
        format (str): Format for the sequence ("fasta").
    """
    base_url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.{format}"
    # URL to fetch the sequence

    try:
        print(f"Fetching protein sequence for ID {protein_id} in {format} "
              f"format...")
        response = requests.get(base_url)  # Sends GET request to download
        #  the sequence
        if response.status_code == 200:
            with open(output_file, "w") as file: # Opens the output file in
                # write mode and saves the sequence data
                file.write(response.text)
            print(f"Protein sequence for ID {protein_id} saved to "
                  f"{output_file}")
        else: # Prints error message if request fails
            print(f"Failed to download protein sequence for ID {protein_id}.")
            print("Response:", response.text)
    except Exception as e: # Handles unexpected errors
        print(f"An error occurred while downloading the sequence: {e}")

# Function to download sequences for multiple species
def download_protein_for_species(protein_name, species_list, output_folder,
                                 format="fasta"):
    """
    Loops through list of species, fetches their protein sequences,
    and saves them to individual files.

    Parameters:
        protein_name (str): Name of the protein ("Heme oxygenase 1").
        species_list (list): List of species names to query.
        output_folder (str): Folder to save the downloaded sequences.
        format (str): Format for the sequences ("fasta").
    """
    if not os.path.exists(output_folder):  # Checks if the output folder exists
        os.makedirs(output_folder)  # Creates the folder if it does not exist

    for species in species_list:  # Loops through each species in the list
        protein_id = fetch_protein_id(protein_name, species) # Fetches the
        # UniProt ID for the protein in the current species (Calls
        # fetch_protein_id function)
        if protein_id:  # If a valid UniProt ID is found
            output_file = os.path.join(output_folder,
            f"{species.replace(' ', '_')}_{protein_id}.{format}")  # Creates
            # the output file path and replaces spaces in the species name with
            # underscores
            download_protein_sequence(protein_id, output_file, format=format)
            # Downloads and saves the protein sequence (Calls
            # download_protein_sequence function)
        else:
            print(f"Could not find protein ID for {protein_name} in "
                  f"{species}.") # Prints an error message if no UniProt ID
            # is found

# Function to summarise alignment results
def summarise_alignment(file_path):
    """
    Reads an alignment file and prints a summary of the alignment.

    Parameters:
        file_path (str): Path to the alignment file.
    """
    file_format = "clustal" if "clustal" in file_path else "fasta"
    # Determines the file format based on its name
    try:
        alignment = AlignIO.read(file_path, file_format)  # Reads the
        # alignment using AlignIO
        print(f"\nSummary of alignment from {file_path}:")
        print(f"  Number of sequences: {len(alignment)}")  # Prints the number
        # of sequences
        print(f"  Alignment length: {alignment.get_alignment_length()}")
        # Prints the alignment length
        alignment_score = calculate_alignment_score(alignment)  # Calculates
        # and prints the alignment score (Calls calculate_alignment_score
        # function)
        print(f"  Alignment score (percent identical positions): "
              f" {alignment_score:.2f}%")
        print(f"  Sequences:")
        for record in alignment:
            print(f"    {record.id}: {record.seq[:50]}"
                  f"{'...' if len(record.seq) > 50 else ''}") # Prints the
            # first 50 characters of each sequence
    except Exception as e: # Handles errors in reading the alignment file
        print(f"Failed to parse alignment file {file_path}. Error: {e}")

# Function to calculate the alignment score
def calculate_alignment_score(alignment):
    """
    Calculates a basic alignment score as the percentage of identical
    positions across all sequences.

    Parameters:
        alignment: The alignment object created by AlignIO.

    Returns:
        float: The percentage of identical positions.
    """
    alignment_length = alignment.get_alignment_length()  # Gets the length of
    # the alignment
    identical_positions = 0  # Initiates counter for identical positions

    for i in range(alignment_length): # Loops through each position in the
        # alignment
        residues_at_position = [record.seq[i] for record in alignment]
        # Collects all residues at this position
        if (len(set(residues_at_position)) == 1 and '-' not in
                residues_at_position): # Checks if all residues are
            # identical and there are no gaps
            identical_positions += 1 # Increments the counter for identical
            # positions

    return (identical_positions / alignment_length) * 100  # Calculates the
    # percentage of identical positions

# Main section of the script
if __name__ == "__main__":
    # Step 1: Defining inputs for the script
    protein_name = "Heme oxygenase 1"  # Name of protein
    species_list = ["Bos taurus", "Homo sapiens",
                    "Sus scrofa", "Rattus norvegicus", "Gallus gallus",
                    "Delphinapterus leucas"]  # List of species
    output_folder = "protein_sequences"  # Folder to save downloaded sequences
    sequence_format = "fasta"  # Sequence format

    # Step 2: Downloading protein sequences
    download_protein_for_species(protein_name, species_list, output_folder,
                                 format=sequence_format) # Calls
    # download_protein_for_species function

    # Step 3: Combining all downloaded sequences into a single file
    msa_output_folder = "msa_protein_sequences"  # Folder for alignment outputs
    os.makedirs(msa_output_folder, exist_ok=True)  # Creates folder if it
    # does not exist
    combined_fasta = os.path.join(msa_output_folder, "combined_sequences.fasta")
    # Path for combined sequences

    with open(combined_fasta, 'w') as outfile:  # Opens combined file in write
        #  mode
        for fasta_file in os.listdir(output_folder):  # Loops through each
            # file  in the output folder
            if fasta_file.endswith('.fasta'):  # Processes only FASTA files
                with open(os.path.join(output_folder, fasta_file)) as infile:
                # Opens each FASTA file
                    outfile.write(infile.read())  # Appends its content to the
                #  combined file
    print(f"Combined all FASTA files into {combined_fasta}.")

    # Step 4: Performing multiple sequence alignment using Clustal Omega
    clustal_output = os.path.join(msa_output_folder, "combined_clustal.aln")
    # Path for Clustal output
    subprocess.run([
        "clustalo",  # Command to call Clustal Omega
        "-i", combined_fasta,  # Input file
        "-o", clustal_output,  # Output file
        "--outfmt", "clu",  # Output format
        "--force"  # Overwrites existing file if necessary
    ])

    # Step 5: Performing multiple sequence alignment using MAFFT
    mafft_output = os.path.join(msa_output_folder, "combined_mafft.aln")
    # Path for MAFFT output
    with open(mafft_output, 'w') as mafft_out:  # Opens output file in write
        # mode
        subprocess.run(["mafft", combined_fasta], stdout=mafft_out)
        # Runs MAFFT and redirects output to file

    # Step 6: Summarising alignment results
    summarise_alignment(clustal_output)  # Summarises Clustal Omega results
    # (Calls summarise_alignment function)
    summarise_alignment(mafft_output)  # Summarises MAFFT results (Calls
    # summarise_alignment function)