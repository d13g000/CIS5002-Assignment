import os # Import os library to interact with operating system
import subprocess # Import subprocess to allow running external commands and
# processes
import requests # Import requests to facilitate HTTP requests to fetch data
from Bio import SeqIO # Import SeqIO part from Biopython library to be able
# to parse and handle biological sequence data

# Define species and UniProt accession numbers (dictionary of species and
# associated Heme oxygenase 1 sequence accession number)
species = {
    "Bos taurus": "Q5E9F2", # Cow
    "Homo sapiens": "P09601", # Human
    "Sus scrofa": "P32394", # Pig
    "Rattus norvegicus": "P06762", # Rat
    "Gallus gallus": "P14791", # Chicken
    "Delphinapterus leucas": "A0A2Y9QJW1" # Beluga whale
}

# Create main directory to store all files/folders downloaded/generated
main_dir = "Task 1"
os.makedirs(main_dir, exist_ok=True) # Create directory if nonexistent

# Create subdirectory to storing downloaded files for sequences from
# different species
sequences_dir = os.path.join(main_dir, "Species sequences")
os.makedirs(sequences_dir, exist_ok=True) # Create subdirectory if nonexistent


# Download sequences from UniProt using species dictionary
def download_sequences():
    print("Downloading sequences...") # Notify start of process
    for species_name, accession in species.items(): # Iterate through each
        # species and associated accession number
        url = f"https://www.uniprot.org/uniprot/{accession}.fasta"
        # Construct UniProt URL for sequence
        response = requests.get(url) # Fetch sequence data from URL
        if response.status_code == 200: # Check if request is successful
            file_path = os.path.join (sequences_dir,
                        f"{species_name.replace(' ', '_')}.fasta")
            # Set file path and name within which to save sequence
            with open(file_path, "w") as file: # Open file for writing "w"
                file.write(response.text) # Write and save fetched sequence
                # into file
            print(f"Successfully downloaded: {species_name}") # Print
            # notification of success
        else:
            print(f"Failed to download sequence for {species_name}.") # Print
            # notification of failure


# Add empty line before each header in FASTA file
def add_empty_line(file_path):
    with open(file_path, "r") as infile: # Open file for reading "r"
        lines = infile.readlines() # Read all lines
    with open(file_path, "w") as outfile:
        for line in lines:
            if line.startswith(">"): # Determine if line is a header
                outfile.write("\n") # Add empty line before header
            outfile.write(line) # Write line into file


# Perform sequence alignment using specified tool (Clustal Omega / MAFFT)
def run_alignment(tool):
    # Combine all individual sequence files into a single file for multiple
    # alignment
    combined_file = os.path.join(sequences_dir, "Combined_sequences.fasta")
    with open(combined_file, "w") as outfile:
        for species_name in species.keys(): # Iterate through species
            infile = os.path.join(sequences_dir,
                        f"{species_name.replace(' ', '_')}.fasta")
            for record in SeqIO.parse(infile, "fasta"): # Parse FASTA file
                # using Biopython
                SeqIO.write(record, outfile, "fasta") # Write each record
                # into combined file
    add_empty_line(combined_file) # Add empty line before each header (Call
    # add_empty_line function)

    # Define output file and run alignment
    aligned_file = os.path.join(main_dir, f"{tool}_alignment.fasta")
    # Clustal Omega
    if tool == "clustalo":
        print("Running Clustal Omega alignment...") # Notify of alignment
        # tool being used
        subprocess.run([
            "clustalo", "-i", combined_file, "-o", aligned_file, "--force"
        ]) # Run Clustal Omega with "combined_file" as input and "aligned_file"
        # as output
    # MAFFT
    elif tool == "mafft":
        print("Running MAFFT alignment...")
        subprocess.run([
            "mafft", combined_file], stdout=open(aligned_file, "w"))

    add_empty_line(aligned_file)
    print(f"Alignment successfully completed using {tool}. Output "
          f"successfully saved to {aligned_file}") # Print notification of
    # successful alignment

def main():
    # Step 1: Download species sequences from UniProt
    download_sequences() # Call download_sequences function

    # Step 2: Run both alignments
    run_alignment("clustalo") # Call run_alignment function using Clustal Omega
    run_alignment("mafft") # Call run_alignment function using MAFFT

if __name__ == "__main__":
    main() # Call main function