import os
import subprocess
import requests
from Bio import SeqIO

# Define the species and UniProt accession numbers (dictionary of
# species name and associated sequence accession number)
species = {
    "Bos taurus": "Q5E9F2",
    "Homo sapiens": "P09601",
    "Sus scrofa": "P32394",
    "Rattus norvegicus": "P06762",
    "Gallus gallus": "P14791",
    "Delphinapterus leucas": "A0A2Y9QJW1"
}

# Create main directory for all files (all files/folders downloaded/generated
# will be saved within this folder)
main_dir = "Task 1"
os.makedirs(main_dir, exist_ok=True)

# Subdirectory for storing downloaded sequences (files for sequences from
# different species will be saved within this folder)
sequences_dir = os.path.join(main_dir, "species sequences")
os.makedirs(sequences_dir, exist_ok=True)

# Download sequences from UniProt (uses species dictionary to download
# sequences)
def download_sequences():
    print("Downloading sequences...")
    for species_name, accession in species.items():
        url = f"https://www.uniprot.org/uniprot/{accession}.fasta"
        response = requests.get(url)
        if response.status_code == 200:
            file_path = os.path.join (sequences_dir,
                        f"{species_name.replace(' ', '_')}.fasta")
            with open(file_path, "w") as file:
                file.write(response.text)
            print(f"Downloaded: {species_name}")
        else:
            print(f"Failed to download sequence for {species_name}.")

def run_alignment(tool):
    # Combine all sequences into a single file (required for multiple alignment)
    combined_file = os.path.join(sequences_dir, "combined.fasta")
    with open(combined_file, "w") as outfile:
        for species_name in species.keys():
            infile = os.path.join(sequences_dir,
                        f"{species_name.replace(' ', '_')}.fasta")
            for record in SeqIO.parse(infile, "fasta"):
                SeqIO.write(record, outfile, "fasta")

    # Run alignment and save outputs in the main directory
    aligned_file = os.path.join(main_dir, f"aligned_{tool}.fasta")

    if tool == "clustalo":
        print("Running Clustal Omega alignment...")
        subprocess.run([
            "clustalo", "-i", combined_file, "-o", aligned_file, "--force"
        ])
    elif tool == "mafft":
        print("Running MAFFT alignment...")
        subprocess.run([
            "mafft", combined_file], stdout=open(aligned_file, "w"))

    print(f"Alignment completed using {tool}. Output: {aligned_file}")

def main():
    download_sequences()

    # Run both alignments
    run_alignment("clustalo")
    run_alignment("mafft")

if __name__ == "__main__":
    main()
