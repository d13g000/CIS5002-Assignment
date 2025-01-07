import os
import subprocess
import requests
from Bio import SeqIO

# Define the species and UniProt accession numbers (replace these IDs with actual ones for Heme oxygenase 1)
species = {
    "Bos taurus": "Q5E9F2",
    "Homo sapiens": "P09601",
    "Sus scrofa": "P32394",
    "Rattus norvegicus": "P06762",
    "Gallus gallus": "P14791",
    "Delphinapterus leucas": "A0A2Y9QJW1"
}

# Directory for storing downloaded sequences
os.makedirs("sequences", exist_ok=True)

# Download sequences from UniProt
def download_sequences():
    print("Downloading sequences...")
    for species_name, accession in species.items():
        url = f"https://www.uniprot.org/uniprot/{accession}.fasta"
        response = requests.get(url)
        if response.status_code == 200:
            with open(f"sequences/{species_name.replace(' ', '_')}.fasta", "w") as file:
                file.write(response.text)
            print(f"Downloaded: {species_name}")
        else:
            print(f"Failed to download sequence for {species_name}.")

def run_alignment(tool):
    # Combine all sequences into a single file
    combined_file = "sequences/combined.fasta"
    with open(combined_file, "w") as outfile:
        for species_name in species.keys():
            infile = f"sequences/{species_name.replace(' ', '_')}.fasta"
            for record in SeqIO.parse(infile, "fasta"):
                SeqIO.write(record, outfile, "fasta")

    # Run alignment
    aligned_file = f"aligned_{tool}.fasta"

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
