import requests # Library used to interact with UniProt and fetch protein IDs
# and sequences
import os # Library used to interact with the operating system


def fetch_protein_id(protein_name, species_name):
    """
    Fetches UniProt ID for a given protein and species.

    Parameters:
        protein_name (str): Name of the protein ("Heme oxygenase 1).
        species_name (str): Name of the species.

    Returns:
        str: UniProt ID of the protein for the species, or error if not found.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search" # Base URL for
    # UniProt website
    query = (f"query={protein_name} AND (organism_name:"
             f"{species_name})&fields=accession") # Adds desired protein and
    # species names (requests UniProt ID/"accession" field)
    url = f"{base_url}?{query}" # Combines base URL to desired protein and
    # species names to generate whole URL (opens UniProt website and
    # navigates to desired protein and species)

    try:
        response = requests.get(url) # Sends get request to UniProt API
        # Handling API response
        if response.status_code == 200: # Checks if API request is successful
            data = response.json() # Converts API response (JSON format) into
            # dictionary
            if "results" in data and data["results"]: # Retrieves results
                # from response
                return data["results"][0]["primaryAccession"] # Returns
                # UniProt ID (if results exist)
        # Error handling
        else:
            print(f"Failed to fetch protein ID for {protein_name} in "
                  f"{species_name}.")
            print("Response:", response.text) # Prints error and server
            # response - UniProt website errors
    except Exception as e: # Prevents programme crashing and terminating due to
        # error - Python error (network/library/JSON parsing)
        print(f"An error occurred while fetching protein ID: {e}")
    return None # Allows programme to continue processing other species
    # within the list without crashing entirely


def download_protein_sequence(protein_id, output_file, format="fasta"):
    """
    Downloads specific protein sequence in FASTA format using UniProt ID.

    Parameters:
        protein_id (str): UniProt ID of the protein.
        output_file (str): Path to save the protein sequence (downloaded file).
        format (str): Desired format for the sequence ("FASTA").

    Returns:
        Downloads protein sequence in FASTA format, or error if not found.
    """
    base_url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.{format}"
    # URL to fetch protein sequence in desired format

    try:
        print(
            f"Fetching protein sequence for ID {protein_id} in {format} format...")
        response = requests.get(base_url)
        if response.status_code == 200:
            with open(output_file, "w") as file: # Opens specified file in
                # 'write' mode and save protein sequence within it
                file.write(response.text)
            print(
                f"Protein sequence for ID {protein_id} saved to {output_file}")
        else:
            print(f"Failed to download protein sequence for ID {protein_id}.")
            print("Response:", response.text)
    except Exception as e:
        print(f"An error occurred while downloading the sequence: {e}")


def download_protein_for_species(protein_name, species_list, output_folder,
                                 format="fasta"):
    """
    Downloads FASTA protein sequences for specified list of species.

    Parameters:
        protein_name (str): Name of the protein.
        species_list (list): List of species names.
        output_folder (str): Folder to save the sequences within.
        format (str): Desired format for the sequences ("FASTA").

    Returns:
        Creates folder and adds downloaded protein sequences into it, or error
        if not found.
    """
    if not os.path.exists(output_folder): # Checks if folder exists
        os.makedirs(output_folder) # Creates folder

    for species in species_list: # Iterates through each species in list
        protein_id = fetch_protein_id(protein_name, species) # Fetches
        # UniProt ID for protein in given sequence (calls fetch_protein_id
        # function)
        if protein_id:
            output_file = os.path.join(output_folder,
                                       f"{species.replace(' ', '_')}_"
                                       f"{protein_id}.{format}") # Creates
            # full path for output file and names it
            download_protein_sequence(protein_id, output_file, format=format)
            # Downloads sequence for given UniProt ID in FASTA format (calls
            # download_protein_sequence function)
        else:
            print(f"Could not find protein ID for {protein_name} in {species}.")


if __name__ == "__main__":
    protein_name = "Heme oxygenase 1" # Target protein
    species_list = ["Bos taurus", "Homo sapiens", "Sus scrofa",
                    "Rattus norvegicus","Gallus gallus",
                    "Delphinapterus leucas"] # Species list
    output_folder = "protein_sequences" # Folder within which protein
    # sequence files are saved
    sequence_format = "fasta" # File format

    download_protein_for_species(protein_name, species_list, output_folder,
                                 format=sequence_format) # Fetches and saves
    # protein sequences (calls main function)