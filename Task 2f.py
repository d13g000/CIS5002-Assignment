import os

# Define the path to the 'task 2' folder
folder_path = os.path.join(os.getcwd(), "Task 2")
file_path = os.path.join(folder_path, "clinvar.txt")

output_file_path = os.path.join(folder_path, "HFE_pathogenic_variants.csv")

# Read the content of the file
with open(file_path, 'r') as file:
    lines = file.readlines()

# Initialize a list to hold the pathogenic variants for the HFE gene
hfe_variants = []

# Parse the file
for line in lines:
    # Ensure the line is properly split by tabs and check its contents
    columns = line.strip().split("\t")
    if len(columns) > 15 and 'HFE' in columns[1] and 'Pathogenic' in columns[15]:
        hfe_variants.append(columns)

# Function to extract mutation details from a line
def extract_variant_details(columns):
    # Extract necessary details
    variation = columns[0] if len(columns) > 0 else "Not Available"
    genes = columns[1] if len(columns) > 1 else "Not Available"
    mutation_type = columns[13].strip().capitalize() if len(columns) > 13 else "Not Available"
    condition = columns[3] if len(columns) > 3 else "Not Available"
    classification = columns[15].strip() if len(columns) > 15 else "Not Available"

    # Extract gene and protein change details
    genes_with_protein_changes = []
    if "p." in variation:
        protein_changes = [f"{gene} ({'p.' + p_change.split()[0]})"
                           for gene, p_change in
                           zip(genes.split(","), variation.split("p.")[1:])]
        genes_with_protein_changes.extend(protein_changes)
    else:
        genes_with_protein_changes = genes.split(",")

    return variation, ", ".join(genes_with_protein_changes), mutation_type, condition, classification

# Count the number of pathogenic variants
num_pathogenic_variants = len(hfe_variants)

# Open the output file for writing
with open(output_file_path, 'w') as output_file:
    # Write the header row
    output_file.write("Variation\tGene (Protein Change)\tType (Consequence)\tCondition\tClassification\n")

    # Print the count of pathogenic variants
    if num_pathogenic_variants > 0:
        print(f"Found {num_pathogenic_variants} pathogenic variants for HFE")
        print("-" * 120)
        # Write details of each variant to the output file and print to console
        for variant_columns in hfe_variants:
            variation, genes_protein, mutation_type, condition, classification = extract_variant_details(variant_columns)
            # Print to console
            print(f"Variant - {variation}")
            print(f"Gene (Protein change) - {genes_protein}")
            print(f"Type (Consequence) - {mutation_type}")
            print(f"Condition - {condition}")
            print(f"Classification - {classification}")
            print("__" * 120)

            # Write to output file
            output_file.write(f"{variation}\t{genes_protein}\t{mutation_type}\t{condition}\t{classification}\n")
    else:
        print("No pathogenic variants for HFE were found in the file.")

# Inform the user where the file is saved
print(f"Pathogenic variant details saved to: {output_file_path}")
