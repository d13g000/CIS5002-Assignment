import os

# Define directory and file paths
script_dir = os.path.dirname(os.path.abspath(__file__))
task_2_dir = os.path.join(script_dir, "Task 2")

clinvar_file = os.path.join(task_2_dir, "Clinvar.txt")  # Path to the ClinVar
# file
output_file = os.path.join(task_2_dir, "2f_HFE_pathogenic_variants.csv")  # Path
# and name of output file


def extract_variant_details(columns):
    """
    Extracts details pertaining to each variant.

    Args:
        columns (list): List of column values from ClinVar input file.

    Returns:
        tuple: Extracted pathogenic details (variation, gene (protein change),
        mutation type (consequence), condition, classification).
    """
    # Extract pathological variant information
    variation = columns[0] if len(columns) > 0 else "Not Available"
    # Extract variant variation or default to N/A if missing (column 1)
    genes = columns[1] if len(columns) > 1 else "Not Available"  # Extract
    # variant gene (protein change) or default to N/A if missing (column 2)
    mutation_type = columns[13].strip().capitalize() if len(columns) > 13 \
        else "Not Available"  # Extract variant mutation type (consequence)
    # or default to N/A if missing (column 14)
    condition = columns[3] if len(columns) > 3 else "Not Available"
    # Extract variant condition or default to N/A if missing (column 4)
    classification = columns[15].strip() if len(columns) > 15 \
        else "Not Available"  # Extract variant classification or default to
    # N/A if missing (column 16)

    # Extract gene and protein change details
    genes_with_protein_changes = [] # Empty list to store gene and
    # corresponding protein changes
    if "p." in variation:  # Check if protein changes exist in variation string
        protein_changes = [f"{gene} (p.{p_change.split()[0]})"
                           for gene, p_change in
                           zip(genes.split(","), variation.split("p.")[1:])]
        # Extract protein changes and genes effected from variant variation
        # column
        genes_with_protein_changes.extend(protein_changes) # Add formatted
        # changes to list
    else:
        genes_with_protein_changes = genes.split(",")  # Use genes effected
        # only if no protein changes are determined

    return (variation, ", ".join(genes_with_protein_changes), mutation_type,
            condition, classification) # Return extracted details


def main():
    try:
        # Step 1: a) Open and read all lines from "Clinvar.txt" file
        with open(clinvar_file, 'r') as file:
            lines = file.readlines()


        hfe_variants = [] # List to hold the pathogenic variants for the HFE
        # gene

        # b) Parse file line by line
        for line in lines:
            columns = line.strip().split("\t")  # Split line into columns by
            # tab
            if len(columns) > 15 and 'HFE' in columns[1] and 'Pathogenic' in columns[15]:
                hfe_variants.append(columns)  # Append pathogenic variants
                # for HFE only to the list

        # Count number of pathogenic variants
        num_pathogenic_variants = len(hfe_variants)

        # c) Write output to file
        with (open(output_file, 'w') as output):
            output.write("Variation\tGene (Protein Change)\tType (Consequence)\t"
                         "Condition\tClassification\n") # Write header
            # representing column titles

            if num_pathogenic_variants > 0:
                print(f"Successfully found {num_pathogenic_variants} "
                      f"pathogenic variants for HFE") # Print number of
                # pathogenic variants determined
                print("__" * 120) # Print separator line

                for variant_columns in hfe_variants:
                    variation, genes_protein, mutation_type, condition,\
                    classification = extract_variant_details(variant_columns)
                    # Process each variant and extract its details (Call
                    # extract_variant_details function)

                    # Print details to console
                    print(f"Variant - {variation}") # Variant
                    print(f"Gene (Protein change) - {genes_protein}") # Genes
                    # affected and protein changes
                    print(f"Type (Consequence) - {mutation_type}") # Type of
                    # mutation
                    print(f"Condition - {condition}") # Condition
                    print(f"Classification - {classification}")
                    # Classification
                    print("__" * 120)

                    # Write details to output file
                    output.write(f"{variation}\t{genes_protein}\t"
                        f"{mutation_type}\t{condition}\t{classification}\n")
            else:
                print("No pathogenic variants for HFE were found in the "
                      "file.")

        print(f"Pathogenic variant details successfully saved to:"
              f" {output_file}")

    except Exception as e:
        print(f"An error occurred: {e}")  # Handle any other errors


if __name__ == "__main__":
    main()