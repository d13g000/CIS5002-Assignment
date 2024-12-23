# CIS5002-Assignment
CIS5002 Assignment question explanation

Task 1a (FASTA download) handles fetching and downloading Heme oxygenase 1 protein sequences for a list of different species
1. For each species the script queries UniProt for the specific protein ID
2. If a matching UniProt ID is found, the corresponding protein ID is downloaded in FASTA format
3. The sequence is saved to a file named after the species and its corresponding UniProt ID

Task 1a (MSA) handles running multiple sequence alignment on the donwloaded FASTA files using Clustal Omega and MAFFT 
(Assumes Clustal Omega, MAFFT and, Biopython are installed)
1. Previously downloaded FASTA files are combined to generate a new combined FASTA file
2. Combined FASTA file is aligned using Clutal Omega and MAFFT as saved in seperate files where:
   Gaps (-) are inserted into the sequence where residues cannot be aligned
   Residues from the sequence of each species are aligned below one another with:
      * indicating a fully conserved position (exact match across all sequences)
      : indicating a conserved position (similar but not identical residues across sequences)
      . indicating a weakly conserved position (residues are somewhat similar but not highly conserved)
3. Alignment summary is generated including alignment score (%age sequence similarity) for each tool 
