# CIS5002-Assignment
CIS5002 Assignment
To run this project, the entire folder ("202425_CIS5002_AGIUS_DIEGO_168499M_assignment_code") must be downloaded from the University of Malta Google Drive and executed as is. 
Each script within the folder relies on the existing directory structure, including specific file names and locations, to function correctly.
Modifying the folder structure or renaming files may result in errors, as the scripts reference predefined paths using the OS library. 
If you wish to use your own files, ensure that the new files are placed within the same directory structure and named exactly as referenced in the scripts. Alternatively, you may modify the scripts accordingly to accommodate different file names or locations. 


---
1a) Download species sequences and MSA
This Python script automates the process of downloading protein sequences from UniProt, preparing them for analysis, and performing multiple sequence alignments using Clustal Omega ("clustalo") and MAFFT. 
The script is designed to handle sequences from multiple species and outputs the aligned sequences into specified files ("1a_{tool}_alignment.fasta").

--> Features:
 - Downloads protein sequences for specified species using their scientific species name and respective Heme oxygenase 1 UniProt accession number. 
 - Combines downloaded sequences into a single file for alignment.
 - Performs multiple sequence alignment using Clustal Omega and MAFFT and saves the aligned sequences in separate output files.

--> Prerequisites:
 - Requests (Library to facilitate HTTP requests and fetching data)
   - Linux download: pip install requests / conda install requests
   - iOS download: brew install requests / conda install requests
 - Clustal Omega (Tool to run multiple sequence alignment)
   - download: http://www.clustal.org/omega/#Download
 - MAFFT (Tool to run multiple sequence alignment)
   - download: https://mafft.cbrc.jp/alignment/software/
 - The following files and directory structure:
    - Task 1 folder/directory 
      - 1a_species_sequences (folder/directory containing downloaded sequences)
 
--> Input files:
 - NONE (script uses defined species and corresponding accession number dictionary)

--> Output files:
 - Within Task 1 directory:
   - 1a_{tool}_alignment.fasta (FASTA files corresponding to Clustal Omega and MAFFT alignment results)
   - The following files and directory structure:
     - Species sequences (folder/directory containing downloaded sequences)
       - 1a_{species_name}.fasta (FASTA files containing sequence from specific species)
 

---
1b) MariaDB upload
This Python script automates the process of creating a MariaDB database (if it does not already exist) and importing an SQL dump file ("go-dump.sql") into the database.

--> Features:
 - Verifies the existence of the specified database (gene_ontology) in MariaDB and creates it if absent.
 - Imports SQL dump files into the database, parsing and executing the SQL statements individually for flexibility and error handling.

--> Prerequisites:
 - mysql-connector (Library to connect and interact with MariaDB)
   - Linux download: pip install mysql-connector-python / conda install mysql-connector-python
   - iOS download: brew install mysql-connector-python / conda install mysql-connector-python
 - sqlparse (Library to parse and split SQL scripts)
   - Linux download: pip install sqlparse / conda install sqlparse
   - iOS download: brew install sqlparse / conda install sqlparse
 - MariaDB (Server to store SQL files)
   - Linux download:
       # Install MariaDB
       - sudo apt update 
       - sudo apt install mariadb-server mariadb-client
       # Start and Enable MariaDB
       - sudo systemctl start mariadb 
       - sudo systemctl enable mariadb
       # Set root password and configure security settings
       - sudo mysql_secure_installation
       # Log in and test installation
       - mysql -u root -p
   - iOS download:
       # Install MariaDB
       - brew install mariadb
       # Start MariaDB server
       - brew services start mariadb
       # Test installation
       - mysql -u root
 - The following files and directory structure:
    - Task 1 folder/directory 
      - go-dump.sql (SQL file to import into the MariaDB database)
 
--> Input files:
 - Gene ontology file ("go-dump.sql")

--> Output files:
 - NONE (script creates a MariaDB database named gene_ontology including the data defined in the go-dump.sql file)


---
2b) HFE gene extraction
This Python script extracts the genomic sequence of the HFE gene from a reference genome ("Reference genome.fna") based on a GTF annotation file ("Annotations.gtf").
It identifies the gene's location using its GeneID and Transcript ID, retrieves the sequence from the reference genome, and writes it to a FASTA file ("2b_HFE_gene.fasta").

--> Features:
 - Extracts the genomic information and location of the HFE gene using its GeneID and a specific Transcript ID.
 - Retrieves the gene sequence (accounting for strand direction) from the reference genome using the extracted genomic information and location.
 - Outputs the sequence in a FASTA file for further analysis.

--> Prerequisites:
 - The following files and directory structure:
    - Task 2 folder/directory 
      - Reference genome.fna file (Reference genome downloaded from NCBI including sequences for all chromosomes)
      - Annotations.gtf file (GTF annotations file downloaded from NCBI containing genomic annotations, including gene and transcript information)

--> Input files:
 - Reference genome ("Reference genome.fna") https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
 - Annotations ("Annotations.gtf") https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/

--> Output files:
 - Within Task 2 directory:
   - 2b_HFE_gene.fasta (FASTA file containing extracted HFE gene sequence)
     - uppercase letters: regions where nucleotides are highly conserved
     - lowercase letters: regions where nucleotides vary with most common nucleotide given in lowercase (soft-masked region) 


---
2c) HFE transcription and mRNA variants
This Python script extracts and analyses mRNA variants for the HFE gene. 
Using reference genome ("Reference genome.fna") and GTF annotations ("Annotations.gtf") files, the script identifies transcript variants, extracts exon sequences, transcribes DNA to mRNA, and outputs a folder ("2c_HFE_mRNA_variants") containing variant mRNA sequences in FASTA format ("2c_{transcript_id}_mRNA.fasta").
Using the DNA sequence ("2b_HFE_gene.fasta") file, the script transcribes DNA to mRNA and outputs the mRNA sequence in FASTA format ("2c_HFE_gene_mRNA.fasta").

--> Features:
 - Identifies transcript variants of the HFE gene using its GeneID.
 - Extracts exon sequences (accounting for strand direction) based on transcript information and combines them into a full gene sequence (DNA).
 - Transcribes DNA to mRNA by replacing thymine (T/t) with uracil (U/u).
 - Outputs mRNA sequences for each transcript variant in individual FASTA files for further analysis.

--> Prerequisites:
 - The following files and directory structure:
    - Task 2 folder/directory 
      - Reference genome.fna file (Reference genome downloaded from NCBI including sequences for all chromosomes)
      - Annotations.gtf file (GTF annotations file downloaded from NCBI containing genomic annotations, including gene and transcript information)
      - 2b_HFE_gene.fasta file (FASTA file containing HFE gene sequence extracted in Task 2b)
 
--> Input files:
 - Reference genome ("Reference genome.fna") https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
 - Annotations ("Annotations.gtf") https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
 - HFE gene sequence ("2b_HFE_gene.fasta")

--> Output files:
 - Within Task 2 directory:
   - 2c_HFE_gene_mRNA.fasta (FASTA file containing transcribed HFE mRNA sequence)
     - uppercase letters: regions where nucleotides are highly conserved
     - lowercase letters: regions where nucleotides vary with most common nucleotide given in lowercase (soft-masked region)
   - The following files and directory structure:
     - 2c_HFE_mRNA_variants (folder/directory containing extracted HFE mRNA variants)
       - 2c_{HFE transcript ID}_mRNA.fasta (FASTA file containing transcribed HFE variants mRNA sequence)


---
2d) HFE translation and ORF
This Python script translates the mRNA sequence of the HFE gene from Task 2c into protein sequences, taking into consideration all six reading frames, generating a folder ("2d_HFE_protein_reading_frames") of FASTA files ("2d_HFE_protein_reading_frame_{frame}.fasta").
It also determines the open reading frame (ORF) and saves it to an individual FASTA file ("2d_HFE_ORF_sequence.fasta").

--> Features:
 - Translates the given mRNA sequence into protein sequences arising from the three forward and three reverse frames.
 - Identifies the ORF.
 - Saves the protein sequences for each frame and the ORF sequence into separate FASTA files for further analysis.

--> Prerequisites:
 - The following files and directory structure:
    - Task 2 folder/directory 
      - 2c_HFE_gene_mrna.fasta file (FASTA file containing HFE mRNA sequence produced in Task 2c)

--> Input files:
 - HFE mRNA sequence ("2c_HFE_gene_mrna.fasta")

--> Output files:
 - Within Task 2 directory:
   - The following files and directory structure:
     - 2d_HFE_protein_reading_frames.fasta (folder/directory containing HFE reading frames and ORF)
       - 2d_HFE_ORF_sequence.fasta (FASTA file containing HFE ORF sequence)
       - 2d_HFE_protein_reading_frame_{frame}.fasta (FASTA file containing translated HFE reading frame)
         - letters: amino acids
         - *: stop codons


---
2e) HFE ORF/RefSeq comparison
This Python script compares the HFE ORF protein sequence extracted in Task 2d ("2d_HFE_ORF_sequence.fasta") with a reference protein sequence extracted from a RefSeq file ("Proteins.faa"). 
The script also identifies any differences between the two sequences and outputs detailed results to a text file ("2e_HFE_protein_comparison.txt").

--> Features:
 - Extracts a specific protein sequence from a RefSeq protein FASTA file based on its protein reference ID.
 - Compares two protein sequences, identifying differences at each position.
 - Handles cases where sequences are of different lengths.
 - Outputs a detailed comparison report to a text file.

--> Prerequisites:
 - The following files and directory structure:
    - Task 2 folder/directory 
      - 2d_HFE_ORF_sequence.fasta file (FASTA file containing HFE ORF sequence)
      - Proteins.faa file (Translated proteins downloaded from NCBI)

--> Input files:
 - Proteins ("Proteins.faa") https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
   - Protein reference ID of protein sequence in RefSeq file ("NP_000401.1")
 - HFE ORF sequence ("2d_HFE_ORF_sequence.fasta")

--> Output files:
 - Within Task 2 directory:
   - 2e_HFE_protein_comparison.txt (text file containing the comparison results, including position-wise differences and extra residues if either of the sequences is longer)


---
2f) HFE pathogenic protein variants
This script extracts and analyses pathogenic variants specific to the HFE gene from a ClinVar file ("Clinvar.txt"). 
It identifies variants classified as pathogenic, formats the information, and saves the output to a CSV file ("HFE_pathogenic_variants.csv") for further analysis. 
The script also provides detailed insights into the extracted variants by printing them to the console.

--> Features:
 - Filters variants for the HFE gene classified as Pathogenic.
 - Extracts key details, including:
   - Variation (variant)
   - Gene name and corresponding protein changes (if available)
   - Mutation type
   - Associated condition/s
   - Classification (pathogenic)
 - Saves the extracted details in a CSV file. 
 - Displays a summary of the results in the console.

--> Prerequisites:
 - The following files and directory structure:
    - Task 2 folder/directory 
      - Clinvar.txt file (text file containing ClinVar information)

--> Input files:
 - ClinVar ("Clinvar.txt") https://www.ncbi.nlm.nih.gov/clinvar/?term=%22hfe%22%5BGENE%5D&redir=gene

--> Output files:
 - Within Task 2 directory:
   - 2f_HFE_pathogenic_variants.csv (CSV file containing the columns for variant, gene (protein change), type (consequence), condition and, classification)