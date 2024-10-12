# BI_Tools

## Repository description

**BI_Tools** is a set of utilities for working with nucleotide sequences (DNA and RNA), filtering fastq-sequences, modifying fasta-files, parsing BLAST output and selecting genes form gbk-files

**Author**: [*Kulenich Viktoriia*](https://github.com/vickulenich)

## Content

### [Installation and launch](#installation-and-launch)
### [Input data](#input-data)
### [Available operations](#available-operations)

## Installation and launch

For installation you need to clone the repository:

    git clone git@github.com:vickulenich/BI_Tools.git

Then you need to go to the root directory of this repository and run the script you need:

    cd BI_Tools
    python bi_tools.py
    python bio_files_processor.py

Enter the required function and its arguments according to the input requirements

## Input data

The program consists of two main scripts. 

bi_tools.py includes two main functions:

1. *run_dna_rna_tools* takes as input any number of nucleotide sequence strings and the last argument with the name of one of the [available operations](#available-operations) and performs a user-selected transformation.

Input data example:

> run_dna_rna_tools("ATGca", "AgTCG", "transcribe")
>
> run_dna_rna_tools("AgTCG", "AcgTcAG", "reverse")
>
> run_dna_rna_tools("AUaG", "CUacG", "AcgTcAG", "complement")

2. *filter_fastq* takes as input:
   1) *input_fastq* - a path to the input fastq file (.fastq)
   2) *output_fastq* - a path to the output fastq file (.fastq)
   3) *length_bounds* - a tuple containing two integers corresponding to the lower and upper bounds of the required fastq-sequence length, by default takes the value (0, 2^32)
   4) *gc_bounds* - a tuple containing two integers corresponding to the lower and upper bounds of the required GC-composition of the fastq-sequence, by default takes the value (0, 100)
   5) *quality_threshold* - an integer corresponding to the lower basecall quality threshold for the fastq-sequence, by default takes the value 0

Input data example:

> filter_fastq(input_fastq = "example_data.fastq", length_bounds=(0,15), gc_bounds=(5,10), quality_threshold=20)

bio_files_processor.py includes three functions:

1. *convert_multiline_fasta_to_oneline* takes as input a path to the input fasta-file (.fasta/.fa) and optionally a a path to the output fasta-file, reads input file where the sequence (DNA/RNA/protein etc.) can be split into several lines and saves it to a new fasta file where each sequence fits on one line

Input data example:

> convert_multiline_fasta_to_oneline(input_fasta = "example_data.fa", output_fasta = "example_data_output.fa")
>
> convert_multiline_fasta_to_oneline(input_fasta = "example_data.fa")

2. *parse_blast_output* takes as input a path to the input txt-file (.txt) and a path to the output txt-file, reads input file, for each QUERY request selects the first row from the Description column and saves the set of obtained proteins in the output file in one column sorted alphabetically

Input data example:

> parse_blast_output(input_file = "example_data.txt", output_file = "example_data.txt")

3. *select_genes_from_gbk_to_fasta* takes as input:
    1) *input_gbk* - a path to the input gbk-file (.gbk)
    2) *genes* - list of the selected genes names
    3) *n_before* - an integer corresponding to the number of genes before each of the genes of interest, by default takes the value 1
    4) *n_after* - an integer corresponding to the number of genes after each of the genes of interest, by default takes the value 1
    5) *output_fasta* - name of the output file

Input data example:

> select_genes_from_gbk_to_fasta(input_gbk = "example_data.gbk", genes = ["gene_1", "gene_2", "gene_3"], output_fasta = "example_data.fasta", n_before = 1, n_after = 1)

## Available operations

The run_dna_rna_tools function allows you to perform the following operations on a nucleotide sequence:

- *reverse* - returns the reverse sequence of the original
- *complement* - returns the complement sequence of the original
- *reverse_complement* - returns the reverse and complement sequence of the original
- *transcribe* - returns the corresponding RNA sequence if the original sequence was DNA, or returns the original sequence if the original sequence was RNA
- *g_c_bound* - calculates the GC content of a sequence in %, rounded to 2 decimal places

The filter_fastq function allows you to select fastq sequences that meet specified requirements for their length, GC composition, and quality level.

The convert_multiline_fasta_to_oneline function allows you to convert the file where the sequence can be split into several lines to a new file where each sequence fits on one line.

The parse_blast_output function allows you to extract the name of the best BLAST match from the database and save all the results in one file - **Unavaliable**

The select_genes_from_gbk_to_fasta function allows you to select a certain number of genes before and after each gene of interest and save their protein sequence to a fasta file that can be sent to the BLAST input - **Unavaliable**