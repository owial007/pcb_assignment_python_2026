# Codes Description

This Directory contains two codes. `bakend_code.py` and `main.py`.

## 1) `bakend_code.py` Code:

As the name suggests, it contains the backend code that is used to analyze the `.fasta` file and the `variatnt_csv` file. Below is description of major classes and functions used in the code. 

### Functions

- `logger(start_time, end_time, update)`: Appends a timestamp, elapsed time, and update message to `log_file.txt`.

### Class `coding_seq`
This handles all the operation regarding the provided reference seq file. Following are its attributes:

- `__init__(dna_seq)`: Stores the input sequence on the object.
- `seq_check(dna_seqs=None, seq_type='Provided')`: Checks for non-ATGC characters, length divisible by 3, start codon `ATG`, and stop codon placement; raises `ValueError` on failure, otherwise returns the uppercase sequence.
- `base_cal()`: Reads the FASTA at `self.dna_seq` and returns lists of sequence lengths and GC percentages.
- `fasta_checker()`: Checks the file extension (`fasta`, `fa`, `fas`) and that it holds one record, then runs `seq_check` and returns the validated sequence.

### Class `variant_analyzer`
This handles all the operation regarding the variants. Following are its attributes:

- `__init__(csv_file)`: Stores the CSV file path.
- `varaint_csv_checker(dna_seq)`: Reads the CSV, checks required columns, and removes rows with duplicate IDs, invalid types, out-of-range positions, or invalid `change_from`/`change_to` values; returns `(removed_rows, valid_rows)`.
- `mutator(variant, dna_seq)`: Creates one mutant sequence per substitution, insertion, and deletion; writes `mutation_summary.csv` and `mutations.fas`.
- `mut_seq_checker(seq, ref_seq_len)`: Returns checks for start codon, internal stops, terminal stop, and whether the length change is in-frame, frameshift, or equal.
- `translations(dna_seq, mutations)`: Translates the reference and mutants, writes `proteins.fasta`, compares each mutant with the reference, and saves `mutation_effect.csv`.
- `protein_prop(protein_fasta)`: Calculates length, molecular weight, pI, GRAVY, and amino acid composition (Biopython and manual), with differences; saves `prot_prop.csv`.

## 2) `main.py` Code:
This file contains the main executable file which is used to perform the analysis. It imports all the attributes from the bakend_code.py and execute all the required tasks. 

It contains the following functions:
### Function:
- `parse_args()`: contains the code for importing the files from the terminal. Both the input files `{HBB_reference_CDS.fasta,variant_cds_new.csv}` files and the output dir.

### Try & Catch Block
- `first block`: Executes the validation for the reference fasta file
- `second block`: Executes the validation and tasks for the variant csv file

## Executing the Code:

To run the file use the following command after activating the conda environment (read: Environment_README.md)

`python main.py -i HBB_reference_CDS.fasta -c variant_cds_new.csv`
--or--
`python main.py -i HBB_reference_CDS.fasta -c variant_cds_new.csv -o C:\Ansh`


## Output files
The files will be stored in the provided directory as:

- `log_file.txt`: Log file of the entire execution
- `mutation_summary.csv`: Summary of all the variants used.
- `mutations.fas`: Contains the mutated sequences based on the selected variants
- `proteins.fasta`: Contains the translated sequences for all the variants and reference seq
- `mutation_effect.csv`: Contains a detailed report on the consequence of each variant
- `prot_prop.csv`: Contains a detailed report on the properties of each translated proteins 
