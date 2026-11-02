from bakend_code import variant_analyzer,logger, time,coding_seq
import argparse
import os
import shutil
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="HBB variant analysis pipeline")
    parser.add_argument("-i", "--fasta", required=True,
                        help="Reference sequence file (.fasta or .FASTA)")
    parser.add_argument("-c", "--csv", required=True,
                        help="Variant file (.csv)")
    parser.add_argument("-o", "--out", required=False,
                        help="Directory where output files are stored")
    args = parser.parse_args()

    if not args.fasta.lower().endswith(".fasta"):
        parser.error(f"-i/--fasta must be a .fasta or .FASTA file, got: {args.fasta}")
    if not args.csv.lower().endswith(".csv"):
        parser.error(f"-c/--csv must be a .csv file, got: {args.csv}")

    return args

##Checking the class on fasta files
args = parse_args()
fasta_path = args.fasta
csv = args.csv
if args.out is not None:
    out_dir = args.out
    out_dir = out_dir+"\//"
else:
    out_dir="./"

if out_dir is None:
    print("FILES WILL BE STORED IN CURRENT DIRECTORY\n\n")

#fasta_path = "HBB_reference_CDS.fasta"
try:

    start_time = time.time()
    rec=coding_seq(fasta_path)
    dna_seq=rec.fasta_checker()
    dna_len_list,gc_frac_list=rec.base_cal()
    for i in range(len(dna_len_list)):
        print("DNA Length: ",dna_len_list[i])
        print(f"GC Frac: {gc_frac_list[i]:.2f}%")
    end_time = time.time()
    logger(start_time,end_time,"Checked the Reference Seq, it passed all the parameters.")
except (ValueError) as e:
    start_time = time.time()
    print(f"Invalid Input: {e}")
    end_time = time.time()
    logger(start_time,end_time,e)



#csv = 'variant_cds_new.csv'
try:
    #checking the csv file
    start_time = time.time()
    csv=variant_analyzer(csv)
    rem_varaint,update_variant=csv.varaint_csv_checker(dna_seq)#('base.fa')
    print("extraced values")
    # display(rem_varaint)
    # display(update_variant)
    end_time = time.time()
    logger(start_time,end_time,"Checked the Varaint CSV, it passed all the parameters.")
    
    #performing mutations
    start_time=time.time()
    #print(dna_seq)
    csv.mutator(update_variant,dna_seq,out_dir)
    end_time = time.time()
    logger(start_time,end_time,"Created `mutations.fasta` file.")

    #performing translations
    start_time=time.time()
    csv.translations(dna_seq,f"{out_dir}mutations.fas",out_dir)
    end_time = time.time()
    logger(start_time,end_time,"Created `mutated_translations.fasta` file.")

    #performing protein analysis
    start_time=time.time()
    csv.protein_prop(f"{out_dir}proteins.fasta",out_dir)
    end_time = time.time()
    logger(start_time,end_time,"Created `prot_prop.fasta` file.")

except (TypeError,IndexError) as e:
    print(f"Invalid Input: {e}")
    end_time = time.time()
    logger(start_time,end_time,f"{e}")