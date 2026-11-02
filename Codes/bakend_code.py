from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq,transcribe,translate
import time
from datetime import datetime
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import warnings
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#Logger
def logger(start_time,end_time,update):

    log_file = open("log_file.txt","a")
    log_file.write(f"Log Status on: {datetime.now().date()} at {datetime.now().time()}\n")
    log_file.write(f"Total Time Taken: {(end_time - start_time):.4f}\n")
    log_file.write(f"Update: {update}\n\n\n")
    log_file.close()
    return None

#Checking the reference seq
class coding_seq:
    #initialising the sequence variable:
    def __init__(self,dna_seq):
        self.dna_seq=dna_seq
        return None
    #Sequence checker
    def seq_check(self,dna_seqs=None, seq_type='Provided'):
        #Extracting the sequence form the BioSeq object
        print(f"\n\nCHECKING THE {seq_type.upper()} SEQUENCE\n\n")
        if dna_seqs is None:
            dna_seq=self.dna_seq.seq.upper()
        else:
            dna_seq=dna_seqs.seq.upper()
        #Extracting the length of the seq
        dna_seq_len = len(dna_seq)
        #Checks if it contains all ATGC values
        pos=[]
        for i in range(dna_seq_len):
            if dna_seq[i] not in ['A','T','C','G']:
                pos.append(i+1)
        if pos !=[]:
            raise ValueError(f"{seq_type} Sequence contains non-DNA Charachters at positions {pos}")

        #Checks if the sequence is a multiple of 3 or not
        if dna_seq_len%3!=0:
            raise ValueError(F"The {seq_type} sequence is not a CDS as it lacks codon multiples.")
            
        #Checks for Start codon and Stop codon
        val=[]
        #Checking for Start Codon
        if dna_seq[:3] !='ATG':
            val.append(False)
        else:
            val.append(True)

        #Checking for Stop codon
        pos=[]
        #Internal stop codon checker
        for i in range(0,dna_seq_len,3):
            if dna_seq[i:i+3] in ['TAA','TGA','TAG']:
                pos.append(i+1)
        
        #Stop codon checker
        if dna_seq[-3:] in ['TAA','TGA','TAG'] and len(pos) == 1:
            val.append(True)
        else:
            val.append(False)
            
        #Value Checker
        #For Start Codon
        if val[0] == False:
            raise ValueError(f"The {seq_type} Seq is missing start codon 'ATG'")
        #For Stop Codon
        elif val[1] == False and len(pos) ==0: 
            raise ValueError(f"The {seq_type} sequence does not contain stop codon 'TAA', 'TGA', or 'TAG'")
        elif val[1] == False and len(pos)!=0:
            #checking if the last value of position list contains ending stop codon or not
            if pos[-1]!=dna_seq_len-2:
                raise ValueError(f"The {seq_type} sequence contains internal stop codon 'TAA', 'TGA', or 'TAG' at position {pos}")
            else:
                raise ValueError(f"The {seq_type} sequence contains internal stop codon 'TAA', 'TGA', or 'TAG' at position {pos[:-1]}")
        else:
            True
        #returns the extracted dna_seq
        return dna_seq

    #base metric calculator for validated sequences
    #1) Length of the given seq
    #2) GC% of the given sequence
    def base_cal(self):#,dna_seq
        #importing necessary 
        from Bio import SeqIO
        from Bio.SeqUtils import gc_fraction
        print("\n\nPERFORMING BASIC PARAMETER CALCULATION\n\n")
        #reading the sequence
        records = list(SeqIO.parse(self.dna_seq,'fasta'))
        #for each record
        dna_len_list=[]
        gc_frac_list=[]
        for record in records:
            dna_len_list.append(len(record.seq))
            gc_frac_list.append(gc_fraction(record.seq)*100)
        return dna_len_list, gc_frac_list

    #Fasta File checker
    def fasta_checker(self):#,dna_seq_path
        #importing necessary packages
        print("\n\nCHECKING THE PROVIDED FASTA FILE\n\n")
        dna_seq_path = self.dna_seq
        #Checking for fasta value
        extension = dna_seq_path.split(".")[1]
        if extension not in ['fasta','fa','fas']:
            raise ValueError(F"The provided file is not a fasta file. It is `.{extension}` file")

        #Checking the no of sequences in fasta file
        records = list(SeqIO.parse(dna_seq_path,'fasta'))
        if len(records)>1:
            raise ValueError(F"More than one reference sequences in the provided fasta file.")
        #Reference Seq Checker:
        chekced_dna_seq=self.seq_check(dna_seqs=records[0],seq_type='Reference')
        #seq_check(dna_seq=records[0])
        return chekced_dna_seq


class variant_analyzer:
    def __init__(self,csv_file):
        self.csv_file = csv_file
    def varaint_csv_checker(self,dna_seq):
        print("\n\nCHECKING THE PROVIDED CSV FILE\n\n")
        # import warnings
        # import pandas as pd
        #reading the csv file
        csv_ext=self.csv_file.split(".")[1]
        if csv_ext not in ['csv','CSV']:
            raise TypeError(f" The provided `{self.csv_file}` file is not a `.csv` file")
        variant = pd.read_csv(self.csv_file)

        #Checking the columns inside the vraint file
        acceptable_cols=['variant_id','type','cds_position','change_from','change_to']
        csv_cols = variant.columns.str.strip().str.lower()
        variant.columns = csv_cols
        
        if len(csv_cols)<len(acceptable_cols):
            missing_cols=[]
            for col in acceptable_cols:
                if col not in csv_cols:
                    missing_cols.append(col)
            raise IndexError(f"The provided `{self.csv_file}` file is missing the following columns {missing_cols}")

        #Checking the unique variants
        updated_variant=None
        rem_variant=None
        if len(variant['variant_id'].unique()) == variant.shape[0]:
            print()
        else:
            #removing the variants with duplicate variant id
            rem_variant= variant[variant.duplicated(subset=variant.columns[0],keep='first')]
            #seperating the unique variant ids
            updated_variant = variant.drop_duplicates(subset=variant.columns[0],keep='first')
            warnings.warn(f"The provided `{self.csv_file}` contains duplicate variant ids. They will be removed from furthere analysis.")
            #return rem_variant,updated_variant

        #Checking for the type of mutations
        
        #checking in the original variant
        if updated_variant is None:
            type_index=[]
            for i in list(variant['cds_position'].index):
                if  variant['type'][i].strip().lower() not in ['substitution','deletion','insertion']:
                    type_index.append(i)

            for i in type_index:
                        rem_variant = pd.concat([rem_variant,variant.loc[[i]]],axis=0)
                        updated_variant=variant.drop(i)
        
        #checking in the updated variant
        else:
            type_index=[]
            for i in list(updated_variant['cds_position'].index):
                if  updated_variant['type'][i].strip().lower() not in ['substitution','deletion','insertion']:
                    type_index.append(i)
            for i in type_index:
                #print(updated_variant.iloc[[i]])
                rem_variant = pd.concat([rem_variant,updated_variant.loc[[i]]],axis=0)
                updated_variant=updated_variant.drop(i)

        if type_index is not []:
            warnings.warn(f"The provided `{self.csv_file}` contains invalid types. They will be removed from furthere analysis.")
        
        
        #Checking for CDS positions
        
        #checking in the original variant
        if updated_variant is None:
            pos_index=[]
            for i in list(variant['cds_position'].index):
                if  variant['cds_position'][i] > len(dna_seq):
                    pos_index.append(i)

            for i in pos_index:
                rem_variant = pd.concat([rem_variant,variant.loc[[i]]],axis=0)
                updated_variant=variant.drop(i)
        
        #checking in the updated variant
        else:
            pos_index=[]
            for i in list(updated_variant['cds_position'].index):
                if  updated_variant['cds_position'][i] > len(dna_seq):
                    pos_index.append(i)

            for i in pos_index:
                #print(updated_variant.iloc[[i]])
                rem_variant = pd.concat([rem_variant,updated_variant.loc[[i]]],axis=0)
                updated_variant=updated_variant.drop(i)

        if pos_index is not []:
            warnings.warn(f"The provided `{self.csv_file}` contains invalid cds positions. They will be removed from furthere analysis.")



        #Checking for changes to be made
        #checking in the original variant
        if updated_variant is None:
            pos_index=[]
            for i in list(variant['change_from'].index):
                if  variant['change_from'][i].strip().lower() not in ['a','t','c','g','delete','insert']:
                    pos_index.append(i)

            for i in pos_index:
                rem_variant = pd.concat([rem_variant,variant.loc[[i]]],axis=0)
                updated_variant=variant.drop(i)
        
        #checking in the updated variant
        else:
            pos_index=[]
            for i in list(updated_variant['change_from'].index):
                if  updated_variant['change_from'][i].strip().lower() not in ['a','t','c','g','delete','insert']:
                    pos_index.append(i)

            for i in pos_index:
                rem_variant = pd.concat([rem_variant,updated_variant.loc[[i]]],axis=0)
                updated_variant=updated_variant.drop(i)

        if pos_index is not []:
            warnings.warn(f"The provided `{self.csv_file}` contains invalid changes. They will be removed from furthere analysis.")

        #Checking for alleles
        #checking in the original variant
        if updated_variant is None:
            pos_index=[]
            for i in list(variant['change_to'].index):
                if  len(variant['change_to'][i].strip().lower())<2 and variant['change_to'][i].strip().lower() not in ['a','t','c','g','n']:
                    pos_index.append(i)
                elif len(variant['change_to'][i].strip().lower())==3 :
                    for val in variant['change_to'][i].strip().lower():
                        if val not in 'atgcn':
                            pos_index.append(i)
                else:
                    continue

            for i in pos_index:
                #print(i)
                rem_variant = pd.concat([rem_variant,variant.loc[[i]]],axis=0)
                updated_variant=variant.drop(i)
            
        #checking in the updated variant
        else:
            pos_index=[]
            for i in list(updated_variant['change_to'].index):
                if  len(variant['change_to'][i].strip().lower())<2 and variant['change_to'][i].strip().lower() not in ['a','t','c','g','n']:
                    
                    pos_index.append(i)
                elif len(variant['change_to'][i].strip().lower())==3 :
                    for val in variant['change_to'][i].strip().lower():
                        if val not in 'atgcn':
                            pos_index.append(i)
                else:
                    continue

            for i in pos_index:
                rem_variant = pd.concat([rem_variant,updated_variant.loc[[i]]],axis=0)
                updated_variant=updated_variant.drop(i)

        if pos_index is not []:
            warnings.warn(f"The provided `{self.csv_file}` contains invalid alleles. They will be removed from furthere analysis.")
        

        #checking if the updated variant is empty or not
        if updated_variant is not None and updated_variant.empty or variant.empty :
            raise ValueError(f" The provided `{self.csv_file}` file does not contain any valid variants.")
        
        #retruning the dropped files, updated files or rectified files
        if updated_variant is None:
            return rem_variant,variant
        else:
            return rem_variant,updated_variant

    #Perfomring Mutation
    def mutator(self,variant,dna_seq,output_dir="./"):
        #creating a global mutation list which will be used to make fasta file
        print("\n\nPERFORMING MUTATIONS\n\n")
        mutations =[]
        var_id=[]
        ref_seq_len = len(dna_seq)
        seq_len = []
        mutation_type=[]

        #performing mutations for substituition:
        sub_var = variant[variant['type']=='Substitution']
        #substituting the given positions:
        for idx in list(sub_var.index):
            pos = sub_var['cds_position'].loc[idx]
            sub = sub_var['change_to'].loc[idx]
            all_to_sub = sub_var['change_from'].loc[idx]
            if all_to_sub != dna_seq[pos-1]:
                warnings.warn(f"The provided allele does not match the ref seq at position {pos}.\n Proceeding with subtitution on the basis CDS Position")
            mutated_seq= dna_seq[:pos-1] + sub + dna_seq[pos:]           
            mutations.append(mutated_seq)
            var_id.append(sub_var['variant_id'].loc[idx])
            mutation_type.append("subtitution")
            seq_len.append(len(mutated_seq))
        #print(mutations)


        #performing mutations for insertions
        ins_var = variant[variant['type']=='Insertion']
        #substituting the given positions:
        for idx in list(ins_var.index):
            pos = ins_var['cds_position'].loc[idx]
            ins = ins_var['change_to'].loc[idx]
            mutated_seq= dna_seq[:pos] + ins + dna_seq[pos+1:]           
            mutations.append(mutated_seq)
            var_id.append(ins_var['variant_id'].loc[idx])
            mutation_type.append("insertion")
            seq_len.append(len(mutated_seq))
        #print(mutations)


        #perfoorming mutations for deletions
        del_var = variant[variant['type']=='Deletion']
        #substituting the given positions:
        for idx in list(del_var.index):
            pos = del_var['cds_position'].loc[idx]
            all_to_del = del_var['change_to'].loc[idx]
            if all_to_del != dna_seq[pos-1:pos+len(all_to_del)-1]:
                warnings.warn(f"The provided allele does not match the ref seq at position {pos}.\n Proceeding with deletion on the basis CDS Position")        
            mutated_seq= dna_seq[:pos] + dna_seq[pos+1:]           
            mutations.append(mutated_seq)
            var_id.append(del_var['variant_id'].loc[idx])
            mutation_type.append("deletion")
            seq_len.append(len(mutated_seq))
        print(mutations)

        #Creating a summary table
        print("\n\nCREATING MUTATION SUMMARY TABLE\n\n")
        sum_table = pd.DataFrame({
            'variant_id':var_id,
            'mutation_type':mutation_type,
            'reference_seq_length':[ref_seq_len for i in range(len(var_id))],
            'mutated_seq_length':seq_len
        })
        print(sum_table)
        sum_table.to_csv(f"{output_dir}mutation_summary.csv",index=False)
        print("\n\nCREATED THE MUTATION SUMMARY TABLE:- stored as `mutation_summary.csv`")

        #Creating the mutated fast file
        print("\n\nCREATING THE MUTATED FASTA FILE\n\n")
        records =[]
        for idx in range(len(mutations)):
            records.append(SeqRecord(seq=Seq(mutations[idx]),id=var_id[idx],description=f"{mutation_type[idx]} with sequence length of {len(mutations[idx])}"))
        SeqIO.write(records,f"{output_dir}mutations.fas","fasta")
        print("\n\nCREATED THE MUTATIONS FASTA FILE:- stored as `mutations.fas`")

    def mut_seq_checker(self,seq,ref_seq_len):
        #checking for early stop codons and no start codons
        
        checks=[]
        val=[]
        #Checking for Start Codon
        if seq[:3] !='ATG':
            checks.append('No Start Codon')
            
        else:
            checks.append("Start Codon Present")


        #Checking for Stop codon
        pos=[]
        #Internal stop codon checker
        for i in range(0,len(seq),3):
            if seq[i:i+3] in ['TAA','TGA','TAG']:
                pos.append(i+1)
        if len(pos)>1:
            checks.append(f"Internal stop codon present:{pos}")
        else:
            checks.append("No Internal Stop Codon")
        
        #Stop codon checker
        if seq[-3:] in ['TAA','TGA','TAG'] and len(pos) == 1:
            checks.append("Stop Codon Present")
            pos.pop()
        else:
            checks.append("No stop codon")

        #Checking for frameshifts
        if len(seq) > ref_seq_len:
            if len(seq)%3==0: 
                checks.append("Inframe")
            else:
                checks.append("Frameshift")

        elif len(seq) < ref_seq_len:
            if len(seq)%3==0:
                checks.append("Inframe")
            else:
                checks.append("Frameshift")
        else:
            checks.append("Equal")
        return checks

    def translations(self,dna_seq,mutations,output_dir="./"):
        #performing translation and getting theri consequences
        print("\n\nPERFORMING TRANSLATION OF THE PROVIDED MUTATED SEQUENCES\n\n")
        seqs=[dna_seq]
        records = list(SeqIO.parse(mutations,"fasta"))
        seqs.extend(records)
        #print(seqs)

        #translations
        trans_seq =[translate(dna_seq)]
        variant_id = ['ref_seq']
        mutation_type=['original']
        ref_seq_len = len(dna_seq)
        seq_len=[ref_seq_len]
        errors=['-',]
        for seq in seqs[1:]:
            #extracting id
            variant_id.append(seq.id)
            #extraction mutation type
            mut_type = seq.description.split(" ")[1]
            mutation_type.append(mut_type)
            #extracting mutation length
            mut_len = seq.description.split(" ")[-1]
            seq_len.append(mut_len)
            #performing translaitonss
            seq = seq.seq.upper()
            checks = self.mut_seq_checker(seq,ref_seq_len)
            errors.append(checks)
            #tranlating the sequence
            trans_seq.append(translate(seq))
        print(trans_seq)
        #print(errors)

        #Creating the protein fasta file
        print("\n\nCREATING THE TRANSLATED FASTA FILE\n\n")
        records =[]
        for idx in range(len(trans_seq)):
            seq = trans_seq[idx]
            id  = seq.find("*")
            seq = seq[:id]
            records.append(SeqRecord(seq=seq,id=variant_id[idx],description=f"{mutation_type[idx]} with sequence length of {seq_len[idx]}"))
        SeqIO.write(records,f"{output_dir}/proteins.fasta","fasta")
        print("\n\nCREATED THE PROTEIN FASTA FILE:- stored as `proteins.fasta`")

        #Creating summary Dictionary
        sum_df ={
            'Variant_Id':variant_id,
            'Mutation_Type':mutation_type,
            'Sequence_Length':seq_len,
            'Ref_Protein_Length':[len(trans_seq[0]) for i in range(len(variant_id))],
            'Mutated_Protein_Length':[len(trans_seq[i]) for i in range(len(variant_id))],
            "Difference_Protein_Length":[abs(len(trans_seq[0])-len(trans_seq[i])) for i in range(len(variant_id))],
            # 'Some_Checks':errors
        }
        #print(pd.DataFrame(sum_df))
        #print(sum_df)
        #Calculating the consequences for the mutations
        print("\n\nCALCULATING CONSEQUENCES\n\n")
        og_protein = trans_seq[0]
        og_protein_len = len(og_protein)
        consequences=[["-","-","-","-","-"]]
        for idx in range(len(trans_seq[1:])):
            #checking changes in the amino acid:
            consequence=[]
            if og_protein[-1]!=trans_seq[idx+1][-1]:
                consequence.append(f"Stop-Loss:{trans_seq[idx+1][-1]} at position {i+1}")
            for i in range(og_protein_len):
                if og_protein[i]=="*":
                    break
                if og_protein[i] != trans_seq[idx+1][i]:
                    if trans_seq[idx][i] !="*":
                        consequence.append(f"Missense:{trans_seq[idx+1][i]} at position {i+1}")
                    else:
                        consequence.append(f"Non-sense:{trans_seq[idx+1][i]} at position {i+1}")
                else:
                    continue
            if consequence ==[]:
                consequence.append("Synonyms at every position.")

            #Insertion and deletion
            error = errors[idx+1]
            for val in error:
                consequence.append(val)
            # print(len(consequence))
            # print(consequence)
            consequences.append(consequence)

        #Creating the consequence dataframe:
        for conseq in consequences:
            sum_df['Substituition_effect'] = conseq[0]
            sum_df['Start_Codon'] = conseq[1]
            sum_df['Internal_Stop_Codon'] = conseq[2]
            sum_df['Stop_Codon'] = conseq[3]
            sum_df['Shifts'] = conseq[4]

        sum_df = pd.DataFrame(sum_df)
        print(sum_df)                      
        print("\n\nSAVING THE CONSEQUENCE DATAFRAME:- stored as mutation_effect.csv\n\n")              
        sum_df.to_csv(f"{output_dir}mutation_effect.csv",index=False)                                

    #Creating protein properites
    def protein_prop(self,protein_fasta,output_dir="./"):
        #getting the records
        print("\n\nCALCULATING THE PROTEIN PROPERTIED\n\n")
        records = list(SeqIO.parse(protein_fasta,"fasta"))
        #analysing
        mol_wt=[]
        length=[]
        pI=[]
        gravy =[]
        composition_cont=[]
        composition_perc=[]
        self_composition_cont=[]
        self_composition_perc=[]
        diff_comp_conts=[]
        diff_comp_percs=[]
        var_id=[]
        amino_acids = "ACDEFGHIKLMNPQRSTVWY"
        for record in records:
            var_id.append(record.id)
            seq =record.seq
            analysis = ProteinAnalysis(seq)
            length.append(len(seq))
            mol_wt.append(analysis.molecular_weight())
            pI.append(analysis.isoelectric_point())
            gravy.append(analysis.gravy())
            comp_cont_bio = analysis.count_amino_acids()
            #print(comp_cont_bio)
            composition_cont.append(analysis.count_amino_acids())
            comp_perc_bio={}
            for i in amino_acids:
                comp_perc_bio[i]=comp_cont_bio[i]/len(seq)*100
            composition_perc.append(comp_perc_bio)
            #self calucaltion of amino acids
            comp={}
            for i in seq:
                if i in comp:
                    comp[i] = comp[i]+1
                else:
                    comp[i] = 1
            self_composition_cont.append(comp)
            comp_per={}
            for i in comp.keys():
                comp_per[i] = comp[i]/len(seq)*100
            self_composition_perc.append(comp_per)
            #differnce in self_and biopytohon cont
            diff_comp_cont={}
            diff_comp_perc={}
            for i in comp.keys():
                diff_comp_cont[i]= comp_cont_bio[i]-comp[i]
                
                diff_comp_perc[i] = comp_perc_bio[i]-comp_per[i]
            diff_comp_conts.append(diff_comp_cont)
            diff_comp_percs.append(diff_comp_perc)
            

        sum_df={
            "Varient_id":var_id,
            "Protein_Length":length,
            "Mol_Wt":mol_wt,
            "pI":pI,
            "gravy":gravy,
            "AA_Composition_Count":composition_cont,
            "Self_AA_Compoisiotn_Count":self_composition_cont,
            "Difference_cont":diff_comp_conts,
            "AA_Composition_Perc":composition_perc,
            "Self_AA_Composition_Perc":self_composition_perc,
            "Differnce_AA_Compositon_Perc":diff_comp_percs
        }
        sum_df=pd.DataFrame(sum_df)
        print(sum_df)
        print("CALCULATED AND SAVING THE PROTIN PROTEINS DF:- stored as prot_prop.csv\n\n")
        sum_df.to_csv(f"{output_dir}prot_prop.csv",index=False)