""""
This script will be generate the design matrix from the genotype files.

@Author: Luis Javier Madrigal Roca & John K. Kelly
@Date: 2024-08-31

"""

import os
import pandas as pd

GEN_DIR = '/home/l338m483/scratch/Cooccurrence_Inv/new_geno_files/' # Directory where the genotype files are located
CSV_DIR = '/home/l338m483/scratch/Cooccurrence_Inv/CSVs/' # Directory where the CSV files will be saved
INV_LIST = '/home/l338m483/scratch/Cooccurrence_Inv/R_directory/list.of.inversions.txt'

# List of the genotype files
gen_files = os.listdir(GEN_DIR)
dictionary_file = INV_LIST

Inversions_per_family = {}

with open(dictionary_file, "r") as file:
    for i, line in enumerate(file):
        
        line = line.strip()
        
        if i > 0:
            inv = line.split('\t', 2)[0]
            families = line.split('\t', 2)[2]
            
            if ";" not in families:
                family = families
            
                if family not in Inversions_per_family.keys():
                    Inversions_per_family[family] = []
            
                Inversions_per_family[family].append(inv)
            
            elif ";" in families:
                families = families.split(";")
                
                for family in families:
                    
                    if family not in Inversions_per_family.keys():
                        Inversions_per_family[family] = []
                    
                    Inversions_per_family[family].append(inv)

            else:
                print("There was an error in the dictionary file. Please check the input file.")                

# 38.genotype.txt

Inversion_status = {}

for family in Inversions_per_family.keys():
    
    Inversion_status_data = pd.DataFrame()

    for file in gen_files:
        
        if "genotype" in file:
            key = file.rsplit(".", 2)[0]
            if key in Inversions_per_family[family]:
                con = open(GEN_DIR + file)
                lines = con.readlines()
                
                Inversion_status = {}
                Inversion_status['Probes'] = []
                Inversion_status['Inv_' + key] = []

                for line in lines:
                    line = line.strip()
                    line = line.split('\t')

                    if line[1].rsplit("-", 1)[0].split("_")[1] == family:
                        # 186	s3_664-P11	664	three	parent	2
                        if line[4] == "good_f2":
                            Inversion_status['Probes'].append(line[1])
                            Inversion_status['Inv_' + key].append(line[5])
                        else:
                            continue
                    else:
                        continue
                
                df = pd.DataFrame(Inversion_status)
                if Inversion_status_data.empty:
                    Inversion_status_data = df
                else:
                    Inversion_status_data = pd.merge(Inversion_status_data, df, on='Probes')

                con.close()

    Inversion_status_data.to_csv(f"{CSV_DIR}Inversion_status_data_{family}.csv", index=False)