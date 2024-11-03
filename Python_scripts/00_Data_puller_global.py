""""
This script will be generate the design matrix from the genotype files, merging all the information in the same data frame.

@Author: Luis Javier Madrigal Roca & John K. Kelly
@Date: 2024-07-30

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

with open(os.path.join(GEN_DIR, dictionary_file), "r") as file:
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

Inversion_status_data = pd.DataFrame()

for file in gen_files:
    if "genotype" in file:
        with open(os.path.join(GEN_DIR, file), "r") as con:
            key = file.rsplit(".", 2)[0]
            lines = con.readlines()
                
            Inversion_status = {}
            Inversion_status['Plant'] = []
            Inversion_status['INV_' + key] = []

            for line in lines:
                line = line.strip()
                line = line.split('\t')
                
                family = line[1].rsplit("-", 1)[0].split("_")[1]

                # 186	s3_664-P11	664	three	parent	2
                Inversion_status['Plant'].append(line[1])
                try:
                    if key in Inversions_per_family[family]: # If the inversion is segregating in that family
                        Inversion_status['INV_' + key].append(line[5]) # Go ahead and add the genotype information
                    else:
                        Inversion_status['INV_' + key].append(0) # If the inversion is not segregating in that family, add a 0 (no copies of the inversion)
                except KeyError:
                    Inversion_status['INV_' + key].append(0) # Catch for the 767 family, which is not in the dictionary file because there is no inversion segregating in that family
            df = pd.DataFrame(Inversion_status)
            if Inversion_status_data.empty:
                Inversion_status_data = df
            else:
                Inversion_status_data = pd.merge(Inversion_status_data, df, on='Plant')

Inversion_status_data.to_csv(f"{CSV_DIR}/Inversion_status_global.csv", index=False)
print(Inversion_status_data)