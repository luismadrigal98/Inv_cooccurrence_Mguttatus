"""
This script will modify the gene names in a set of text files based on a mapping provided in an Excel file.

"""

import pandas as pd
import os

# Read the Excel file
df = pd.read_excel('/home/l338m483/scratch/Corrected_inv_genotype_files/new_name_old_name.genes.xlsx')

# Create an structure for saving the new names, the inversion ID, and other information
result = {}

# Create a dictionary mapping old names to new names
name_mapping = df.set_index('old_name')['new_name'].to_dict()

# List the files in the directory

geno_dir = '/home/l338m483/scratch/Corrected_inv_genotype_files/'

# Example of a first line
# 47      Chr_13:1188673,1234110  17

# Example of a second line:
# 47      genes_included:['MiIM7v11033888m.g', 'MiIM7v11033889m.g', 'MiIM7v11033891m.g', 'MiIM7v11033892m.g', 'MiIM7v11033894m.g', 'MiIM7v11033896m.g', 
#'MiIM7v11033897m.g', 'MiIM7v11033898m.g', 'MiIM7v11033899m.g', 'MiIM7v11033901m.g', 'MiIM7v11033904m.g', 'MiIM7v11033909m.g', 'MiIM7v11033910m.g', 
#'MiIM7v11033911m.g', 'MiIM7v11033912m.g', 'MiIM7v11033913m.g', 'MiIM7v11033915m.g']

for file in os.listdir(geno_dir):
    if file.endswith('genotype.txt'):
        # Read the file
        with open(os.path.join(geno_dir, file), 'r') as f:
            lines = f.readlines()
        
        first_line = [element.strip() for element in lines[0].split('\t')]
        INV_ID = first_line[0]
        INV_loc = first_line[1]
        INV_gene_number = first_line[2]

        # Get the gene names
        gene_names = [name.replace("'", "") for name in lines[1].split('[')[1].split(']')[0].split(', ')]

        # Modify the gene names
        new_gene_names = []
        for gene in gene_names:
            new_gene_names.append(name_mapping[gene])

        # Join the new gene names into a single string
        new_gene_names_str = ', '.join(new_gene_names)

        # Save the new names and the complementary information
        result[INV_ID] = {'INV_ID': INV_ID, 'INV_loc': INV_loc, 'INV_gene_number': INV_gene_number, 'new_gene_names': new_gene_names_str}

# Save the result to a csv file
df_result = pd.DataFrame(result).T
df_result.to_csv('/home/l338m483/scratch/Corrected_inv_genotype_files/renamed_genes.csv', index=False)