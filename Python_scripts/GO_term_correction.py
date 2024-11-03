"""
This script will change the genes name in the annotation file that will be used for the GO analysis

"""

## Line in info file from the annotation of MG (tab-delimited):
## PAC:64874943	MgIM767.01G000300	MgIM767.01G000300.1	MgIM767.01G000300.1.p		PTHR12663 PTHR12663:SF3				GO:0005488	AT4G31880	(1 of 3) PTHR12663:SF3 - ASPARTYL BETA-HYDROXYLASE N-TERMINAL REGION DOMAIN-CONTAINING PROTEIN-RELATED			LOC_Os04g42320	AT hook motif family protein, expressed (LOC_Os04g42320), mRNA.	LOC_Os04g42320

## Line in info file from the annotation of MG (tab-delimited):
## GO:0034309	primary alcohol biosynthetic process	THI1	MIMGU_mgv1a002581mg	MIMGU_mgv1a004117mg	MIMGU_mgv1a002082mg	MIMGU_mgv1a002684mg	MIMGU_mgv1a002040mg	MIMGU_mgv1a002143mg	MIMGU_mgv1a021684mg

ANN_FILE = '/home/l338m483/scratch/Corrected_inv_genotype_files/Directory/GO_annotations/BACKUP/MY_GO_terms_build/GO_annotation_MG.gmt'
INFO_FILE = '/home/l338m483/scratch/Corrected_inv_genotype_files/Directory/GO_annotations/BACKUP/MY_GO_terms_build/Mguttatusvar_IM767_887_v2.1.annotation_info.txt'
OUT_FILE = '/home/l338m483/scratch/Corrected_inv_genotype_files/Directory/GO_annotations/BACKUP/MY_GO_terms_build/GO_annotation_MG_corrected.gmt'

## Read the info file and store the GO terms as keys of a dictionary and the genes as values

GO_genes = {}
global GO_column

with open(INFO_FILE, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        line = line.strip().split('\t')
        gene = line[1]
        
        for item in line:
            if "GO:" in item:
                GO_terms = item.split(' ')
                
                for GO in GO_terms:
                    if GO.startswith("GO:"):  # Check if the term starts with "GO:"
                        if GO not in GO_genes:
                            GO_genes[GO] = []

                        GO_genes[GO].append(gene)

## Read the annotation file and change the gene names

with open(ANN_FILE, 'r') as f, open(OUT_FILE, 'w') as out:
    for line in f:
        line = line.strip().split('\t', 2)
        GO = line[0]
        description = line[1]
        if GO in GO_genes:
            genes = GO_genes[GO]
        else:
            genes = []
        out.write(GO + '\t' + description + '\t' + '\t'.join(genes) + '\n')