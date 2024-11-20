# Inv_cooccurrence_Mguttatus

# Title: Are you with me? Co-occurrence tests from community ecology can identify positive and negative epistasis between inversions in Mimulus guttatus 

# Abstract

Chromosomal inversions are structural genetic variants that can play a crucial role in adaptive evolution and speciation. Patterns of attraction and repulsion among inversion orientations (alleles)—whether they tend to be carried by the same or different individuals—can indicate how selection is acting on these polymorphisms. In this study, we compare analytical techniques using data from 64 inversions that segregate among 1373 F2 plants of the yellow monkeyflower Mimulus guttatus.  Mendelian assortment provides a strong null hypothesis for χ^2 contingency tests.  Here, we show how co-occurrence metrics used in community ecology can provide additional insight regarding attraction and repulsion. The Jaccard-Tanimoto index and the affinity score describe the specific way that inversions interact to generate epistasis for plant survival. We further explore the use of network analysis to visualize inversion interactions and to identify essential nodes based on centrality measures, and the relation between the segregation distortion of each inversion and the likelihood of being involved in potentially repulsions or attractions. We suggest that a combination of different statistics will provide the most complete characterization of the fitness effects, both for inversions and other polymorphisms essential to adaptation and speciation.

# Cloning the repository

Generic link: https://github.com/luismadrigal98/Inv_cooccurrence_Mguttatus

ssh protocol: git clone git@github.com:luismadrigal98/Inv_cooccurrence_Mguttatus

This repository contains the code and data analysis scripts for the study on chromosomal inversion o-occurrence patterns in Mimulus guttatus. The study compares analytical techniques to investigate patterns of attraction and repulsion among inversion orientations.

In order to replicate the analysis, the only thing you need to do is defining the work directory in the main script, or running such script from the command line by typing `Rscript ./main.r` from the cloned directory.

# Repository structure:

.
├── Data
│   ├── CSVs
│   │   ├── Inversion_status_data_1034.csv
│   │   ├── Inversion_status_data_1192.csv
│   │   ├── Inversion_status_data_155.csv
│   │   ├── Inversion_status_data_444.csv
│   │   ├── Inversion_status_data_502.csv
│   │   ├── Inversion_status_data_541.csv
│   │   ├── Inversion_status_data_62.csv
│   │   ├── Inversion_status_data_664.csv
│   │   └── Inversion_status_data_909.csv
│   └── inv_and_gene_metadata.csv
├── LICENSE
├── Python_scripts
│   ├── 00_Data_puller.py
│   ├── 00_Data_puller_global.py
│   ├── GO_term_correction.py
│   └── Gene_renamer.py
├── README.md
├── R_scripts
│   ├── Aux1_Contingency_tables_simulation.R
│   ├── Aux2_Segregation_distortion_analysis.R
│   ├── Aux3_Co-occurence_per_line.R
│   ├── Aux4_Individual_effects_and_co-occurrence_patterns.R
│   ├── Dictionary_of_inversions.R
│   └── Extreme_but_not_significant_exploration_alpha.R
├── Results
│   └── Plots
│       └── Networks
├── main.R
└── src
    ├── G_calculator.R
    ├── G_per_contrast_per_line.R
    ├── X2_sub1_calculator.R
    ├── X2_sub2_calculator.R
    ├── affinity_filter.R
    ├── alpha_values_extractor.R
    ├── are_different_chr.R
    ├── check_zero_marginals.R
    ├── compare_all_pairs.R
    ├── compare_networks.R
    ├── contingency_expanded.R
    ├── contingency_filter.R
    ├── count_inversions.R
    ├── create_chr_lookup.R
    ├── create_combined_dataframe.R
    ├── deviant_signal_summarizer.R
    ├── dosage_splitter.R
    ├── expected_mod.R
    ├── is.valid.R
    ├── jaccard_tanimoto_condenser.R
    ├── marginal_compliance_checking.R
    ├── mask_chromosomes.R
    ├── observed_count_extractor.R
    ├── p_corrector.R
    ├── p_value_extractor.R
    ├── perform_heterogeneity_analysis.R
    ├── plot_network.R
    ├── plot_network_lite.R
    ├── process_cols.R
    ├── process_submatrix_tests.R
    ├── pull_aux_functions.R
    ├── replicated_G_test.R
    ├── select_frequent_inversions.R
    ├── set_environment.R
    ├── sig_network_builder.R
    ├── sig_network_builder_lite.R
    ├── significant_X2_counter.R
    ├── submatrix_extractor.R
    ├── support_counter.R
    ├── table_sim_null.R
    ├── two_vectors_from_submatrix.R
    ├── x2_calculator.R
    ├── x2_posthoc_calculator.R
    └── x2p_to_square.R

8 directories, 67 files