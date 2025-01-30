"""
This program consolidates genotypic information of inversions analyzed across multiple crosses into a single table.
It combines individual CSV files per cross and assigns 0 to inversions not present in a given cross.

Final result is a tab-delimited file with the following columns:
- Cross
- Probe
- Columns for each inversion, coded as Inv_XX
"""

import os
import pandas as pd
import logging
import sys

def make_genotype_matrix(input_folder, output_file):
    """
    Combines genotype data from multiple CSV files into a single matrix.
    
    Args:
        input_folder (str): Path to folder containing input CSV files
        output_file (str): Path where output file will be saved
    
    Returns:
        pd.DataFrame: Combined genotype matrix
    """
    logging.basicConfig(level=logging.INFO)
    
    # Check if input folder exists
    if not os.path.exists(input_folder):
        raise FileNotFoundError(f"Input folder {input_folder} does not exist")
    
    # Get CSV files from input folder
    files = [f for f in os.listdir(input_folder) if f.endswith('.csv')]
    if not files:
        raise ValueError(f"No CSV files found in {input_folder}")
    
    # Initialize containers
    data_frames = []
    all_inversions = set()
    
    # Process each file
    for file in files:
        try:
            # Read data
            file_path = os.path.join(input_folder, file)
            df = pd.read_csv(file_path)
            
            # Extract cross name from filename
            try:
                cross = file.rsplit('_')[0].split('.')[0]  # Adjust splitting logic based on your filenames
                logging.info(f"Processing cross: {cross}")
            except IndexError:
                logging.error(f"Could not extract cross name from {file}")
                continue
            
            # Add Cross and ensure Probe column exists
            df['Cross'] = cross
            if 'Probes' not in df.columns:
                df['Probe'] = df.index
            else:
                df.insert(0, 'Probe', df['Probes'])
            
            # Update set of all inversions
            inversion_cols = [col for col in df.columns if col not in ['Cross', 'Probe']]
            all_inversions.update(inversion_cols)
            
            data_frames.append(df)
            
        except Exception as e:
            logging.error(f"Error processing file {file}: {str(e)}")
            continue
    
    if not data_frames:
        raise ValueError("No valid data frames were created")
    
    # Combine all data frames
    logging.info("Combining data frames...")
    results = pd.DataFrame(columns=['Cross', 'Probe'] + sorted(list(all_inversions)))
    
    for df in data_frames:
        # Ensure all inversion columns exist
        for inv in all_inversions:
            if inv not in df.columns:
                df[inv] = pd.NA
        
        # Reorder columns to match final format
        df = df[['Cross', 'Probe'] + sorted(list(all_inversions))]
        results = pd.concat([results, df], ignore_index=True)
    
    # Save results
    logging.info(f"Saving results to {output_file}")
    results.to_csv(output_file, sep="\t", index=False)

# Example usage
if __name__ == "__main__":
    input_folder = sys.argv[1]
    output_file = sys.argv[2]
    matrix = make_genotype_matrix(input_folder, output_file)