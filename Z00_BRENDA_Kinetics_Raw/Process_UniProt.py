import os
import pandas as pd
import csv
import re

# File list
files = ["uniprot_annotation_3_001_200.tsv", "uniprot_annotation_3_201_300.tsv", "uniprot_annotation_3_301_400.tsv", "uniprot_annotation_3_401_600.tsv", "uniprot_annotation_3_601_800.tsv", "uniprot_annotation_3_801_inf.tsv", "uniprot_annotation_4.tsv", "uniprot_annotation_5.tsv"]

# Create output file
with open('output_file.tsv', 'w', newline='') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    # Write header
    tsv_output.writerow(['Sequence', 'Reaction', 'Rhea'])

    for file in files:
        if os.path.exists(file):
            chunksize = 10 ** 6  # adjust chunksize according to your available memory
            for i, chunk in enumerate(pd.read_csv(file, delimiter='\t', chunksize=chunksize)):
                print(f'Processing chunk {i+1} from {file}...')
                for index, row in chunk.iterrows():
                    # Convert to string before regex
                    catalytic_activity = str(row['Catalytic activity'])

                    sequence = row['Sequence']

                    # Make sure catalytic_activity is a string before searching with regex
                    if isinstance(catalytic_activity, str):
                        # Check if reaction exists in the row
                        reactions = re.findall('Reaction=(.*?);', catalytic_activity)
                        rheas = re.findall('Xref=Rhea:(.*?);', catalytic_activity)
                        


                        if reactions:
                            for reaction in reactions:
                                
                                rhea = next((r for r in rheas), None) # find the corresponding Rhea, if not available, None will be used
                                # Write to the output file
                                tsv_output.writerow([sequence, reaction, rhea])
                print(f'Finished chunk {i+1} from {file}...')
            print(f'Finished processing {file}!')
