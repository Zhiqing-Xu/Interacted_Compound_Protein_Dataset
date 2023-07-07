import csv
import os

input_file = 'Z03B_CompoundName_LigandID_SMILES.csv'
temp_file = 'Z03B_CompoundName_LigandID_SMILES_1.csv'

with open(input_file, 'r', newline='', encoding='utf-8') as infile:
    reader = csv.reader(infile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL, skipinitialspace=True)

    with open(temp_file, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.writer(outfile, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        for row in reader:
            writer.writerow(row[1:])  # Remove the first column (index)

# Replace the input file with the updated temp_file
os.remove(input_file)
os.rename(temp_file, input_file)

print(f"First column removed from {input_file} and separator updated to ';'.")
