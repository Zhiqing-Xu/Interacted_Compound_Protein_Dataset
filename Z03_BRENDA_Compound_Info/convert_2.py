import csv

input_file  = 'Z03B_CompoundName_LigandID_SMILES.csv'
output_file = 'Z03B_CompoundName_LigandID_SMILES_updated.csv'

from rdkit import Chem

def process_second_column(ligandID_x):

    try:
        ligandID_x = int(ligandID_x)
    except:
        return ""
    
    subfolder = f"Z02_All_Ligands_Molfiles_{ ( (ligandID_x - 1) // 50000 ) * 50000 + 1 }_{ ( 1 + (ligandID_x - 1) // 50000 ) * 50000 }"
    
    try:
        mol    = Chem.MolFromMolFile(f'../Z02_All_Ligands_Molfiles/{subfolder}/{ligandID_x}.mol')
        smiles = Chem.MolToSmiles(mol)
    except:
        smiles = ""

    return smiles

with open(input_file, 'r', newline='', encoding='utf-8') as infile:
    reader = csv.reader(infile, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL, skipinitialspace=True)

    with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.writer(outfile, delimiter=';', quotechar='"', quoting=csv.QUOTE_MINIMAL)

        for row in reader:
            new_value   = process_second_column(row[1])
            updated_row = row + [new_value]
            writer.writerow(updated_row)

print(f"New column added to the {output_file} based on processing the second column.")








