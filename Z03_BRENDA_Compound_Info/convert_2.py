#!/usr/bin/env python
# coding: utf-8
#====================================================================================================#
# The following code ensures the code work properly in 
# MS VS, MS VS CODE and jupyter notebook on both Linux and Windows.
import os 
import sys
from os import path
from sys import platform
from pathlib import Path

if __name__ == "__main__":
    print("="*80)
    if os.name == 'nt' or platform == 'win32':
        print("Running on Windows")
        if 'ptvsd' in sys.modules:
            print("Running in Visual Studio")

    if os.name != 'nt' and platform != 'win32':
        print("Not Running on Windows")

    if "__file__" in globals().keys():
        print('CurrentDir: ', os.getcwd())
        try:
            os.chdir(os.path.dirname(__file__))
        except:
            print("Problems with navigating to the file dir.")
        print('CurrentDir: ', os.getcwd())
    else:
        print("Running in python jupyter notebook.")
        try:
            if not 'workbookDir' in globals():
                workbookDir = os.getcwd()
                print('workbookDir: ' + workbookDir)
                os.chdir(workbookDir)
        except:
            print("Problems with navigating to the workbook dir.")
#====================================================================================================#
# Imports
import csv
from rdkit import Chem

#====================================================================================================#
# Globals
input_file  = 'Z03B_CompoundName_LigandID_SMILES.csv'
output_file = 'Z03B_CompoundName_LigandID_SMILES_updated.csv'






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

        for i, row in enumerate(reader):
            print(f"Processing row {i+1} ...")
            new_value    = process_second_column(row[1])    # compound SMILES. 
            updated_row  = row[0:2] + [new_value]           # 
            writer.writerow(updated_row)                    # 





print(f"New column added to the {output_file} based on processing the second column.")





































