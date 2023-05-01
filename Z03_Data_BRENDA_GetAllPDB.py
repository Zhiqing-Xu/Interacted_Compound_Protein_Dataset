#!/usr/bin/env python
# coding: utf-8
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
# The following code ensures the code work properly in 
# MS VS, MS VS CODE and jupyter notebook on both Linux and Windows.
#--------------------------------------------------#
import os 
import sys
import os.path
from sys import platform
from pathlib import Path
#--------------------------------------------------#
if __name__ == "__main__":
    print("="*80)
    if os.name == 'nt' or platform == 'win32':
        print("Running on Windows")
        if 'ptvsd' in sys.modules:
            print("Running in Visual Studio")
#--------------------------------------------------#
    if os.name != 'nt' and platform != 'win32':
        print("Not Running on Windows")
#--------------------------------------------------#
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
#--------------------------------------------------#

###################################################################################################################
###################################################################################################################
# Imports
#--------------------------------------------------#
import re
import copy
import pickle
import argparse
import numpy as np
import pandas as pd
#--------------------------------------------------#
from urllib.parse import quote
from urllib.request import urlopen
from unidecode import unidecode
###################################################################################################################
###################################################################################################################
# Basic Functions
# Print the DataFrame obtained.
def beautiful_print(df):
    # Print the dataset in a well-organized format.
    with pd.option_context('display.max_rows'      , 20, 
                           'display.min_rows'      , 20, 
                           'display.max_columns'   , 7, 
                           #"display.max_colwidth" , None,
                           "display.width"         , None,
                           "expand_frame_repr"     , True,
                           "max_seq_items"         , None,):  # more options can be specified
        # Once the display.max_rows is exceeded, 
        # the display.min_rows options determines 
        # how many rows are shown in the truncated repr.
        print(df)
    return 

#--------------------------------------------------#
def get_longest_smiles_from_aggregates(one_smiles_string):
    return max(one_smiles_string.split("."), key = len)

#--------------------------------------------------#
def isfloat(num): 
    # Return True if the input can be converted to a float or is a float already.
    try:
        float(num)
        return True
    except ValueError:
        return False

###################################################################################################################
###################################################################################################################
# For Accessing BRENDA.
import hashlib
from zeep import Client

wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
password = hashlib.sha256("2-Oxoglutarate".encode("utf-8")).hexdigest()
client = Client(wsdl)















#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#                     `7MMF'`7MN.   `7MF'`7MM"""Mq. `7MMF'   `7MF'MMP""MM""YMM  .M"""bgd                                                               #
#                       MM    MMN.    M    MM   `MM.  MM       M  P'   MM   `7 ,MI    "Y                                                               #
#    ,pP""Yq.           MM    M YMb   M    MM   ,M9   MM       M       MM      `MMb.                                                                   #
#   6W'    `Wb          MM    M  `MN. M    MMmmdM9    MM       M       MM        `YMMNq.                                                               #
#   8M      M8          MM    M   `MM.M    MM         MM       M       MM      .     `MM                                                               #
#   YA.    ,A9 ,,       MM    M     YMM    MM         YM.     ,M       MM      Mb     dM                                                               #
#    `Ybmmd9'  db     .JMML..JML.    YM  .JMML.        `bmmmmd"'     .JMML.    P"Ybmmd"                                                                #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 

def Step_0_get_raw_df_0(dataset_idx = 1): # Shall be 1, 2 and 3 only. 

    # Input Arguments

    data_folder_name    = ["Ki_BRENDA"        , "KM_BRENDA"        , "kcat_BRENDA"        , "kcat_KM_BRENDA"        ][dataset_idx]
    data_file_name      = ["brenda_Ki_raw.csv", "brenda_KM_raw.csv", "brenda_kcat_raw.csv", "brenda_kcat_KM_raw.csv"][dataset_idx]
    output_file_name    = ["Ki_BRENDA.csv"    , "KM_BRENDA.csv"    , "kcat_BRENDA.csv"    , "kcat_KM_BRENDA.csv"    ][dataset_idx]
    data_name           = ["Ki"               , "Km"               , "kcat"               , "kcat_KM"               ][dataset_idx]
    data_name_0         = ["max_Ki"           , "max_Km"           , "max_kcat"           , "max_kcat_KM"           ][dataset_idx]

    data_folder         = Path("Z00_BRENDA_Kinetics_Raw/" + data_folder_name)
    data_file           = data_file_name

    saving_file_suffix  = "_extended.csv"

    output_folder       = Path("Z04_DataPreprocessing_Savings/")
    output_file         = output_file_name

    output_1_suffix     = '_core.csv' # 
    output_2_suffix     = '_scrn.csv' # 
    output_3_suffix     = '_fine.csv' # 



    ################################################################################################################### 
    #         `7MMM.     ,MMF'      db      `7MMF'`7MN.   `7MF'     `7MM"""YMM `7MMF'`7MMF'      `7MM"""YMM           #
    #           MMMb    dPMM       ;MM:       MM    MMN.    M         MM    `7   MM    MM          MM    `7           #
    #           M YM   ,M MM      ,V^MM.      MM    M YMb   M         MM   d     MM    MM          MM   d             #
    #           M  Mb  M' MM     ,M  `MM      MM    M  `MN. M         MM""MM     MM    MM          MMmmMM             #
    #           M  YM.P'  MM     AbmmmqMA     MM    M   `MM.M         MM   Y     MM    MM      ,   MM   Y  ,          #
    #           M  `YM'   MM    A'     VML    MM    M     YMM         MM         MM    MM     ,M   MM     ,M          #
    #         .JML. `'  .JMML..AMA.   .AMMA..JMML..JML.    YM       .JMML.     .JMML..JMMmmmmMMM .JMMmmmmMMM          #
    ###################################################################################################################

    print("\n\n" + "#"*100 + "\n" + "#"*100 +"\n>>> Part 0: Read the raw dataset and perform rough cleaning... ")

    # Read the main data file.
    raw_df_0 = pd.read_csv(filepath_or_buffer   =   data_folder / data_file , 
                           on_bad_lines         =   'skip'                  , 
                           index_col            =   None                    , 
                           #names               =   ["", "", ]              , 
                           header               =   None                    , 
                           sep                  =   '	'                   , 
                           encoding             =   'cp1252'                ,)

    raw_df_0.rename(columns = {0 : 'EC'          , 
                            1 : data_name     , 
                            2 : data_name_0   , 
                            3 : 'cmpd'        , 
                            4 : 'condition'   , 
                            5 : 'organism'    , 
                            6 : 'uniprot_id'  , 
                            7 : 'publication' , 
                            8 : 'other'       , 
                            }, 
                    inplace = True)

    # The cmpd column is 

    print("\n\n"+"-"*90+"\n# Step 0.0 Input: Raw data read from the csv file, raw_df_0: ")
    beautiful_print(raw_df_0)
    print("len(raw_df_0): ", len(raw_df_0))


    #====================================================================================================#
    # Rough Prescreening.
    raw_df_0_1 = copy.deepcopy(raw_df_0) # Make a copy as usual.


    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Rough Prescreening #1:
    raw_df_0_1 = raw_df_0_1[ raw_df_0_1[data_name].str.contains("additional information") == False ] #1
    raw_df_0_1 = raw_df_0_1[ raw_df_0_1[data_name].str.strip() != "-" ]                              #1
    #print(len(raw_df_0_1)) # [out]: 162216.


    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Rough Prescreening #2:
    raw_df_0_2 = raw_df_0_1[ raw_df_0_1["cmpd"].str.contains("additional information") == False]     #2
    #print(len(raw_df_0_2)) # [out]: 162208.

    # Some comments on Rough Prescreening #2
    """
    In Rough Prescreening #2, there is NOT a row of 'additional information',
    yet some rows of the dataframe are still removed.

    This is because there are rows with too few fields in the dataset.
    The function 'read_csv()' read those rows and the value of those fields as NaN.

    The #2 filtering remove those rows as the value of those fields are requested.
    # [in]: print(raw_df_0[["cmpd",]].iloc[38737]) [out]: NaN
    """


    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Remove those rows with a max_[data_name] value (optional).
    raw_df_0_3  = raw_df_0_2[ raw_df_0_2[data_name_0].str.contains('-') == True ] #3
    #print("len(raw_df_0): ", len(raw_df_0_3))

    raw_df_0 = raw_df_0_2


    #====================================================================================================#
    # Remove useless columnms.
    try:
        raw_df_0 = raw_df_0.drop(columns=[data_name_0, "condition", "publication", "other"])
    except:
        raw_df_0 = raw_df_0.drop(columns=[data_name_0, "condition", "publication", ])


    #====================================================================================================#
    # Change the order of the columns.
    raw_df_0 = raw_df_0.iloc[:, [-2, 0, -1, 2, 1]]


    #====================================================================================================#
    # Simple Preprocessing.
    raw_df_0['uniprot_id'] = raw_df_0['uniprot_id'].str.replace('and', ',')
    raw_df_0['uniprot_id'] = raw_df_0['uniprot_id'].str.replace('AND', ',')
    raw_df_0['uniprot_id'] = raw_df_0['uniprot_id'].str.replace(' ', '')

    #====================================================================================================#
    # Part #0 Checkpoint. 
    raw_df_0.reset_index(inplace=True)
    raw_df_0.drop(columns = ["index", ], inplace = True)
    print("\n\n"+"-"*90+"\n# Step 0.1 Prescreen and simple clean, raw_df_0: ")
    beautiful_print(raw_df_0)
    print("len(raw_df_0): ", len(raw_df_0))

    return raw_df_0









#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#                `7MMF'   `7MF'`7MN.   `7MF'`7MMF'`7MM"""Mq. `7MM"""Mq.    .g8""8q.   MMP""MM""YMM     `7MMF'`7MM"""Yb.                                #
#    __,           MM       M    MMN.    M    MM    MM   `MM.  MM   `MM. .dP'    `YM. P'   MM   `7       MM    MM    `Yb.                              #
#   `7MM           MM       M    M YMb   M    MM    MM   ,M9   MM   ,M9  dM'      `MM      MM            MM    MM     `Mb ,pP"Ybd                      #
#     MM           MM       M    M  `MN. M    MM    MMmmdM9    MMmmdM9   MM        MM      MM            MM    MM      MM 8I   `"                      #
#     MM           MM       M    M   `MM.M    MM    MM         MM  YM.   MM.      ,MP      MM            MM    MM     ,MP `YMMMa.                      #
#     MM  ,,       YM.     ,M    M     YMM    MM    MM         MM   `Mb. `Mb.    ,dP'      MM            MM    MM    ,dP' L.   I8                      #
#   .JMML.db        `bmmmmd"'  .JML.    YM  .JMML..JMML.     .JMML. .JMM.  `"bmmd"'      .JMML.        .JMML..JMMmmmdP'   M9mmmP'                      #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 

def Step_1_get_data_wi_unip_df(raw_df_0 = None): # Shall be 1, 2 and 3 only. 

    print("\n\n" + "#"*100 + "\n" + "#"*100 +"\n>>> Part 1: Split the dataset (raw_df_0) into two df's, one with uniprot_id, one without... ")
    print("    ( raw_df_0 -> data_wi_unip_df + data_wo_unip_df ) ")

    #====================================================================================================#
    # Separating data with vs without Uniprot IDs
    # Check UniProt IDs for each record in the dataset and split the data accordingly.
    data_wi_unip_df = data_wo_unip_df = copy.deepcopy(raw_df_0)
    data_wi_unip_df = data_wi_unip_df[ data_wi_unip_df['uniprot_id'].str.contains('-') == False ]
    data_wo_unip_df = data_wo_unip_df[ data_wo_unip_df['uniprot_id'].str.contains('-') == True  ]

    #====================================================================================================#
    # Exploding data points with multiple Uniprot IDs listed
    print("\n\n"+"-"*90+"\n# Step 1.0 data_wi_unip_df, data_wo_unip_df: ")
    print("len(data_wi_unip_df) before/after exploding with multiple Uniprot IDs: ", len(data_wi_unip_df), end = ", ")
    data_wi_unip_df['uniprot_id'] = data_wi_unip_df['uniprot_id'].str.split(',')
    data_wi_unip_df = data_wi_unip_df.explode('uniprot_id')
    print(len(data_wi_unip_df))

    #====================================================================================================#
    # Create the Sequence column.
    data_wi_unip_df['Sequence'] = np.nan
    data_wo_unip_df['Sequence'] = np.nan

    #====================================================================================================#
    # Part #1 Checkpoint. 
    print("len(data_wo_unip_df): ", len(data_wo_unip_df))
    print("len(raw_df_0): ", len(raw_df_0))
    print("Confirmed that: `` len(raw_df_0) == len(data_wo_unip_df) + len(data_wi_unip_df) `` (before explosion)")

    data_wi_unip_df.reset_index(inplace=True)
    data_wi_unip_df.drop(columns = ["index", ], inplace = True)
    print("\n\n"+"-"*90+"\n# Step 1.1 Dataset with Uniprot IDs: data_wi_unip_df: ")
    beautiful_print(data_wi_unip_df)
    print("len(data_wi_unip_df): ", len(data_wi_unip_df))


    data_wo_unip_df.reset_index(inplace=True)
    data_wo_unip_df.drop(columns = ["index", ], inplace = True)
    print("\n\n"+"-"*90+"\n# Step 1.2 Dataset without Uniprot IDs: data_wo_unip_df: ")
    beautiful_print(data_wo_unip_df)
    print("len(data_wo_unip_df): ", len(data_wo_unip_df))

    return data_wi_unip_df









#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#                                     `7MM"""YMM  `YMM'   `MP'`7MM"""Mq.   db     `7MN.   `7MF'`7MM"""Yb.                                              #
#                                       MM    `7    VMb.  ,P    MM   `MM. ;MM:      MMN.    M    MM    `Yb.                                            #
#                                       MM   d       `MM.M'     MM   ,M9 ,V^MM.     M YMb   M    MM     `Mb                                            #
#                                       MMmmMM         MMb      MMmmdM9 ,M  `MM     M  `MN. M    MM      MM                                            #
#                                       MM   Y  ,    ,M'`Mb.    MM      AbmmmqMA    M   `MM.M    MM     ,MP                                            #
#                                       MM     ,M   ,P   `MM.   MM     A'     VML   M     YMM    MM    ,dP'                                            #
#                                     .JMMmmmmMMM .MM:.  .:MMa.JMML. .AMA.   .AMMA.JML.    YM  .JMMmmmdP'                                              #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


data_folder         = Path("Z00_BRENDA_Kinetics_Raw/")
data_file           = "Z00_Query_Result.p"


raw_df_1 = Step_0_get_raw_df_0(dataset_idx = 1)
raw_df_2 = Step_0_get_raw_df_0(dataset_idx = 2)
raw_df_3 = Step_0_get_raw_df_0(dataset_idx = 3)

data_wi_unip_df_1 = Step_1_get_data_wi_unip_df(raw_df_1)
data_wi_unip_df_2 = Step_1_get_data_wi_unip_df(raw_df_2)
data_wi_unip_df_3 = Step_1_get_data_wi_unip_df(raw_df_3)

uniprotID_list_1 = data_wi_unip_df_1["uniprot_id"].tolist()
uniprotID_list_2 = data_wi_unip_df_2["uniprot_id"].tolist()
uniprotID_list_3 = data_wi_unip_df_3["uniprot_id"].tolist()

uniprotID_list_all = list(set(uniprotID_list_1 + uniprotID_list_2 + uniprotID_list_3))

print("Number of UniProt ID's in KM       dataset: ", len(list(set(uniprotID_list_1))))
print("Number of UniProt ID's in Kcat     dataset: ", len(list(set(uniprotID_list_2))))
print("Number of UniProt ID's in Kcat/KM  dataset: ", len(list(set(uniprotID_list_3))))
print("Number of UniProt ID's in Combined dataset: ", len(uniprotID_list_all))



import requests
from Bio import SeqIO
from io import StringIO


def get_pdb_id(uniprot_id):
    url = f'https://www.uniprot.org/uniprot/{uniprot_id}.fasta'
    response = requests.get(url)
    fasta_str = response.text
    
    try:
        fasta = next(SeqIO.parse(StringIO(fasta_str), "fasta"))
        print(fasta)
        # Extract PDB IDs from the FASTA header
        pdb_ids = [x.split(":")[1] for x in fasta.description.split() if "PDB" in x]
        print(pdb_ids)
        return pdb_ids
    except:
        return []
    

def download_pdb(pdb_id, file_path):
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    response = requests.get(url)

    with open(file_path, 'w') as f:
        f.write(response.text)




def download_alphafold_pdb(uniprot_id, file_path):
    url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb'
    response = requests.get(url)

    with open(file_path, 'w') as f:
        f.write(response.text)

for idx, uniprot_id in enumerate(uniprotID_list_all):
    print(f'Downloading {idx} out of 9965...')
    download_alphafold_pdb(uniprot_id, f'Z03_PDB_files/alphafold/{uniprot_id}.pdb')






#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#


#       M              M              M              M              M               M              M              M              M              M      #
#       M              M              M              M              M               M              M              M              M              M      #
#       M              M              M              M              M               M              M              M              M              M      #
#   `7M'M`MF'      `7M'M`MF'      `7M'M`MF'      `7M'M`MF'      `7M'M`MF'       `7M'M`MF'      `7M'M`MF'      `7M'M`MF'      `7M'M`MF'      `7M'M`MF'  #
#     VAMAV          VAMAV          VAMAV          VAMAV          VAMAV           VAMAV          VAMAV          VAMAV          VAMAV          VAMAV    #
#      VVV            VVV            VVV            VVV            VVV             VVV            VVV            VVV            VVV            VVV     #
#       V              V              V              V              V               V              V              V              V              V      #

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
###################################################################################################################
###################################################################################################################
#====================================================================================================#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#--------------------------------------------------#
#------------------------------

#                                                                                                                                                          
#      `MM.              `MM.             `MM.             `MM.             `MM.             `MM.             `MM.             `MM.             `MM.       
#        `Mb.              `Mb.             `Mb.             `Mb.             `Mb.             `Mb.             `Mb.             `Mb.             `Mb.     
# MMMMMMMMMMMMD     MMMMMMMMMMMMD    MMMMMMMMMMMMD    MMMMMMMMMMMMD    MMMMMMMMMMMMD    MMMMMMMMMMMMD    MMMMMMMMMMMMD    MMMMMMMMMMMMD    MMMMMMMMMMMMD   
#         ,M'               ,M'              ,M'              ,M'              ,M'              ,M'              ,M'              ,M'              ,M'     
#       .M'               .M'              .M'              .M'              .M'              .M'              .M'              .M'              .M'       
#                                                                                                                                                          

#------------------------------
#--------------------------------------------------#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#====================================================================================================#
###################################################################################################################
###################################################################################################################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

#       A              A              A              A              A               A              A              A              A              A      #
#      MMM            MMM            MMM            MMM            MMM             MMM            MMM            MMM            MMM            MMM     #
#     MMMMM          MMMMM          MMMMM          MMMMM          MMMMM           MMMMM          MMMMM          MMMMM          MMMMM          MMMMM    #
#   ,MA:M:AM.      ,MA:M:AM.      ,MA:M:AM.      ,MA:M:AM.      ,MA:M:AM.       ,MA:M:AM.      ,MA:M:AM.      ,MA:M:AM.      ,MA:M:AM.      ,MA:M:AM.  #
#       M              M              M              M              M               M              M              M              M              M      #
#       M              M              M              M              M               M              M              M              M              M      #
#       M              M              M              M              M               M              M              M              M              M      #




