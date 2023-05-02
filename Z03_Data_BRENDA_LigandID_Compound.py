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
    print("\n" + "#"*100 + "\n" + "#"*100)
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
import time
import copy
import pickle
import argparse
import numpy as np
import pandas as pd
#--------------------------------------------------#
import requests
import xmltodict
#--------------------------------------------------#
from timeit import timeit
#--------------------------------------------------#
import urllib
import xml.etree.ElementTree as ET
from urllib.request import urlopen


###################################################################################################################
###################################################################################################################
# Basic Functions

def beautiful_print(df): # DataFrame Printing.
    # Print the dataset in a well-organized format.
    with pd.option_context('display.max_rows'      , 20  , 
                           'display.min_rows'      , 20  , 
                           'display.max_columns'   , 7   , 
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
#                    `7MMF'`7MN.   `7MF'`7MM"""Mq. `7MMF'   `7MF'MMP""MM""YMM  .M"""bgd                                                                #
#                      MM    MMN.    M    MM   `MM.  MM       M  P'   MM   `7 ,MI    "Y                                                                #
#   ,pP""Yq.           MM    M YMb   M    MM   ,M9   MM       M       MM      `MMb.                                                                    #
#  6W'    `Wb          MM    M  `MN. M    MMmmdM9    MM       M       MM        `YMMNq.                                                                #
#  8M      M8          MM    M   `MM.M    MM         MM       M       MM      .     `MM                                                                #
#  YA.    ,A9 ,,       MM    M     YMM    MM         YM.     ,M       MM      Mb     dM                                                                #
#   `Ybmmd9'  db     .JMML..JML.    YM  .JMML.        `bmmmmd"'     .JMML.    P"Ybmmd"                                                                 #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 

# Input Arguments
data_folder_name    = ["Ki_BRENDA"        , "KM_BRENDA"        , "kcat_BRENDA"        , "kcat_KM_BRENDA"        ][2]
data_file_name      = ["brenda_Ki_raw.csv", "brenda_KM_raw.csv", "brenda_kcat_raw.csv", "brenda_kcat_KM_raw.csv"][2]
output_file_name    = ["Ki_BRENDA.csv"    , "KM_BRENDA.csv"    , "kcat_BRENDA.csv"    , "kcat_KM_BRENDA.csv"    ][2]
data_name           = ["Ki"               , "Km"               , "kcat"               , "kcat_KM"               ][2]
data_name_0         = ["max_Ki"           , "max_Km"           , "max_kcat"           , "max_kcat_KM"           ][2]

data_folder         = Path("Z00_BRENDA_Kinetics_Raw/" + data_folder_name)
data_file           = data_file_name

saving_file_suffix  = "_extended.csv"

output_folder       = Path("Z03_BRENDA_Compound_Info/")
output_file         = output_file_name

output_1_suffix     = '_core.csv' # 
output_2_suffix     = '_scrn.csv' # 
output_3_suffix     = '_fine.csv' # 


# Set Threshold for Screening.
if data_name in ["Ki", "Km", "kcat", ]:
    val_SD_threshold = 1
    
if data_name in ["kcat_KM", ]:
    val_SD_threshold = 10


################################################################################################################### 
#         `7MMM.     ,MMF'      db      `7MMF'`7MN.   `7MF'     `7MM"""YMM `7MMF'`7MMF'      `7MM"""YMM           #
#           MMMb    dPMM       ;MM:       MM    MMN.    M         MM    `7   MM    MM          MM    `7           #
#           M YM   ,M MM      ,V^MM.      MM    M YMb   M         MM   d     MM    MM          MM   d             #
#           M  Mb  M' MM     ,M  `MM      MM    M  `MN. M         MM""MM     MM    MM          MMmmMM             #
#           M  YM.P'  MM     AbmmmqMA     MM    M   `MM.M         MM   Y     MM    MM      ,   MM   Y  ,          #
#           M  `YM'   MM    A'     VML    MM    M     YMM         MM         MM    MM     ,M   MM     ,M          #
#         .JML. `'  .JMML..AMA.   .AMMA..JMML..JML.    YM       .JMML.     .JMML..JMMmmmmMMM .JMMmmmmMMM          #
###################################################################################################################

print("\n\n\n\n" + "#"*100 + "\n" + "#"*100 +"\n>>> Part 0: Read the raw dataset and perform rough cleaning... ")

# Read the main data file.
raw_df_0 = pd.read_csv(filepath_or_buffer   =   data_folder / data_file, 
                       on_bad_lines         =   'skip', 
                       index_col            =   None, 
                       #names               =   ["", "", ], 
                       header               =   None, 
                       sep                  =   '	', 
                       encoding             =   'cp1252')

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











################################################################################################################### 
#                   `7MM"""YMM  `YMM'   `MP'`7MM"""Mq.   db     `7MN.   `7MF'`7MM"""Yb.                           #
#                     MM    `7    VMb.  ,P    MM   `MM. ;MM:      MMN.    M    MM    `Yb.                         #
#                     MM   d       `MM.M'     MM   ,M9 ,V^MM.     M YMb   M    MM     `Mb                         #
#                     MMmmMM         MMb      MMmmdM9 ,M  `MM     M  `MN. M    MM      MM                         #
#                     MM   Y  ,    ,M'`Mb.    MM      AbmmmqMA    M   `MM.M    MM     ,MP                         #
#                     MM     ,M   ,P   `MM.   MM     A'     VML   M     YMM    MM    ,dP'                         #
#                   .JMMmmmmMMM .MM:.  .:MMa.JMML. .AMA.   .AMMA.JML.    YM  .JMMmmmdP'                           #
###################################################################################################################

raw_df_0_extended = copy.deepcopy(raw_df_0)

from Z01_Data_BRENDA_GetReaction import *

# Expand the dataset through searching BRENDA for ``reaction partners`` and ``ligand IDs``.
if os.path.exists(data_folder / "../../Z01_BRENDA_QueryResults/Z01_Query_Result.p"):
    Z01_Query_Result_df = pd.read_pickle  (data_folder / "../../Z01_BRENDA_QueryResults/Z01_Query_Result.p" )
    Unresponded_Query   = pickle.load(open(data_folder / "../../Z01_BRENDA_QueryResults/Z01_Unresponded_Query.p", 'rb'))

else:
    print("Saving File NOT Found. You have to run Z01_Data_BRENDA_GetReaction.py")
    Z01_Get_Reaction_Main()


print("\n\n"+"-"*90+"\n# Step 0.2 Reaction Info, Z01_Query_Result_df: ")
beautiful_print(Z01_Query_Result_df[["organism","substrate","EC","reaction"]])
print("len(Z01_Query_Result_df): ", len(Z01_Query_Result_df))
count = Z01_Query_Result_df['reaction'].apply(lambda x: x != []).sum()
print("Number of non-empty query result: ", count)

Z01_Query_Result_df = Z01_Query_Result_df.rename( columns = {"substrate": "cmpd",} )












################################################################################################################### 
#                                                                                                                 #
#                          mm        .M"""bgd `7MMM.     ,MMF'`7MMF'`7MMF'      `7MM"""YMM   .M"""bgd             #
#                          MM       ,MI    "Y   MMMb    dPMM    MM    MM          MM    `7  ,MI    "Y             #
#        .P"Ybmmm .gP"Ya mmMMmm     `MMb.       M YM   ,M MM    MM    MM          MM   d    `MMb.                 #
#       :MI  I8  ,M'   Yb  MM         `YMMNq.   M  Mb  M' MM    MM    MM          MMmmMM      `YMMNq.             #
#        WmmmP"  8M""""""  MM       .     `MM   M  YM.P'  MM    MM    MM      ,   MM   Y  , .     `MM             #
#       8M       YM.    ,  MM       Mb     dM   M  `YM'   MM    MM    MM     ,M   MM     ,M Mb     dM             #
#        YMMMMMb  `Mbmmd'  `Mbmo    P"Ybmmd"  .JML. `'  .JMML..JMML..JMMmmmmMMM .JMMmmmmMMM P"Ybmmd"              #
#       6'     dP                                                                                                 #
#       Ybmmmd'                                                                                                   #
###################################################################################################################
# Based on the ligand ID's and compound name obtained from the query results, we are able to get a dictionary.
# The dictionary allows mapping from `compound name` -> `ligand ID` -> `mol file` -> `SMILES`.
# Here, generate a three-column csv file and a compound_name -> SMILES dictionary (saved as .p file).


Z01_Query_Result_df_subs_ligd = Z01_Query_Result_df[["cmpd", "ligandID"]].dropna()
Z01_Query_Result_df_subs_ligd = Z01_Query_Result_df_subs_ligd[Z01_Query_Result_df_subs_ligd['ligandID'] != 0]
Z01_Query_Result_df_subs_ligd = Z01_Query_Result_df_subs_ligd.drop_duplicates()

# Get mapping from Compound Name to Ligand IDs, ~22212 mappings exists.
Z01_Query_Result_df_subs_ligd.reset_index(inplace=True)
Z01_Query_Result_df_subs_ligd.drop(columns = ["index", ], inplace = True)
print("\n\n"+"-"*90+"\n# Step 0.3 compound_name -> SMILES, Z01_Query_Result_df_subs_ligd: ")
beautiful_print(Z01_Query_Result_df_subs_ligd)
print("len(Z01_Query_Result_df_subs_ligd): ", len(Z01_Query_Result_df_subs_ligd))


# Get All SMILES.
if os.path.exists(output_folder / "Z03_CompoundName_LigandID_SMILES.p"):
    Z03_CompoundName_LigandID_SMILES = pd.read_pickle(output_folder / "Z03_CompoundName_LigandID_SMILES.p")

else:
    all_ligandIDs = Z01_Query_Result_df_subs_ligd["ligandID"].tolist()
    all_SMILES    = []

    from rdkit import Chem
    for idx, ligandID_x in enumerate(all_ligandIDs):

        subfolder = f"Z02_All_Ligands_Molfiles_{ ( (ligandID_x - 1) // 50000 ) * 50000 + 1 }_{ ( 1 + (ligandID_x - 1) // 50000 ) * 50000 }"
        try:
            mol = Chem.MolFromMolFile(f'Z02_All_Ligands_Molfiles/{subfolder}/{ligandID_x}.mol')
            smiles = Chem.MolToSmiles(mol)
        except:
            smiles = ""
        print(idx, ligandID_x, smiles)
        all_SMILES.append(smiles)

    Z01_Query_Result_df_subs_ligd["smiles"] = all_SMILES
    Z03_CompoundName_LigandID_SMILES = copy.deepcopy(Z01_Query_Result_df_subs_ligd)
    Z03_CompoundName_LigandID_SMILES.to_pickle(output_folder / "Z03_CompoundName_LigandID_SMILES.p")


print("\n\n"+"-"*90+"\n# Step 0.4 compound_name -> SMILES, Z03_CompoundName_LigandID_SMILES: ")
beautiful_print(Z03_CompoundName_LigandID_SMILES)
print("len(Z01_CompoundName_LigandID_SMILES): ", len(Z03_CompoundName_LigandID_SMILES))

Z03_CompoundName_LigandID_SMILES.rename(columns = {"cmpd"     : "CMPD"           , 
                                                   "ligandID" : "ligandID"       , 
                                                   "smiles"   : "CMPD_SMILES"    , 
                                                   }, 
                                        inplace = True)

Z03_CompoundName_LigandID_SMILES.to_csv(output_folder / "Z03_CompoundName_LigandID_SMILES.csv", sep=';', index = False)











################################################################################################################### 
#          `7MM"""Mq. `7MM"""YMM        db       .g8"""bgd MMP""MM""YMM `7MMF' .g8""8q. `7MN.   `7MF'             #
#            MM   `MM.  MM    `7       ;MM:    .dP'     `M P'   MM   `7   MM .dP'    `YM. MMN.    M               #
#            MM   ,M9   MM   d        ,V^MM.   dM'       `      MM        MM dM'      `MM M YMb   M               #
#            MMmmdM9    MMmmMM       ,M  `MM   MM               MM        MM MM        MM M  `MN. M               #
#            MM  YM.    MM   Y  ,    AbmmmqMA  MM.              MM        MM MM.      ,MP M   `MM.M               #
#            MM   `Mb.  MM     ,M   A'     VML `Mb.     ,'      MM        MM `Mb.    ,dP' M     YMM               #
#          .JMML. .JMM.JMMmmmmMMM .AMA.   .AMMA. `"bmmmd'     .JMML.    .JMML. `"bmmd"' .JML.    YM               #
###################################################################################################################
# Get All Compound Names from the reaction strings returned by the queries. 


Z01_All_Compounds = copy.deepcopy(Z01_Query_Result_df.explode('reaction'))

raw_df_0.reset_index(inplace=True)
raw_df_0.drop(columns = ["index", ], inplace = True)
print("\n\n"+"-"*90+"\n# Step 0.5 Exploded Z01_Query_Result_df, Z01_All_Compounds: ")
beautiful_print(Z01_All_Compounds)
print("len(Z01_CompoundName_LigandID_SMILES): ", len(Z01_All_Compounds))



subs_list = Z01_All_Compounds["cmpd"     ].tolist()
rctn_list = Z01_All_Compounds["reaction" ].tolist()

all_cmpds_list = []

for one_rctn_info in rctn_list:
    if one_rctn_info is not None:
        if not (isinstance(one_rctn_info, (int, float)) and np.isnan(one_rctn_info)):
            #print(one_rctn_info)
            for one_group in one_rctn_info.split(" = "):
                for one_cmpd in one_group.split(" + "):
                    all_cmpds_list.append(one_cmpd)


all_cmpds_list += subs_list

all_cmpds_list = list(set(all_cmpds_list))

print("\n\n"+"-"*90+"\n# Step 0.6 Get all compound names : ")
# print(all_cmpds_list)
print("Number of all compound names found : ", len(all_cmpds_list) )

































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




