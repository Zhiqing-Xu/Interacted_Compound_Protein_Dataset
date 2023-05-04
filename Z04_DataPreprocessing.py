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

output_folder       = Path("Z04_DataPreprocessing_Savings/")
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

#====================================================================================================#
# Use the only unresponded query as an example here, show how an example of querying BRENDA
# Unresponded_Query = [('Pyrococcus abyssi', 'S-adenosyl-L-methionine', '2.1.1.B114')]
subs, orgm, ecnm = Unresponded_Query[0]

parameters = ("palladium.xu@gmail.com"              , 
              password                              , 
              "ecNumber*"           + ecnm          , 
              "substrate*"          + subs          , 
              "reactionPartners*"                   , 
              "organism*"           + orgm          , 
              "ligandStructureId*"                  ,
              )
try:
    query_result = client.service.getSubstrate(*parameters)
except:
    query_result = []
    print("bad one found! ")

print("Previously unresponded query retried: ", query_result) # Empty Result shows NO reaction info.

#====================================================================================================#
# Get a dataframe that works as a organism + EC + Substrate -> reaction dictionary.
OG_EC_CP__RXN_grouped_df = Z01_Query_Result_df.groupby(['organism', 'EC', 'cmpd'])['reaction'].apply(list)
OG_EC_CP__RXN_grouped_df = OG_EC_CP__RXN_grouped_df.to_frame().reset_index(inplace = False)

OG_EC_CP__LGD_grouped_df = Z01_Query_Result_df.groupby(['organism', 'EC', 'cmpd'])['ligandID'].apply(list)
OG_EC_CP__LGD_grouped_df = OG_EC_CP__LGD_grouped_df.to_frame().reset_index(inplace = False)

# Merging the df with the sequence df to get a combined df.
raw_df_0_extended = raw_df_0_extended.merge(OG_EC_CP__RXN_grouped_df, on = ['organism', 'EC', 'cmpd'], how = 'left')
raw_df_0_extended = raw_df_0_extended.explode('reaction') # Exploding the df.

raw_df_0_extended = raw_df_0_extended.merge(OG_EC_CP__LGD_grouped_df, on = ['organism', 'EC', 'cmpd'], how = 'left')
raw_df_0_extended = raw_df_0_extended.explode('ligandID') # Exploding the df.



#====================================================================================================#
# Part #0.3 Checkpoint. 
print("\n\n"+"="*90+"\n#0.3 Extended raw_df_0, raw_df_0_extended: ")
beautiful_print(raw_df_0_extended)
print("len(raw_df_0_extended): ", len(raw_df_0_extended))

count = raw_df_0_extended['reaction'].apply(lambda x: x != []).sum()
print("Number of non-empty reaction & ligandID: ", count)

raw_df_0 = copy.deepcopy(raw_df_0_extended)














#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#                `7MMF'   `7MF'`7MN.   `7MF'`7MMF'`7MM"""Mq. `7MM"""Mq.    .g8""8q.   MMP""MM""YMM     `7MMF'`7MM"""Yb.                                #
#    __,           MM       M    MMN.    M    MM    MM   `MM.  MM   `MM. .dP'    `YM. P'   MM   `7       MM    MM    `Yb.                              #
#   `7MM           MM       M    M YMb   M    MM    MM   ,M9   MM   ,M9  dM'      `MM      MM            MM    MM     `Mb ,pP"Ybd                      #
#     MM           MM       M    M  `MN. M    MM    MMmmdM9    MMmmdM9   MM        MM      MM            MM    MM      MM 8I   `"                      #
#     MM           MM       M    M   `MM.M    MM    MM         MM  YM.   MM.      ,MP      MM            MM    MM     ,MP `YMMMa.                      #
#     MM  ,,       YM.     ,M    M     YMM    MM    MM         MM   `Mb. `Mb.    ,dP'      MM            MM    MM    ,dP' L.   I8                      #
#   .JMML.db        `bmmmmd"'  .JML.    YM  .JMML..JMML.     .JMML. .JMM.  `"bmmd"'      .JMML.        .JMML..JMMmmmdP'   M9mmmP'                      #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 

print("\n\n\n\n" + "#"*100 + "\n" + "#"*100 +"\n>>> Part 1: Split the dataset (raw_df_0) into two df's, one with uniprot_id, one without... ")
print("    ( raw_df_0 -> data_wi_unip_df + data_wo_unip_df ) ")

#====================================================================================================#
#%% Looking at dataset with UniProt Available Only.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Get a dataframe that contains data with uniprot ids only. (Rewritten below)
'''
data_wi_unip_df = copy.deepcopy(raw_df_0)
data_wi_unip_df = data_wi_unip_df[ data_wi_unip_df['uniprot_id'].str.contains('-')                      == False ]
data_wi_unip_df = data_wi_unip_df[ data_wi_unip_df[data_name]   .str.contains("additional information") == False ] 
'''

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Extend the dataset for some test. (Skipped)
'''
# Get combinations of compounds and uniprot ids.
data_wi_unip_df['cmpd_unip_pair'] = data_wi_unip_df['cmpd'] + ',' + data_wi_unip_df['uniprot_id']
data_wi_unip_df['cmpd_unip_orgm'] = data_wi_unip_df['cmpd'] + ',' + data_wi_unip_df['uniprot_id'] + ',' + str(data_wi_unip_df['organism'])

# Get arrays of unique combinations of compounds and uniprot ids.
uniq_cmpd_unip_pair_array = data_wi_unip_df['cmpd_unip_pair'].unique()
uniq_cmpd_unip_orgm_array = data_wi_unip_df['cmpd_unip_orgm'].unique()
'''

#====================================================================================================#
#%%  Rana's extra work of getting a list of UniProt IDs.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Save the dataset that contains data with uniprot ids only. (Skipped)
'''
# Get a list of UniProt IDs.
uniprot_ids_col = raw_df_0['uniprot_id']
uniprot_ids_col = uniprot_ids_col.fillna("")
print("Number of rows of raw data, len(uniprot_ids_col): ", len(uniprot_ids_col))

uniprot_ids_list = [one_uniprot_id for one_uniprot_id in uniprot_ids_col if '-' not in one_uniprot_id]
print("Number of rows of data w/ Uniprot IDs, len(uniprot_ids_list): ", len(uniprot_ids_list))


# Recursively flatten nested list/tuple.
def flatten(container):
    for i in container:
        if isinstance(i, (list,tuple)):
            for j in flatten(i):
                yield j
        else:
            yield i


# ? Why did Rana flatten uniprot_ids_list? 
# ? It doesnt seem to contain nested tuples/lists.
uniprot_ids_list_2 = list(flatten(uniprot_ids_list)) 
print("len(uniprot_ids_list): ", len(uniprot_ids_list_2))
print(uniprot_ids_list_2 == uniprot_ids_list)

# Saving Uniprot IDs
np.savetxt(data_folder / (data_name + '_data_wi_unip.csv'), uniprot_ids_list, delimiter = ",", fmt = '%s')
'''

#====================================================================================================#
#%% Separating data with vs without Uniprot IDs
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Check UniProt IDs for each record in the dataset and split the data accordingly.
data_wi_unip_df = data_wo_unip_df = copy.deepcopy(raw_df_0)
data_wi_unip_df = data_wi_unip_df[ data_wi_unip_df['uniprot_id'].str.contains('-') == False ]
data_wo_unip_df = data_wo_unip_df[ data_wo_unip_df['uniprot_id'].str.contains('-') == True  ]


'''
# Process in a different way if there are "additional information" in the uniprot_id column.
# It turns out that there is NO "additional information" in the uniprot_id column.
data_wi_unip_df = data_wo_unip_df = copy.deepcopy(raw_df_0)

# Data with UniProt IDs. (First, remove "-". Then get "addi info" for data_wo_unip_df. Finally remove "addi info".)
data_wi_unip_df   = data_wi_unip_df[ data_wi_unip_df['uniprot_id'].str.contains('-')                      == False ]

data_wo_unip_df_1 = data_wi_unip_df[ data_wi_unip_df[data_name]   .str.contains("additional information") == True  ]
data_wi_unip_df   = data_wi_unip_df[ data_wi_unip_df[data_name]   .str.contains("additional information") == False ] 

# Data without UniProt IDs.
data_wo_unip_df_2 = data_wo_unip_df[ data_wo_unip_df['uniprot_id'].str.contains('-') == True ]
print("len(data_wo_unip_df_1): ", len(data_wo_unip_df_1), ", those with `additional information`. ")
print("len(data_wo_unip_df_2): ", len(data_wo_unip_df_2), ", those with `-`. ")

data_wo_unip_df   = pd.concat([data_wo_unip_df_1, data_wo_unip_df_2], axis=0)
print("len(data_wo_unip_df): ", len(data_wo_unip_df), ", combination of two. ")
'''

#====================================================================================================#
# Exploding data points with multiple Uniprot IDs listed
print("\n\n"+"-"*90+"\n# Step 1.0 data_wi_unip_df, data_wo_unip_df: ")
print("len(data_wi_unip_df) before/after exploding with multiple Uniprot IDs: ", len(data_wi_unip_df), end = ", ")
data_wi_unip_df['uniprot_id'] = data_wi_unip_df['uniprot_id'].str.split(',')
data_wi_unip_df['uniprot_id_list'] = data_wi_unip_df['uniprot_id']
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
print("Number of unique uniprot ID: ", len(list(sorted(data_wi_unip_df['uniprot_id'].unique()))))



data_wo_unip_df.reset_index(inplace=True)
data_wo_unip_df.drop(columns = ["index", ], inplace = True)
print("\n\n"+"-"*90+"\n# Step 1.2 Dataset without Uniprot IDs: data_wo_unip_df: ")
beautiful_print(data_wo_unip_df)
print("len(data_wo_unip_df): ", len(data_wo_unip_df))















#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#                   .M"""bgd `7MM"""YMM    .g8""8q.   `7MMF'   `7MF'`7MM"""YMM  `7MN.   `7MF'  .g8"""bgd `7MM"""YMM                                    #
#                  ,MI    "Y   MM    `7  .dP'    `YM.   MM       M    MM    `7    MMN.    M  .dP'     `M   MM    `7                                    #
#   pd*"*b.        `MMb.       MM   d    dM'      `MM   MM       M    MM   d      M YMb   M  dM'       `   MM   d                                      #
#  (O)   j8          `YMMNq.   MMmmMM    MM        MM   MM       M    MMmmMM      M  `MN. M  MM            MMmmMM                                      #
#      ,;j9        .     `MM   MM   Y  , MM.      ,MP   MM       M    MM   Y  ,   M   `MM.M  MM.           MM   Y  ,                                   #
#   ,-='    ,,     Mb     dM   MM     ,M `Mb.    ,dP'   YM.     ,M    MM     ,M   M     YMM  `Mb.     ,'   MM     ,M                                   #
#  Ammmmmmm db     P"Ybmmd"  .JMMmmmmMMM   `"bmmd"'      `bmmmmd"'  .JMMmmmmMMM .JML.    YM    `"bmmmd'  .JMMmmmmMMM                                   #
#                                           MMb                                                                                                        #
#                                            `bood'                                                                                                    #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 


print("\n\n\n\n" + "#"*100 + "\n" + "#"*100 +"\n>>> Part 2: Get a small dict (uniprotID -> sequence) for this dataset... ")
print("    ( mainly for data with UniProt ID (data_wi_unip_df) )")

# (1) Uniprot Sequences and associated properties has been downloaded as tsv files (super large files).
#     These tsv files can be used to look up Uniprot IDs and get the sequence we need (a local solution).
# (2) If these Uniprot database files are unavailable, use online solution (APIs, ...) which is slow.
# (3) Save the look up results.



###################################################################################################################
#     `7MMF'        .g8""8q.     .g8"""bgd     db     `7MMF'         `7MM"""Yb. `7MMF' .g8"""bgd MMP""MM""YMM     #
#       MM        .dP'    `YM. .dP'     `M    ;MM:      MM             MM    `Yb. MM .dP'     `M P'   MM   `7     #
#       MM        dM'      `MM dM'       `   ,V^MM.     MM             MM     `Mb MM dM'       `      MM          #
#       MM        MM        MM MM           ,M  `MM     MM             MM      MM MM MM               MM          #
#       MM      , MM.      ,MP MM.          AbmmmqMA    MM      ,      MM     ,MP MM MM.              MM          #
#       MM     ,M `Mb.    ,dP' `Mb.     ,' A'     VML   MM     ,M      MM    ,dP' MM `Mb.     ,'      MM          #
#     .JMMmmmmMMM   `"bmmd"'     `"bmmmd'.AMA.   .AMMA.JMMmmmmMMM    .JMMmmmdP' .JMML. `"bmmmd'     .JMML.        #
###################################################################################################################


#====================================================================================================#
# Read a table contains a large number of AA_seqs info.
# Return a dictionary of UniProt IDs -> sequences.
def look_up_seqs_dict(dictionary         =  "uniprot_annotation_3_201_300.tsv" ,
                      search_unip_list   =  []                                 ,
                      ):

    #--------------------------------------------------#
    # Initialize dicts.
    # Dict to be output, contains the uniprot_ids for those in the dataset only.
    unip_seqs_dict_x = dict([]) 
    # List to be output, contains the uniprot_ids that are NOT found in the dict file.
    missing_unip_list = []
    # Dict that used to store UniProt IDs -> sequences through reading the dict file.
    all_unip_seqs_dict = dict([])
    #--------------------------------------------------#
    # Read the dict file.
    begin_time = time.time()
    
    if dictionary.find(".tsv") != -1:
        uniprot_info_all_df = pd.read_csv(filepath_or_buffer   =   data_folder / ".." /  dictionary,
                                          header               =   0       , 
                                          sep                  =   '\t'    , 
                                          encoding             =   'cp1252', 
                                          low_memory           =   False   , )
    elif dictionary.find(".csv") != -1:
        uniprot_info_all_df = pd.read_csv(filepath_or_buffer   =   data_folder / ".." /  dictionary,
                                          header               =   0        , 
                                          sep                  =   ','      , 
                                          encoding             =   'cp1252' , 
                                          low_memory           =   True    , )
    else:
        print("File Name Error. Expect a csv or tsv file.")

    seqs_list = uniprot_info_all_df["Sequence"].tolist()
    unip_list = uniprot_info_all_df["Entry"].tolist()

    for i, (seqs, unip) in enumerate(zip(seqs_list, unip_list)):
        all_unip_seqs_dict[unip] = seqs
    del seqs_list, unip_list
    #--------------------------------------------------#
    print("len(all_unip_seqs_dict): ", len(all_unip_seqs_dict), ". # of mapping: UniProt IDs -> sequences. ")
    print("read_csv takes ", time.time() - begin_time, "seconds.")

    #--------------------------------------------------#
    # Use all_unip_seqs_dict to look up the list of uniprot_ids to be searched.
    # Get a dictionary for the list of uniprot_ids to be searched.
    for unip in search_unip_list:
        if unip in all_unip_seqs_dict:
            unip_seqs_dict_x[unip] = all_unip_seqs_dict[unip]
        else:
            missing_unip_list.append(unip)

    return unip_seqs_dict_x, missing_unip_list, len(missing_unip_list)


###################################################################################################################
#                `7MMF'        .g8""8q.     .g8""8q. `7MMF' `YMM'     `7MMF'   `7MF'`7MM"""Mq.                    #
#                  MM        .dP'    `YM. .dP'    `YM. MM   .M'         MM       M    MM   `MM.                   #
#                  MM        dM'      `MM dM'      `MM MM .d"           MM       M    MM   ,M9                    #
#                  MM        MM        MM MM        MM MMMMM.           MM       M    MMmmdM9                     #
#                  MM      , MM.      ,MP MM.      ,MP MM  VMA          MM       M    MM                          #
#                  MM     ,M `Mb.    ,dP' `Mb.    ,dP' MM   `MM.        YM.     ,M    MM                          #
#                .JMMmmmmMMM   `"bmmd"'     `"bmmd"' .JMML.   MMb.       `bmmmmd"'  .JMML.                        #
###################################################################################################################

#====================================================================================================#
# Look up Uniprot tsv files dictionaries and Generate a small dictionary for dataset final_unip_seqs_dict.
# Save final_unip_seqs_dict to a csv file after generating it for the first time.
# After the final_unip_seqs_dict.csv is generated for the first time, use it directly.

# Full unip -> seqs dictionary downloaded from UniProtKB is very large. 
# Once a unip -> seqs dict is generated, do not check those large files.
# This is totally different from processing cmpd_name -> smiles, 
#  where dictionaries are NOT that large.

if not Path.exists(data_folder / "final_unip_seqs_dict.csv"):
    # If there is NOT a final dictionary, get one.
    final_unip_seqs_dict = dict([])
    # Get a list of unique uniprot ids.
    unique_unip_list = list(set(data_wi_unip_df['uniprot_id'].tolist()))
    print("len(unique_unip_list): (Number of uniprot IDs in dataset)", len(unique_unip_list))
    # Start with a small dict file (AC -> SQ small dict version 2.1). 
    final_unip_seqs_dict, _, _ = look_up_seqs_dict(dictionary = "uniprot_AC_SQ_2_1.csv"   , search_unip_list = unique_unip_list )
    #--------------------------------------------------#
    # The following dictionaries are more than 25 GB combined. Avoid this step as long as final_unip_seqs_dict.csv is generated.
    try:
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_5.tsv"         , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_4.tsv"         , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_3_001_200.tsv" , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_3_201_300.tsv" , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_3_301_400.tsv" , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_3_401_600.tsv" , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_3_601_800.tsv" , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
        unip_seqs_dict_x, _, _ = look_up_seqs_dict(dictionary = "uniprot_annotation_3_801_inf.tsv" , search_unip_list = unique_unip_list )
        final_unip_seqs_dict = {**final_unip_seqs_dict, **unip_seqs_dict_x}
    except:
        pass
    #--------------------------------------------------#
    # Collect those uniprot ids that are left out by the dictionaries.
    uniprot_IDs_missing_list = []
    for one_unip in unique_unip_list:
        if one_unip not in final_unip_seqs_dict.keys():
            uniprot_IDs_missing_list.append(one_unip)
    print("Number of uniprot IDs not found: ", len(uniprot_IDs_missing_list))

    #--------------------------------------------------#
    # Search for sequences directly in UniProt DB.
    for i, one_unip in enumerate(uniprot_IDs_missing_list):
        print(i, "out of", len(uniprot_IDs_missing_list)) if i % 100 == 0 else 0
        try:
            xml_url = 'https://rest.uniprot.org/uniprotkb/' + str(one_unip) + '.xml'
            xml_response = requests.get(xml_url)
            xml_data = xmltodict.parse(xml_response.content)
            seqs = xml_data["uniprot"]["entry"]["sequence"]['#text']
            final_unip_seqs_dict[one_unip] = seqs
        except:
            print("Unidentified UniProt ID: ", one_unip)

    #--------------------------------------------------#
    # Finally get the final_unip_seqs_dict and save it to a csv file.
    final_unip_seqs_dict_df_data = {'Entry': list(final_unip_seqs_dict.keys()), 'Sequence': list(final_unip_seqs_dict.values())}
    final_unip_seqs_dict_df = pd.DataFrame.from_dict(final_unip_seqs_dict_df_data)
    final_unip_seqs_dict_df.to_csv( data_folder / 'final_unip_seqs_dict.csv' )

#====================================================================================================#
# Use the saved final_unip_seqs_dict.csv to look up seqs directly if there is one.
else: 
    #--------------------------------------------------#
    # Read final_unip_seqs_dict.csv
    final_unip_seqs_dict_df = pd.read_csv(data_folder / 'final_unip_seqs_dict.csv')

    seqs_list = final_unip_seqs_dict_df["Sequence"].tolist()
    unip_list = final_unip_seqs_dict_df["Entry"].tolist()
    final_unip_seqs_dict = dict([])

    for i, (seqs, unip) in enumerate(zip(seqs_list, unip_list)):
        final_unip_seqs_dict[unip] = seqs

    del seqs_list, unip_list
    #print(final_unip_seqs_dict)



#====================================================================================================#
# Part #2 Checkpoint. 
print("\n\n"+"-"*90+"\n# Step 2.1 A small UniprotID -> Sequence dictionary: final_unip_seqs_dict_df: ")
beautiful_print(final_unip_seqs_dict_df)
print("len(final_unip_seqs_dict_df): ", len(final_unip_seqs_dict_df))

del final_unip_seqs_dict_df

print("\n\n"+"-"*90+"\nDone getting the uniprotID -> sequence dict,  final_unip_seqs_dict.")














#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#   pd""b.             .g8"""bgd   .g8""8q.   `7MMM.     ,MMF'`7MM"""Mq.   .g8""8q.   `7MMF'   `7MF'`7MN.   `7MF'`7MM"""Yb.                            #
#  (O)  `8b          .dP'     `M .dP'    `YM.   MMMb    dPMM    MM   `MM..dP'    `YM.   MM       M    MMN.    M    MM    `Yb.                          #
#       ,89          dM'       ` dM'      `MM   M YM   ,M MM    MM   ,M9 dM'      `MM   MM       M    M YMb   M    MM     `Mb                          #
#     ""Yb.          MM          MM        MM   M  Mb  M' MM    MMmmdM9  MM        MM   MM       M    M  `MN. M    MM      MM                          #
#        88          MM.         MM.      ,MP   M  YM.P'  MM    MM       MM.      ,MP   MM       M    M   `MM.M    MM     ,MP                          #
#  (O)  .M'  ,,      `Mb.     ,' `Mb.    ,dP'   M  `YM'   MM    MM       `Mb.    ,dP'   YM.     ,M    M     YMM    MM    ,dP'                          #
#   bmmmd'   db        `"bmmmd'    `"bmmd"'   .JML. `'  .JMML..JMML.       `"bmmd"'      `bmmmmd"'  .JML.    YM  .JMMmmmdP'                            #
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 


print("\n\n\n\n" + "#"*100 + "\n" + "#"*100 +"\n>>> Part 3.1: Get a small dict (Compound Name -> SMILES) for this dataset... ")
print("    ( for BOTH data_wi_unip_df AND data_wo_unip_df )")

# Processing Compounds ()
# This is totally different from UniProt IDs -> sequences.
# Since the compounds -> SMILES dictionaries are relatively small.
# Import all dictionaries and get a comprehensive one.
# Use the combined dictionary to look up SMILES.

###################################################################################################################
#                .g8"""bgd  `7MM"""YMM  MMP""MM""YMM     `7MM"""Yb.   `7MMF'  .g8"""bgd MMP""MM""YMM              #
#              .dP'     `M    MM    `7  P'   MM   `7       MM    `Yb.   MM  .dP'     `M P'   MM   `7              #
#              dM'       `    MM   d         MM            MM     `Mb   MM  dM'       `      MM                   #
#              MM             MMmmMM         MM            MM      MM   MM  MM               MM                   #
#              MM.    `7MMF'  MM   Y  ,      MM            MM     ,MP   MM  MM.              MM                   #
#              `Mb.     MM    MM     ,M      MM            MM    ,dP'   MM  `Mb.     ,'      MM                   #
#                `"bmmmdPY  .JMMmmmmMMM    .JMML.        .JMMmmmdP'   .JMML.  `"bmmmd'     .JMML.                 #
###################################################################################################################


#====================================================================================================#
# Return a small dictionary of Compound Name -> SMILES for this dataset (only contained compounds in this dataset). 
# This approach has used the following resources to search compound names and map it to SMILES.
#   - Dictionary #1: KEGG Database      (KEGG_name_id_dict, KEGG_id_SMILES_by_mol_dict)
#   - Dictionary #2: MetaNetX Database  (nme_smiles_MNXid_dict)
#   - additional_cmpd_smls_dict:
#       - PubChem Database              (pubchem_cmpd_smls_dict)
#       - NIH Cactus                    (cactus_cmpd_smls_dict)

def look_up_cmpd_dict(additional_dict_list  =  [
                                                "cmpd_smls_pubchem.csv" , 
                                                "cmpd_smls_cactus.csv"  , 
                                                #"cmpd_smls_all.csv"    , 
                                               ] ,
                      search_cmpd_list      =  [] ,
                      data_folder           =  data_folder,
                      ):

    # Import dictionaries previously extracted in a separate program.
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Dictionary #1: KEGG Database (Compound Name --> KEGG ID --> SMILES)
    # Get dict for reading ALL_KEGG_DATA_DICT. (Compound Name --> KEGG ID)
    folder = Path("./CMPD_RXN_tools/KEGGScrapSavings")
    filename = "CPD02_ALL_KEGG_DATA_DICT_bk.p"
    ALL_KEGG_DATA_DICT = pickle.load(open( folder / filename, "rb"))

    KEGG_name_id_dict = dict([])
    for one_KEGG_id in ALL_KEGG_DATA_DICT:
        if "NAME" in ALL_KEGG_DATA_DICT[one_KEGG_id]:
            for one_name in ALL_KEGG_DATA_DICT[one_KEGG_id]["NAME"]:
                KEGG_name_id_dict[one_name.lower()] = one_KEGG_id
    #--------------------------------------------------#
    # Get dict for reading KEGG compound IDs. (KEGG ID --> SMILES)
    KEGG_id_SMILES_by_mol_df = pd.read_csv("./CMPD_RXN_tools/KEGGScrapSavings/CPD02_KEGG_id_SMILES_by_mol_file_no_uniq.csv", header = 0)

    KEGGid_list     = KEGG_id_SMILES_by_mol_df["KEGGid"].tolist()
    SMILES_by_mol   = KEGG_id_SMILES_by_mol_df["SMILES_by_mol"].tolist()

    KEGG_id_SMILES_by_mol_dict = dict([])
    for i, id in enumerate(KEGGid_list):
        KEGG_id_SMILES_by_mol_dict[id] = SMILES_by_mol[i]

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Dictionary #2: MetaNetX Database (Compound Name --> SMILES)
    folder = Path("./CMPD_RXN_tools/MNXScrapSavings")
    filename = "CPD01_nme_smiles_MNXid_dict.p"
    nme_smiles_MNXid_dict = pickle.load(open( folder / filename, "rb"))

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Dictionary #3: Import addtional dictionary.
    additional_cmpd_smls_dict = dict([])
    for one_dict_nme in additional_dict_list:
        additional_cmpd_smls_df = pd.read_csv(data_folder / one_dict_nme,
                                              encoding = "ISO-8859-1",
                                              index_col = None, 
                                              header = 0, 
                                              sep = ';',
                                              engine = 'python')

        cmpd_name_list_additional = additional_cmpd_smls_df["CMPD"].tolist()
        cmpd_smls_list_additional = additional_cmpd_smls_df["CMPD_SMILES"].tolist()

        for i, one_cmpd_name in enumerate(cmpd_name_list_additional):
            additional_cmpd_smls_dict[str(one_cmpd_name).lower()] = cmpd_smls_list_additional[i]

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Imported:
    # pubchem_cmpd_smls_dict
    # cactus_cmpd_smls_dict
    # nme_smiles_MNXid_dict
    # KEGG_name_id_dict, KEGG_id_SMILES_by_mol_dict
    # additional_cmpd_smls_dict

    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
    # Look up dictionaries in a certain order.
    # Output a final dictionary for the list of compound names input.
    cmpd_smls_dict_x = dict([])
    missing_cmpd_list = []

    for i, one_cmpd in enumerate(search_cmpd_list):
        append_bool = False
        one_cmpd = str(one_cmpd)
        one_cmpd = one_cmpd.replace(";", "")
        one_cmpd = one_cmpd.lower()

        if (not append_bool) and one_cmpd in additional_cmpd_smls_dict:
            cmpd_smls = additional_cmpd_smls_dict[one_cmpd]
            if cmpd_smls not in ["NA", "", "nan", "None"]:
                cmpd_smls_dict_x[one_cmpd] = cmpd_smls
                append_bool = True

        
        if (not append_bool) and one_cmpd in nme_smiles_MNXid_dict:
            cmpd_smls = nme_smiles_MNXid_dict[one_cmpd][0]
            if cmpd_smls not in ["NA", "", "nan", "None"]:
                cmpd_smls_dict_x[one_cmpd] = cmpd_smls
                append_bool = True

        if (not append_bool) and one_cmpd in KEGG_name_id_dict:
            if KEGG_name_id_dict[one_cmpd] not in ["NA", "", "nan", "None"] and KEGG_name_id_dict[one_cmpd] in KEGG_id_SMILES_by_mol_dict:
                if KEGG_id_SMILES_by_mol_dict[KEGG_name_id_dict[one_cmpd]] != "None":
                    cmpd_smls_dict_x[one_cmpd] = KEGG_id_SMILES_by_mol_dict[KEGG_name_id_dict[one_cmpd]]
                    append_bool = True
                            
        if (not append_bool):
            cmpd_smls_dict_x[one_cmpd] = "None"
            missing_cmpd_list.append(one_cmpd)
            #print("cmpd:"+one_cmpd+", smls: None")
            
    return cmpd_smls_dict_x, missing_cmpd_list, len(missing_cmpd_list)


###################################################################################################################
#                 `7MMF'        .g8""8q.     .g8""8q.   `7MMF' `YMM'     `7MMF'   `7MF'`7MM"""Mq.                 #
#                   MM        .dP'    `YM. .dP'    `YM.   MM   .M'         MM       M    MM   `MM.                #
#                   MM        dM'      `MM dM'      `MM   MM .d"           MM       M    MM   ,M9                 #
#                   MM        MM        MM MM        MM   MMMMM.           MM       M    MMmmdM9                  #
#                   MM      , MM.      ,MP MM.      ,MP   MM  VMA          MM       M    MM                       #
#                   MM     ,M `Mb.    ,dP' `Mb.    ,dP'   MM   `MM.        YM.     ,M    MM                       #
#                 .JMMmmmmMMM   `"bmmd"'     `"bmmd"'   .JMML.   MMb.       `bmmmmd"'  .JMML.                     #
###################################################################################################################
# Getting cmpd_name -> SMILES dictionary.
#====================================================================================================#
# A list of all compounds without duplicates.
uniq_cmpd_list = list(raw_df_0['cmpd'].unique())



# Getting compound name to SMILES dictionary
cmpd_smls_dict, missing_cmpd_list, missing_cmpd_list_len = \
    look_up_cmpd_dict(additional_dict_list  =  ["cmpd_smls_cactus.csv"                                                  ,
                                                "cmpd_smls_pubchem.csv"                                                 , 
                                                "../../Z03_BRENDA_Compound_Info/Z03_CompoundName_LigandID_SMILES.csv"   , 
                                               ]                                                                        ,
                        search_cmpd_list    =  uniq_cmpd_list                                                           ,
                        data_folder         =  data_folder                                                              ,
                        )


#====================================================================================================#
# Get a dataframe with compound names and SMILES.
cmpd_smls_dict_df = pd.DataFrame.from_dict({"CMPD": list(cmpd_smls_dict.keys()), "CMPD_SMILES": cmpd_smls_dict.values()})
cmpd_smls_dict_df = cmpd_smls_dict_df.drop_duplicates(['CMPD'], keep='first')
#cmpd_smls_dict_df.rename(columns={'CMPD': 'Substrate', 'CMPD_SMILES': 'SMILES'}, inplace=True)


#====================================================================================================#
# Part #3 Checkpoint. 

print("\n\n"+"-"*90+"\n# Step 3.1 A small Compound Name -> SMILES dictionary: cmpd_smls_dict_df: ")
beautiful_print(cmpd_smls_dict_df)
print("len(cmpd_smls_dict_df): ", len(cmpd_smls_dict_df))


print("\n\n"+"-"*90+"\nDone getting the {Compound Name -> SMILES} dict,  cmpd_smls_dict.")
print("Number of different compounds in the set: ", len(uniq_cmpd_list))
print("Number of compounds with smiles found:    ", len(cmpd_smls_dict) - missing_cmpd_list_len)
print("Number of compounds with no smiles found: ", missing_cmpd_list_len)















#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 
#         ,AM              .g8""8q.   `7MMF'   `7MF'MMP""MM""YMM `7MM"""Mq. `7MMF'   `7MF'MMP""MM""YMM         ,M' dP         
#        AVMM            .dP'    `YM.   MM       M  P'   MM   `7   MM   `MM.  MM       M  P'   MM   `7         dP .M'    __,  
#      ,W' MM            dM'      `MM   MM       M       MM        MM   ,M9   MM       M       MM           mmmMmmMmm   `7MM  
#    ,W'   MM            MM        MM   MM       M       MM        MMmmdM9    MM       M       MM             MP dP       MM  
#    AmmmmmMMmm          MM.      ,MP   MM       M       MM        MM         MM       M       MM          mmdMmmMmmm     MM  
#          MM     ,,     `Mb.    ,dP'   YM.     ,M       MM        MM         YM.     ,M       MM           ,M' dP        MM  
#          MM     db       `"bmmd"'      `bmmmmd"'     .JMML.    .JMML.        `bmmmmd"'     .JMML.         dP ,M'      .JMML.
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$# 

print("\n\n\n\n" + "#"*100 + "\n" + "#"*100 +"\n>>> Part 4: Process data_wi_unip_df and generate outputs... ")




###################################################################################################################
#         `7MM           mm                                    db                              db                 #
#           MM           MM                                                                                       #
#      ,M""bMM   ,6"Yb.mmMMmm  ,6"Yb.      `7M'    ,A    `MF'`7MM      `7MM  `7MM `7MMpMMMb. `7MM `7MMpdMAo.      #
#    ,AP    MM  8)   MM  MM   8)   MM        VA   ,VAA   ,V    MM        MM    MM   MM    MM   MM   MM   `Wb      #
#    8MI    MM   ,pm9MM  MM    ,pm9MM         VA ,V  VA ,V     MM        MM    MM   MM    MM   MM   MM    M8      #
#    `Mb    MM  8M   MM  MM   8M   MM          VVV    VVV      MM        MM    MM   MM    MM   MM   MM   ,AP      #
#     `Wbmd"MML.`Moo9^Yo.`Mbmo`Moo9^Yo.         W      W     .JMML.      `Mbod"YML.JMML  JMML.JMML. MMbmmd'       #
#                                                                                                   MM            #
#                                                                                                 .JMML.          #
###################################################################################################################
# Look up seqs, and add sequence column, (for those that are NOT identified, assign a np.nan.)
data_wi_unip_unip_list = data_wi_unip_df['uniprot_id'].tolist()

data_wi_unip_seqs_list = []
for unip in data_wi_unip_unip_list:
    if unip in final_unip_seqs_dict:
        data_wi_unip_seqs_list.append(final_unip_seqs_dict[unip])
    else:
        data_wi_unip_seqs_list.append(np.nan)

data_wi_unip_df["Sequence"] = data_wi_unip_seqs_list



data_wi_unip_df.reset_index(inplace=True)
data_wi_unip_df.drop(columns = ["index", ], inplace = True)
print("\n\n"+"-"*90+"\n# Step 4.0 Dataframe after adding sequence column, (only for those rows with a uniprot ID,) data_wi_unip_df: ")
beautiful_print(data_wi_unip_df[["reaction", "ligandID", "Sequence"]])
print("len(data_wi_unip_df): ", len(data_wi_unip_df))



#====================================================================================================#
# Separate data with found sequences and data w/out sequences mapped yet.

print("\n\n"+"-"*90+"\n# Step 4.0.1 Get a dataframe of those data points with uniprot IDs but without a sequence found through searching UniProt DB.: ")
# Get a dataframe of those data points with uniprot IDs but without a sequence found through searching UniProt DB.
data_wi_unip_wo_seqs_df = data_wi_unip_df[data_wi_unip_df['Sequence'].isnull()]
print("len(data_wi_unip_wo_seqs_df): ", len(data_wi_unip_wo_seqs_df), "# data points with uniprot IDs but without a sequence found") 
# len(data_wi_unip_wo_seqs_df): 217

# Get data_wi_seqs_df & data_wo_seqs_df.
data_wo_seqs_df = pd.concat([data_wo_unip_df, data_wi_unip_wo_seqs_df], axis=0)
data_wi_seqs_df = data_wi_unip_df[~data_wi_unip_df['Sequence'].isnull()]

# After processing uniprot_ids, redefine data_wi_seqs_df and data_wo_seqs_df.
print("After processing uniprot_ids, redefine data_wi_seqs_df and data_wo_seqs_df.")
print("len(data_wi_seqs_df): ", len(data_wi_seqs_df)) # len(data_wi_seqs_df): 68316
print("len(data_wo_seqs_df): ", len(data_wo_seqs_df)) # len(data_wo_seqs_df): 107420



###################################################################################################################
#            .g8""8q. `7MMF'   `7MF'MMP""MM""YMM `7MM"""Mq.`7MMF'   `7MF'MMP""MM""YMM                             #
#          .dP'    `YM. MM       M  P'   MM   `7   MM   `MM. MM       M  P'   MM   `7     __,                     #
#          dM'      `MM MM       M       MM        MM   ,M9  MM       M       MM         `7MM                     #
#          MM        MM MM       M       MM        MMmmdM9   MM       M       MM           MM                     #
#          MM.      ,MP MM       M       MM        MM        MM       M       MM           MM                     #
#          `Mb.    ,dP' YM.     ,M       MM        MM        YM.     ,M       MM           MM                     #
#            `"bmmd"'    `bmmmmd"'     .JMML.    .JMML.       `bmmmmd"'     .JMML.       .JMML.                   #
###################################################################################################################
# Output #1 Protein (with known UniProt ID) and Substrate (with known SMILES, most of them has Molfiles available)


#====================================================================================================#
# Add smiles to the dataframe using the dictionary.
data_wi_seqs_df_0 = copy.deepcopy(data_wi_seqs_df)
data_wi_seqs_df_0['smiles'] = data_wi_seqs_df_0["cmpd"].apply(lambda x: cmpd_smls_dict[x.lower()])
data_wi_seqs_df_0 = data_wi_seqs_df_0[data_wi_seqs_df_0['smiles'].str.contains('None') == False]



#====================================================================================================#
# Checkpoint 4.1
print("\n\n"+"-"*90+"\n# Step 4.1 Dataframe after adding smiles string through looking up compound-smiles dictionary, data_wi_seqs_df_0: ")
data_wi_seqs_df_0.reset_index(inplace=True)
data_wi_seqs_df_0.drop(columns = ["index", ], inplace = True)
beautiful_print(data_wi_seqs_df_0[["smiles", "Sequence"]])
print("len(data_wi_seqs_df_0): ", len(data_wi_seqs_df_0))




#====================================================================================================#
# Prepare Output #1
print("\n\n" + "#"*100 + "\n" + "OUTPUT #1 Protein (with known UniProt ID) and Substrate (with known SMILES, most of them has Molfiles available)\n" + "#"*100 + "\n")

def beautiful_print_(df): # DataFrame Printing.
    # Print the dataset in a well-organized format.
    with pd.option_context('display.max_rows'      , 20  , 
                           'display.min_rows'      , 20  , 
                           'display.max_columns'   , 8   , 
                           #"display.max_colwidth" , None,
                           "display.width"         , None,
                           "expand_frame_repr"     , True,
                           "max_seq_items"         , None,):  # more options can be specified
        # Once the display.max_rows is exceeded, 
        # the display.min_rows options determines 
        # how many rows are shown in the truncated repr.
        print(df)
    return 


data_wi_seqs_df_0['uniprot_id_tuple'] = data_wi_seqs_df_0['uniprot_id_list'].apply(tuple)
data_wi_seqs_df_0_1 = copy.deepcopy(data_wi_seqs_df_0.drop_duplicates(subset=['smiles', 'uniprot_id_tuple'], keep = "first"))
data_wi_seqs_df_0_1.drop(columns = ["uniprot_id_list", "uniprot_id", ], inplace = True)


data_wi_seqs_df_0_1.reset_index(inplace=True)
data_wi_seqs_df_0_1.drop(columns = ["index", ], inplace = True)
beautiful_print_(data_wi_seqs_df_0_1)
print("Number of pairs of identified Protein Structures & Substrate SMILES, len(data_wi_seqs_df_0_1): ", len(data_wi_seqs_df_0_1))

data_wi_seqs_df_0_1.to_csv(output_folder / ("BRENDA_" + data_name + "_substrate" + output_1_suffix))
data_wi_seqs_df_0_1.to_pickle(output_folder / ("BRENDA_" + data_name + "_substrate" + output_1_suffix.replace(".csv", ".p")))




###################################################################################################################
#            .g8""8q. `7MMF'   `7MF'MMP""MM""YMM `7MM"""Mq.`7MMF'   `7MF'MMP""MM""YMM                             #
#          .dP'    `YM. MM       M  P'   MM   `7   MM   `MM. MM       M  P'   MM   `7                             #
#          dM'      `MM MM       M       MM        MM   ,M9  MM       M       MM          pd*"*b.                 #
#          MM        MM MM       M       MM        MMmmdM9   MM       M       MM         (O)   j8                 #
#          MM.      ,MP MM       M       MM        MM        MM       M       MM             ,;j9                 #
#          `Mb.    ,dP' YM.     ,M       MM        MM        YM.     ,M       MM          ,-='                    #
#            `"bmmd"'    `bmmmmd"'     .JMML.    .JMML.       `bmmmmd"'     .JMML.       Ammmmmmm                 #
###################################################################################################################
# Output #2 Protein (with known UniProt ID) and Substrate (with known SMILES, most of them has Molfiles available)

data_wi_seqs_df_1 = copy.deepcopy(data_wi_seqs_df_0)
data_wi_seqs_df_1 = copy.deepcopy(data_wi_seqs_df_1.explode('reaction'))
"""
When you use the explode function on a DataFrame column that contains empty lists, 
pandas will replace the empty lists with None.
"""


# Drop a few useless columns
data_wi_seqs_df_1.drop(columns = ["uniprot_id_list", "uniprot_id", ], inplace = True)
data_wi_seqs_df_1.dropna(subset=['reaction'], inplace = True)


#====================================================================================================#
# Checkpoint 4.2
print("\n\n"+"-"*90+"\n# Step 4.2 Dataframe after explode reaction, data_wi_seqs_df_1: ")
data_wi_seqs_df_1.reset_index(inplace=True)
data_wi_seqs_df_1.drop(columns = ["index", ], inplace = True)
beautiful_print(data_wi_seqs_df_1[["reaction", "uniprot_id_tuple"]])
print("len(data_wi_seqs_df_1): ", len(data_wi_seqs_df_1))



#====================================================================================================#
# Checkpoint 4.3
data_wi_seqs_df_1_1 = copy.deepcopy(data_wi_seqs_df_1.drop_duplicates(subset = ['reaction', 'uniprot_id_tuple'], keep = "first"))

print("\n\n"+"-"*90+"\n# Step 4.3 Dataframe after getting unique pairs of `reaction` & `uniprot_id_tuple`, data_wi_seqs_df_1_1: ")
data_wi_seqs_df_1_1.reset_index(inplace=True)
data_wi_seqs_df_1_1.drop(columns = ["index", ], inplace = True)
beautiful_print(data_wi_seqs_df_1_1)
print("Number of pairs of identified Protein Structures & Reaction Strings, len(data_wi_seqs_df_1_1): ", len(data_wi_seqs_df_1_1))




#====================================================================================================#
# 
rctn_list = data_wi_seqs_df_1_1["reaction"].tolist()

rctt_list = []
prod_list = []
identify_bool_list = []

for idx, rctn in enumerate(rctn_list):
    rctt = rctn.split(" = ")[0].split(" + ")
    prod = rctn.split(" = ")[1].split(" + ")

    all_cmpd_of_rctn = rctt + prod
    all_smiles_found_bool = True
    for cp in all_cmpd_of_rctn: 
        if cp.lower() not in cmpd_smls_dict.keys():
            all_smiles_found_bool = False
            break
    
    rctt_list.append(rctt)
    prod_list.append(prod)
    identify_bool_list.append(all_smiles_found_bool)


data_wi_seqs_df_1_1["reactants"]       = rctt_list
data_wi_seqs_df_1_1["products" ]       = prod_list
data_wi_seqs_df_1_1["all_identified"]  = identify_bool_list


#====================================================================================================#
# Checkpoint 4.4


print("\n\n"+"-"*90+"\n# Step 4.4 Dataframe after parsing the reaction string: ")
data_wi_seqs_df_1_1.reset_index(inplace=True)
data_wi_seqs_df_1_1.drop(columns = ["index", ], inplace = True)
beautiful_print(data_wi_seqs_df_1_1)
print("Number of pairs of identified Protein Structures & Reaction Strings, len(data_wi_seqs_df_1_1): ", len(data_wi_seqs_df_1_1))





#====================================================================================================#
# Prepare Output #1
print("\n\n" + "#"*100 + "\n" + "OUTPUT #2 Protein (with known UniProt ID) and Reaction (with known SMILES)\n" + "#"*100 + "\n")
data_wi_seqs_df_1_1.to_csv(output_folder / ("BRENDA_" + data_name + "_reaction" + output_1_suffix))
data_wi_seqs_df_1_1.to_pickle(output_folder / ("BRENDA_" + data_name + "_reaction" + output_1_suffix.replace(".csv", ".p")))




#====================================================================================================#
# Checkpoint 4.5

print("\n\n"+"-"*90+"\n# Step 4.5 Dataframe with all compounds identified: ")
data_wi_seqs_df_1_1 = data_wi_seqs_df_1_1[data_wi_seqs_df_1_1['all_identified'] == True]
data_wi_seqs_df_1_1.reset_index(inplace=True)
data_wi_seqs_df_1_1.drop(columns = ["index", ], inplace = True)
beautiful_print(data_wi_seqs_df_1_1)
print("Number of pairs of identified Protein Structures & Reaction Strings, len(data_wi_seqs_df_1_1): ", len(data_wi_seqs_df_1_1))































