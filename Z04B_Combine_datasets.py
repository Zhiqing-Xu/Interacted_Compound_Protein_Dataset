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
# Hard Coding for getting combined IPC datasets (Interacted Proteins and Compounds(Reactions) datasets)

folder        = Path("Z04A_DataPreprocessing_Savings/")
output_folder = Path("Z04B_DataPreprocessing_Savings/")

df_substrate_1 = pd.read_pickle(folder / "BRENDA_kcat_substrate_core.p")    .drop(columns=["kcat",])
df_substrate_2 = pd.read_pickle(folder / "BRENDA_Km_substrate_core.p")      .drop(columns=["Km",])
df_substrate_3 = pd.read_pickle(folder / "BRENDA_kcat_KM_substrate_core.p") .drop(columns=["kcat_KM",])
df_reaction_1  = pd.read_pickle(folder / "BRENDA_kcat_reaction_core.p")     .drop(columns=["kcat",])
df_reaction_2  = pd.read_pickle(folder / "BRENDA_Km_reaction_core.p")       .drop(columns=["Km",])
df_reaction_3  = pd.read_pickle(folder / "BRENDA_kcat_KM_reaction_core.p")  .drop(columns=["kcat_KM",])
   

df_substrate = copy.deepcopy(df_substrate_1)
df_substrate = df_substrate.append(df_substrate_2, ignore_index=True)
df_substrate = df_substrate.append(df_substrate_3, ignore_index=True)
df_substrate = df_substrate.drop_duplicates(subset=['smiles', 'uniprot_id_tuple'], keep = "first")



df_reaction = copy.deepcopy(df_reaction_1)
df_reaction = df_reaction.append(df_reaction_2, ignore_index=True)
df_reaction = df_reaction.append(df_reaction_3, ignore_index=True)
df_reaction = df_reaction.drop_duplicates(subset=['reaction', 'uniprot_id_tuple'], keep = "first")



print("\n\n"+"-"*90+"\n# Step 1 df_substrate: ")
df_substrate.reset_index(inplace=True)
df_substrate.drop(columns = ["index", ], inplace = True)
beautiful_print(df_substrate[["uniprot_id_tuple", "smiles"]])
print("Number of pairs of identified Protein Structures & Substrates, len(df_substrate): ", len(df_substrate))

df_substrate.to_csv(output_folder / ("Z04B_df_substrate.csv"))
df_substrate.to_pickle(output_folder / ("Z04B_df_substrate.p"))


print("\n\n"+"-"*90+"\n# Step 1 df_reaction: ")
df_reaction.reset_index(inplace=True)
df_reaction.drop(columns = ["index", ], inplace = True)
beautiful_print(df_reaction[["reactants", "products", "uniprot_id_tuple"]])
print("Number of pairs of identified Protein Structures & Reactions, len(df_reaction): ", len(df_reaction))


df_reaction.to_csv(output_folder / ("Z04B_df_reaction.csv"))
df_reaction.to_pickle(output_folder / ("Z04B_df_reaction.p"))



























































































