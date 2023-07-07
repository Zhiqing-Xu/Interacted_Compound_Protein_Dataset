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
    with pd.option_context('display.max_rows'      , 40  , 
                           'display.min_rows'      , 40  , 
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
output_folder = Path("Z04C_DataPreprocessing_Savings/")

data_name           = ["Ki"          , "Km"         , "kcat"        , "kcat_KM"    ][1]
output_suffix       = ['_core.csv'   , '_scrn.csv'  , '_fine.csv'                  ][0] 


#df_CPI = pd.read_csv(folder / "CPI_Km_core.csv", index_col=0, header = 0)
df_CPI = pd.read_csv(folder / f"CPI_{data_name}{output_suffix}", index_col=0, header = 0)
df_CPI['organism'] = df_CPI['organism'].str.lower()



print("\n\n"+"-"*90+"\n# Step 1 df_CPI: ")
df_CPI.reset_index(inplace=True)
df_CPI.drop(columns = ["index", ], inplace = True)
beautiful_print(df_CPI)
print("Number of pairs of identified Protein Structures & Substrates, len(df_CPI): ", len(df_CPI))


df_CPI_count = df_CPI['organism'].value_counts()
print("\n\n"+"-"*90+"\n# Step 1.1 df_CPI value count: ")
df_CPI_count = df_CPI_count.reset_index()
#df_CPI_count.drop(columns = ["index", ], inplace = True)
beautiful_print(df_CPI_count)
print("Number of pairs of identified Protein Structures & Substrates, len(df_CPI): ", len(df_CPI_count))


###################################################################################################################
###################################################################################################################
# Hard Coding for getting combined IPC datasets (Interacted Proteins and Compounds(Reactions) datasets)

df_CPI_sacc = df_CPI[df_CPI['organism'] == "escherichia coli"]

print("\n\n"+"-"*90+"\n# Step 1.2 df_CPI_sacc: ")
df_CPI_sacc.reset_index(inplace=True)
df_CPI_sacc.drop(columns = ["index", ], inplace = True)
beautiful_print(df_CPI_sacc)
print("Number of pairs of identified Protein Structures & Substrates, len(df_CPI_sacc): ", len(df_CPI_sacc))
df_CPI_sacc.to_csv(output_folder / (f"CPI_" + f"sacc_" + f"{data_name}{output_suffix}" ))









































































