#!/usr/bin/env python
# coding: utf-8
# author: Zhiqing Xu

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
# Print the DataFrame obtained.
def beautiful_print(df):
    # Print the dataset in a well-organized format.
    with pd.option_context('display.max_rows'      , 20, 
                           'display.min_rows'      , 20, 
                           'display.max_columns'   , 10, 
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
from zeep import Client

import hashlib

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

data_folder      = Path("Z00_BRENDA_Kinetics_Raw/" + data_folder_name)
data_file        = data_file_name

output_folder    = Path("Z04A_DataPreprocessing_Savings/")
output_file      = output_file_name


#====================================================================================================#
# Read the main data file.
raw_df_0 = pd.read_csv(filepath_or_buffer   =   data_folder / data_file      , 
                       on_bad_lines         =   'skip'                       , 
                       index_col            =   None                         , 
                       #names               =   ["EC", data_name, data_name_0, "cmpd", "condition", "organism", "uniprot_id", "publication", "other"], 
                       header               =   None                         , 
                       sep                  =   '	'                        , 
                       encoding             =   'cp1252'                     ,)


raw_df_0.rename(columns = {0 : 'EC'           , 
                           1 : data_name      , 
                           2 : data_name_0    , 
                           3 : 'cmpd'         , 
                           4 : 'condition'    , 
                           5 : 'organism'     , 
                           6 : 'uniprot_id'   , 
                           7 : 'publication'  , 
                           8 : 'other'        , 
                           }                  , 
                inplace = True                ,)


cmpd_list = raw_df_0["cmpd"      ].tolist()
EC_list   = raw_df_0["EC"        ].tolist()
orgm_list = raw_df_0["organism"  ].tolist()


###################################################################################################################
###################################################################################################################
# beautiful_print(raw_df_0)

# unique_cmpd_list = list(set(cmpd_list))
# print("len(unique_cmpd_list): ", len(unique_cmpd_list))

# new_cmpd_list = []
# mew_lgid_list = []

# for idx, cmpd in enumerate(unique_cmpd_list[0:10]):
#     parameters = ("palladium.xu@gmail.com", password, "pyruvic acid")
#     resultString = client.service.getLigandStructureIdByCompoundName(*parameters)
#     print(idx, cmpd, resultString)
#     new_cmpd_list.append(cmpd)
#     new_lgid_list.append(resultString)

# cmpd_lgid_dict = {c : l for (c, l) in zip(new_cmpd_list, new_lgid_list)}
# pickle.dump(cmpd_lgid_dict, open(output_folder / output_file.replace(".csv", "_cmpd_lgid_dict.p"), 'wb'))


###################################################################################################################
###################################################################################################################
# beautiful_print(raw_df_0)


# new_cmpd_list = []
# new_EC_list   = []
# new_orgm_list = []
# new_lgid_list = []



# for idx, (EC, orgm, cmpd) in enumerate(zip(EC_list, orgm_list, cmpd_list)):
#     parameters = ("palladium.xu@gmail.com", password, "ecNumber*", "role*", "ligand*"+"pyruvic acid", "organism*", "ligandStructureId*")
#     resultString = client.service.getLigands(*parameters)
#     print(idx, cmpd, resultString)
#     new_cmpd_list.append(cmpd)
#     new_lgid_list.append(resultString)

# cmpd_lgid_dict = {c : l for (c, l) in zip(new_cmpd_list, new_lgid_list)}
# pickle.dump(cmpd_lgid_dict, open(output_folder / output_file.replace(".csv", "_cmpd_lgid_dict.p"), 'wb'))


###################################################################################################################
###################################################################################################################
beautiful_print(raw_df_0)

new_cmpd_list = []
new_lgid_list = []



for idx in range(1600):
    idx += 10000
    parameters = ("palladium.xu@gmail.com", password, "ecNumber*", "role*", "ligand*", "organism*", "ligandStructureId*"+str(idx))
    resultString = client.service.getLigands(*parameters)
    print(idx, resultString)
    new_cmpd_list.append(resultString)
    new_lgid_list.append(idx)

cmpd_lgid_dict = {c : l for (c, l) in zip(new_cmpd_list, new_lgid_list)}
pickle.dump(cmpd_lgid_dict, open(output_folder / output_file.replace(".csv", "_cmpd_lgid_dict.p"), 'wb'))














































































