# -*- coding: utf-8 -*-
"""
Created on Thu May  4 03:26:23 2023

@author: Rana
"""

import pandas as pd
import numpy as np 
import re 
import hashlib
import urllib.request
import ssl
from zeep import Client

km = pd.read_csv('Raw_Data/brenda_km_results_04212023.csv',  header=None, sep='	', encoding='cp1252')
kcatkm = pd.read_csv('Raw_Data/brenda_kcatkm_results_04212023.csv',  header=None, sep='	', encoding='cp1252')
kcat = pd.read_csv('Raw_Data/brenda_kcat_results_04242023.csv',  header=None, sep='	', encoding='cp1252')

df_list = [kcat, kcatkm, km]
df = pd.concat(df_list)
df = df.iloc[:, :-1]
col_names = {0:'EC_Number', 1:'Enzyme', 4:'Substrate', 6:'Organism', 7:'UniProt_ID', 8:'Substrate_ID'}
df = df.rename(columns=col_names)
df = df.drop(columns=[2, 3, 5, 9])
df = df[df['Substrate'].str.contains("additional information")==False] 
df = df.drop_duplicates()
#%% Getting Reactions from EC Number

wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
password = hashlib.sha256("proteins2022_5160!".encode("utf-8")).hexdigest()
client = Client(wsdl)

parameters = ("rana.barghout@mail.utoronto.ca", password, "ecNumber*1.1.1.1", "substrate*13-cis-retinal", "reactionPartners*", "organism*Equus caballus", "ligandStructureId*")
resultString = client.service.getSubstrate(*parameters)
print(resultString)

        
#%% Obtaining unique EC number, substrate, organism combinations to reduce query space

cols = ['EC_Number', 'Substrate', 'Organism']
df['check_combo'] = df[cols].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
query_subject = df['check_combo'].unique()
query = []
for item in query_subject:
    separated_string = item.split("_")
    print(separated_string)
    query.append(separated_string)

del query[57421]
query_subject = np.delete(query_subject,57421)

unique_entries = {inner_list[0] for inner_list in query}

unique_entries = list(unique_entries)

del unique_entries[876]
del unique_entries[4011]

#%% 

all_substrates = []
all_substrate_ids = []
all_products = []
all_product_ids = []

for i in range(len(unique_entries)):
    ec = unique_entries[i]
    subs_params = ("rana.barghout@mail.utoronto.ca", password, f"ecNumber*{ec}", "substrate*", "reactionPartners*", "organism*", "ligandStructureId*")
    prods_params = ("rana.barghout@mail.utoronto.ca", password, f"ecNumber*{ec}", "product*", "reactionPartners*", "organism*", "ligandStructureId*")
    sub_result = client.service.getSubstrate(*subs_params)
    prod_result = client.service.getProduct(*prods_params)
    print("Data number ", i)
    print("----------------")
    
    ## Starting with the substrate results
    print("Starting substrate query....")
    if len(sub_result)>1:
        subs, ids = [], []
        for j in range(len(sub_result)):
            sub = sub_result[j]['substrate']
            sub_id = sub_result[j]['ligandStructureId']
            subs.append(sub)
            ids.append(sub_id)
        all_substrates.append(subs)
        all_substrate_ids.append(ids)
    elif len(sub_result)==1:
        all_substrates.append(sub_result[0]['substrate'])
        all_substrate_ids.append(sub_result[0]['ligandStructureId'])
    else:
        all_substrates.append('None')
        all_substrate_ids.append('None')
        
    ## Now we do the same for all the products
    print("Starting product query....")
    if len(prod_result)>1:
        prods, ids = [], []
        for j in range(len(prod_result)):
            prod = prod_result[j]['product']
            prod_id = prod_result[j]['ligandStructureId']
            prods.append(prod)
            ids.append(prod_id)
        all_products.append(prods)
        all_product_ids.append(ids)
    elif len(prod_result)==1:
        all_products.append(prod_result[0]['product'])
        all_product_ids.append(prod_result[0]['ligandStructureId'])
    else:
        all_products.append('None')
        all_product_ids.append('None')
        
    print("----------------")

#%% Expanding and cleaning nested lists

# Define a function to flatten a nested list
def flatten_list(nested_list):
    # Initialize an empty list to hold the flattened list
    flattened_list = []
    # Loop through each element in the list
    for element in nested_list:
        # If the element is a list, recursively call the function
        if isinstance(element, list):
            flattened_list.extend(flatten_list(element))
        # Otherwise, add the element to the flattened list
        else:
            flattened_list.append(element)
    # Return the flattened list
    return flattened_list

# Call the flatten_list function on the nested list
all_products = flatten_list(all_products)
all_product_ids = flatten_list(all_product_ids)
all_substrates = flatten_list(all_substrates)
all_substrate_ids = flatten_list(all_substrate_ids)

all_compounds = all_substrates + all_products
all_ids = all_substrate_ids + all_product_ids 


for i in range(len(all_compounds)):
    if re.match('^\d+\s', all_compounds[i]):
        all_compounds[i] = re.sub('^\d+\s', '', all_compounds[i])


comps_ids = pd.DataFrame({'Compound_Name': all_compounds, 'Compound_ID': all_ids})
comps_ids = comps_ids.drop_duplicates()

comps_ids.to_csv('compound_to_id.csv')
#compounds_ids = comps_ids.dropna()
 #%% Creating dataframe for compound names and corresponding IDs

compounds_ids = comps_ids.dropna()
unique_compounds = compounds_ids['Compound_ID'].unique()
#%% Loop that Queries BRENDA and downloads all .mol files for the compounds

ssl._create_default_https_context = ssl._create_unverified_context

for i in range(len(unique_compounds)):
    print('----------------')
    print('Data number: ', i)
    compound_id = unique_compounds[i]
    url = f"https://www.brenda-enzymes.org/molfile.php?LigandID={compound_id}"
    filename = f"Molecule_Files/{compound_id}.mol"
    urllib.request.urlretrieve(url, filename)
    print('----------------')

