# Interacted Compound Protein Dataset

This repository contains a data processing pipeline for obtaining a `Interacted Compound Protein Dataset` based on kinetic parameters datasets retrieved from BRENDA database. The final processed dataset contains ~40,000 datapoints with compound structures, protein structures, reaction information (including SMILES representations of reactants and products). BRENDA has some anti-scraping and users cannot access all the data conveniently (for example, CANNOT get ligand ID from ligand name). 


## Z00_Data_BRENDA_cmpd_smiles_dict.py
- Obtain a list of compound names in the dataset we are processing.
- Generata a `Compound Name -> SMILES` dictionary based on (1) NIH Cactus and (2) PubChem.
- The generated dictionary is used in later steps for converting compound names to SMILES.

## Z01_Data_BRENDA_GetReaction.py (~4.5 hours)
- Use BRENDA online queries (SOAP access) to retrieve reaction information (missing in the kinetic parameters datasets).
    - [in]:  ~87000 queries of `organism + EC number + Substrate` 
    - [out]: ~75000 results of `reactions` associated with each combination + `ligand ID's` of the substrates
- Problems:
    - Ligand ID's are identified for ONLY the substrates
        - Solution: Currently trying to use other BRENDA online queries to get more mappings of `compound names` -> `ligand ID`.
        - Solution: Use multiple database to look up compound names (since there are a lot of synonyms.)
    - Multiple molecules in the reaction string 
        - Solution : Identifying the main reactant and product can be helpful.
        - Solution : Identify the Co-factors (i.e., NADPH+, ATP, etc.)
    - Multiple reactions for one group of `organism + EC number + Substrate`
        - `explode()` the dataframe, include these reactions as different instances
        - ~75000 results of `reactions` -> ~130000 single reactions


## Z01_Data_BRENDA_GetSequence_OPTIONAL.py (In Progress)
- Use BRENDA online queries (SOAP access) to retrieve sequences for these instances without a UniProt ID. 
    - Previously use KINMOD to do this step.


## Z02_Data_BRENDA_GetAllPDB.py
- Currently all PDBs are predicted PDB's obtained from https://alphafold.ebi.ac.uk/
    - Currently UniProt ID -> AlphaFold predicted PDB file (single file)
    - Can also do UniProt ID -> PDB ID's (mutliple PDB ID's) -> PDB files (multiple)
        - Problem: Need to decide which PDB file to use.
    - All PDB files are zipped and uploaded, download and run the unzip.py will get all the PDB files in a folder.

## Z02_Download_All_Molfiles.py (~8 hours)
- BRENDA has some anti-scraping to prevent users from downloading using computer program.
- Download all ~260,000 Molfiles on all the pages of https://www.brenda-enzymes.org/ligand.php?brenda_ligand_id={}
- Results are all saved to Z02_All_Ligands_Molfiles
- All MOL files are zipped and uploaded, download and run the unzip.py will get all the MOL files in a folder.

## Z02_Scrape_All_Ligand_Name.py (TAKES AT LEAST 100+ hours)
- BRENDA has some anti-scraping to prevent users from frequently request the ligand webpages.
- Scrape compound names and synonyms on all the pages of https://www.brenda-enzymes.org/ligand.php?brenda_ligand_id={}
- Results will be saved to Z02_All_Ligands_Names

## Z03_Data_BRENDA_LigandID_Compound.py
- Based on the query results, get a dictionary of `compound names` -> `ligand ID` (currently ~22,000 records)
    - Only (~35%) of the compounds in the reaction strings are covered. 
- Get a list of unidentified compound names from the reaction strings (which covers all compounds we are interested)
- Update: Now, we got a dictionary of `compound names` -> `ligand ID` (currently ~110,000 records)
    - 95% of the compounds in the reaction strings are covered. 


## Z04_Data_BRENDA_LigandID_Compound.py
- Obtain a dataset with pairs of interacted protein and reaction.


## Z05_.py (In Progress)
- Fill in those instances where reactions are NOT identified.
    - `EC Number` + `substrate` -> `reaction` (previous idea, only implemented on the KM dataset)
    - `UniProt` + `substrate` -> `reaction` (previous idea, only implemented on the KM dataset)
    - `substrate` + `reaction rules` associated with `EC Number` -> `reaction` (In Progress)
        - This step can be unnecessary since it gets a small number of reactions that are NOT as reliable AND requires using lots of different tools.

































