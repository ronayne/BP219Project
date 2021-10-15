#IMPORTANT NOTE: you will need to provide an absolute path to the directory you want your pdb files stored in as well as an absolute path to the directory your GetArea job submission 
#script is IN THIS SCRIPT.

#To successfully run this script, you will need to ensure the following: 1) you have made an empty directory to hold your fetched pdb files and 2) you are in the directory that contains
#the aforementioned directory holding your pdb files AND your GetArea job submission script AND subcell_location_cytosol.tsv. You may run this script using the following 
command: python ID_cleavage_site.py  


#IMPORTS
import requests
import urllib.request
import urllib.parse
import pandas as pd
import numpy as np
import os

######################################################################### FUNCTIONS

#####FETCH AMINO ACID SEQUENCES #####

#This function takes in a list of UniProt IDs. From there, it pulls down FASTA amino acid sequences corresponding to a set of UniProt IDs from the internet. The output is a list of tuples.
#These tuples are of the form: (UniProt ID, amino acid sequence).

def fetch_seqs(uniprot_IDs):
    seqs = []
    for i in uniprot_IDs:
        URL = 'https://www.uniprot.org/uniprot/' + i + '.fasta'
        page = requests.get(URL)
        fullstring = page.text
        if fullstring[0:4] == '<!DO':
            print('Sequence information for ' + i + ' could not be found')
        else:
            fullstring_split = fullstring.splitlines()
            aa = ''.join(fullstring_split[1:])
            seqs.append((i,aa))
            #print(i + ' sequence extracted')
    return seqs

##### IDENTIFY POTENTIAL PLP SUBSTRATES #####

#This function takes in a list of tuples, as outputted by 'fetch_seqs'. The function then scans the amino acid sequence and determines if a plp cleavage sequence is present. 
#Finally, the function appends the UniProt ID corresponding to a plp substrate to 'ids', which is ultimately outputted.

def id_cleav_site(tuple_list):
    ids = []
    for i in tuple_list:
        ID,seq = i
        for j in range(0,len(seq) - 3):
            if seq[j] + seq[j+2] + seq[j+3] == 'LGG':
                ids.append(ID)
                break
    return ids

##### MAP UNIPROT IDS TO PDB IDS #####

#This function takes a list of uniprot ids that correspond to potential plp substrates (as outputted by id_cleav_site). For all uniprot IDs that have a corresponding PDB ID, 
# said PDB ID is appended to an output vector. 

#NOTE: this code was adapted from UniProt's ID mapping tool (see https://www.uniprot.org/help/api%5Fidmapping)

def uniprot_2_pdb(uniprot_ids):
    url = 'https://www.uniprot.org/uploadlists/'
    
    uniprot_ids_rf = ''
    for i in uniprot_ids:
        uniprot_ids_rf += i + ' '
    uniprot_ids_rf2 = uniprot_ids_rf[0:len(uniprot_ids_rf)-1]

    params = {
    'from': 'ACC+ID',
    'to': 'PDB_ID',
    'format': 'tab',
    'query': uniprot_ids_rf2
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
       response = f.read()
    uniprot_pdbs = response.decode('utf-8')
    uniprot_pdbs_linesplit = uniprot_pdbs.splitlines( )
    pdbs = []
    for i in uniprot_ids:
        for j in uniprot_pdbs_linesplit:
            if j[0:6] == i:
                pdbs.append(j[7:])
                #print('PDB ID found for ' + i)
                break
    return pdbs

##### FETCH PDB FILES #####

#This function takes a list of pdb files (outputted by uniprot_2_pdb) as well as a path to a preexisting directory. With pdb IDs in hand, the function fetches pdb files corresponding to 
#the input pdb IDs (if they exist). These pdb files are held in the specified directory. PDB IDs that have a corresponding PDB file are appended to a new list, which is outputted. 

def fetch_pdbs(pdbs, path2directory):
    #pdbs is an array of PDB identifiers
    pdbs_updated = []
    for i in pdbs:
        URL = 'http://files.rcsb.org/download/' + i + '.pdb'
        try:
            path = urllib.request.urlretrieve(URL)[0]
            new_path =  path2directory + i + '.pdb'
            os.replace(path, new_path)
            #file = open(new_path, 'r')        
            #pdb_file = file.read()
            pdbs_updated.append(i)
            
        except:
            print('No PDB file exists for ' + i)
    return pdbs_updated

##### FILTER PUTATIVE PLP SUBSTRATES BY SOLVENT ACESSIBILITY #####

#This function parses through processed GetArea output and determines if potential plp substrates have solvent exposed recognition sequences. PDB IDs corresponding to sequences with solvent
#accessible cleavage sites are returned. 

def process_SASA_output(pathname,pdb_ids):
    filtered_pdbs_ids = []
    for i in pdb_ids:
        file = open(pathname + i + '.pdb.txt.processed', 'r') 
        output = file.read()
        output_split = output.splitlines()
        for j in range(0, len(output_split)-3):
            if output_split[j] + output_split[j+2] + output_split[j+3] == 'LEUoGLYoGLYo':
                filtered_pdbs_ids.append(i)
                break
    return filtered_pdbs_ids

######################################################################### INTEGRATED PIPELINE

cytosol = pd.read_csv("subcell_location_cytosol.tsv", sep='\t')
uniprot = cytosol["Uniprot"]
df = uniprot.dropna()
df1 = df.reset_index(drop=True)

seqs = fetch_seqs(df1)
putative_substrates = id_cleav_site(seqs)
putative_substrates_pdb_ids = uniprot_2_pdb(putative_substrates)

###DEFINE YOUR ABSOLUTE PATH TO PDB STORAGE DIRECTORY AS 'path2pdbd' HERE (must be a string)

####DEFINE YOUR ABSOLUTE PATH TO GETAREA SCRIPT AS 'path2ga' HERE (must be a string)

updated_ids = fetch_pdbs(putative_substrates_pdb_ids,path2pdbdest)

#call job submission script
for i in updated_ids:
    path = path2pdbdest + i + '.pdb'
    !perl path2ga $path

#call job processing script
for i in updated_ids:
    path = path2pdbdest + i + '.pdb.txt'
    new_file = path2pdbdest + i + '.pdb.txt.processed'
    !cut -b 11-13,71 $path > $new_file

plp_substrates = process_SASA_output(path2pdb, updated_ids)
print(plp_substrates)
