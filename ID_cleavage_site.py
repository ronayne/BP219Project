#To run this script, you will need to ensure subcell_location_cytosol.tsv file (can download this from https://www.proteinatlas.org/search/subcell_location%3AAggresome%2CCytosol%2CCytoplasmic+bodies%2CRods+%26+Rings)
#Additionally, this script takes one input: the absolute path to a directory you want to store pdb files of plp substrates. 

#Run this script as follows: python ID_cleavage_site.py pathname=<path/to/directory>

#IMPORTS
import requests
import urllib.request
import urllib.parse
import pandas as pd
import numpy as np

########################FUNCTIONS##############

#Pulls down FASTA sequences of proteins with a given UniProtID. Takes in a list of UniProt IDs as input, outputs a list of tuples. 
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
            print(i + ' sequence extracted')
    return seqs

#Identify PLP recognition sequence

#takes in a list of tuples (ID,seq) as input. This should be outputted from fetch_seqs

def id_cleav_site(tuple_list):
    ids = []
    for i in tuple_list:
        ID,seq = i
        for j in range(0,len(seq) - 3):
            if seq[j] + seq[j+2] + seq[j+3] == 'LGG':
                ids.append(ID)
                break
    return ids

#Get PDB IDs for all proteins with a recognition sequence from UniProt (see https://www.uniprot.org/help/api%5Fidmapping)
#Takes a list of uniprot ids of proteins with a plp cleavage sequence. Outputs a list of pdb identifiers. 

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
                print('PDB ID found for ' + i)
                break
    return pdbs

#Grab pdb files. Input for this function should be both a list of pdb identifiers for plp substrates as well as 
#an absolute path to a directory you want these pdb files downloaded to.

def fetch_pdbs(pdbs, path2directory):
    #pdbs is an array of PDB identifiers
    for i in pdbs:
        URL = 'http://files.rcsb.org/download/' + i + '.pdb'
        try:
            path = urllib.request.urlretrieve(URL)[0]
            new_path =  path2directory + i + '.pdb'
            os.replace(path, new_path)
            #file = open(new_path, 'r')        
            #pdb_file = file.read()cd .
    
        except:
            print('No PDB file exists for ' + i)
    return 

#########IMPLEMENTATION OF FUNCTIONS########

##Read in subcell_location_cytosol.tsv file and pull out UniProt identifiers
cytosol = pd.read_csv("subcell_location_cytosol.tsv", sep='\t')
uniprot = cytosol["Uniprot"]
df = uniprot.dropna()
df1 = df.reset_index(drop=True)

##Use fetch_seqs funtion to retrieve FASTA sequences from UniProt
seqs = fetch_seqs(df1)

##Identify proteins that have a plp cleavage sequences
putative_substrates = id_cleav_site(seqs)

##Find pdb identifiers of plp substrates, if they exist
pdb_ids = uniprot_2_pdb(putative_substrates)

##Retrieve pdb files and place them in an existing directory on your local device

fetch_pdbs(pdb_ids,pathname)


