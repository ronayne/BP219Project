#!/usr/bin/python

import pandas as pd
import numpy as np

#read tsv file of cytosolic proteins
cytosol = pd.read_csv("subcell_location_cytosol.tsv", sep='\t')

#extract Uniprot IDs
uniprot = cytosol["Uniprot"]

#delete rows without Uniprot ID
df = uniprot.dropna()

#reindex
df1 = df.reset_index(drop=True)
print(df1)
