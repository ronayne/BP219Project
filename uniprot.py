#!/usr/bin/python

import pandas as pd
import numpy as np

cytosol = pd.read_csv("subcell_location_cytosol.tsv", sep='\t')
uniprot = cytosol["Uniprot"]
df = uniprot.dropna()
df1 = df.reset_index(drop=True)
print(df1)
