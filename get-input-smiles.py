#!/usr/bin/env python

from rdkit import Chem
import pandas
import sys

df = pandas.read_csv("Tautomer_database_release_3a.csv")

index = int(sys.argv[1])
smiles = df.iloc[index]['SMILES_1']
mol = Chem.MolFromSmiles(smiles)
print(Chem.MolToSmiles(mol))
