#!/usr/bin/env python

from rdkit import Chem
import pandas
import sys

df = pandas.read_csv("Tautomer_database_release_3a.csv")
index = int(sys.argv[1])

#smiles = df.iloc[index]['SMILES_1']
#mol = Chem.MolFromSmiles(smiles)
#print(Chem.MolToSmiles(mol))

row = df.iloc[index]
n = row['Size']
tautomers = [Chem.MolFromSmiles(row["SMILES_%d" % (i+1)]) for i in range(n)]
transforms = []
preferences = []
smiles = []
for i in range(n):
    k = row["Quantitative_ratio_%d" % (i+1)]
    q = row["Qualitative_prevalence_%d" % (i+1)]
    c = row["Prevalence_Category_%d" % (i+1)]
    preference = "%s/%s/%s" % (k, q, c)
    preferences.append(preference)
    if i > 0:
        transform = row["Transf_1_%d" % (i+1)]
        smiles_ = Chem.MolToSmiles(tautomers[i])
        transforms.append(transform)
        smiles.append(smiles_)

maxlen_transforms = max([len(t) for t in transforms])
maxlen_preferences = max([len(p) for p in preferences])
transform_str = "%%%ds" % maxlen_transforms
preference_str = "%%%ds" % maxlen_preferences
tautomer_string = transform_str + ' ' + preference_str + ' ' + '%s'
input_string = '%%%ds' % (maxlen_transforms) % 'input' + ' ' + preference_str % preferences[0] + ' ' + '%s'
preferences.pop(0)
print(input_string % (Chem.MolToSmiles(tautomers[0])))
for i in range(n-1):
    print(tautomer_string % (transforms[i], preferences[i], smiles[i]))

