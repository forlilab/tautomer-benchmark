#!/usr/bin/env python

import pandas
from tqdm import tqdm
from rdkit import Chem

df = pandas.read_csv("Tautomer_database_release_3a.csv")

wanted_transforms =  ["PT_06_00", "PT_07_00", "PT_09_00"]

count_total = 0
count_wanted = 0
count_none = 0

hist = {}

for index, row in tqdm(df.iterrows(), total=len(df)):

    for i in range(row["Size"]):
        #k = row["Quantitative_ratio_%d" % (i+1)]
        #q = row["Qualitative_prevalence_%d" % (i+1)]
        if i <= 0:
            continue

        tautomers = [Chem.MolFromSmiles(row["SMILES_%d" % (i+1)]) for i in range(row["Size"])]
        if None in tautomers:
            count_none += row["Size"]

        transform = row["Transf_1_%d" % (i+1)]
        is_wanted = False
        for t in wanted_transforms:
           if t in transform:
               is_wanted = True
               break
        if is_wanted:
            count_wanted += 1
        count_total += 1
        hist.setdefault(transform, 0)
        hist[transform] += 1

fraction = count_wanted / count_total
print(f"{count_wanted=} / {count_total=} ({fraction:.3f})")

sorted_dict = dict(sorted(hist.items(), key=lambda item: item[1], reverse=True))
checksum = 0
i = 0
for key, value in sorted_dict.items():
    i += 1
    print(key, value)
    checksum += value
    if i == 10:
         break
print(f"{checksum=} {count_none=}")
