#!/usr/bin/env python

from tautomerizer import Tautomerizer
from rdkit import Chem
from rdkit.Chem import Draw
import pandas
from PIL import ImageDraw
from PIL import ImageFont
import os
import sys

if len(sys.argv) == 2:
    smirks_fn = sys.argv[1]
else:
    smirks_fn = "smirks2.txt"
tautomerizer = Tautomerizer(smirks_fn)

font = ImageFont.truetype("FreeMono.ttf", 16)
font_small = ImageFont.truetype("FreeMono.ttf", 16)

df = pandas.read_csv("Tautomer_database_release_3a.csv")


imgsize = 300
count = 0

hits = {}
wanted_transforms = ["PT_07_00"]

success_count = 0
transform_count = 0

for index, row  in df.iterrows():
    tautomers = [Chem.MolFromSmiles(row["SMILES_%d" % (i+1)]) for i in range(row["Size"])]
    if None in tautomers:
        print("Mol is None, index=%d" % index)
        continue
    labels = []
    transforms = []
    selected_mask = [True] # first one is input, so we keep it

    for i in range(row["Size"]):
        k = row["Quantitative_ratio_%d" % (i+1)]
        q = row["Qualitative_prevalence_%d" % (i+1)]
        c = row["Prevalence_Category_%d" % (i+1)]
        labels.append("%s/%s/%s" % (k, q, c))
        if i > 0:
            transform = row["Transf_1_%d" % (i+1)]
            transforms.append(transform)
            is_wanted = transform in wanted_transforms
            selected_mask.append(is_wanted)

    if sum(selected_mask) < 2: continue

    print('index = %d, SMILES_1: %s' % (index, Chem.MolToSmiles(tautomers[0])))
    #mols = [tautomers[i] for i in range(len(tautomers)) if selected_mask[i]]
    #texts = [labels[i] for i in range(len(tautomers)) if selected_mask[i]]
    #img = Draw.MolsToGridImage(mols, molsPerRow=sum(selected_mask),
    #        legends=texts, subImgSize=(imgsize, imgsize))
    #draw = ImageDraw.Draw(img)
    #draw.text((10, 10), row["Solvent"], (0, 0, 0), font)
    ##print(index, len(transforms))
    #for i in range(len(transforms)):
    #    draw.text((imgsize*i+int(.8*imgsize), int(0.7*imgsize)), transforms[i], (0, 0, 0), font_small)
    #    print((int(0.85*imgsize), imgsize*i+int(.8*imgsize)), transforms[i], (0, 0, 0), font_small)

    #img2, ts = tautomerizer.show_transforms(tautomers[0])
    #img.save("png-wanted/%04d-database.png" % (index))
    #img2.save("png-wanted/%04d-rdkit.png" % (index))

    hit_mask = [False for i in range(1, len(selected_mask))]
    predicted_tautomers = tautomerizer.apply_rules(tautomers[0]) 
    for rule_id in predicted_tautomers:
        for t in predicted_tautomers[rule_id]:
            pred_smiles = Chem.MolToSmiles(Chem.RemoveHs(t), isomericSmiles=False)
            for i in range(1, len(selected_mask)):
                if selected_mask[i]:
                    smiles = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(tautomers[i])), isomericSmiles=False)
                    hit_mask[i-1] = hit_mask[i-1] or (smiles == pred_smiles)
                    #print(smiles, pred_smiles)

    #print(index, row["Size"], hit_mask, sum(hit_mask), sum(selected_mask)-1)
    success_count += sum(hit_mask)
    transform_count += sum(selected_mask) - 1

    if index > 2000: break
    #if len(ts) == 0: continue

    #for i in range(len(transforms)):
    #    folder = transforms[i].replace(' ', '').replace('/', '--')
    #    if not os.path.isdir('png/' + folder):
    #        os.mkdir('png/' + folder) 
    #    img.save("png/%s/%04d.png" % (folder, index))
    #print(tautomers, row["Size"])
            
    #if mol.HasSubstructMatch(p1):# or mol.HasSubstructMatch(p2):
    #if mol.HasSubstructMatch(p2):
    #    print(row["Solvent"], row["pH"], row["Transf_1_2"], row["SMILES_1"], row["Quantitative_ratio_1"], row["Qualitative_prevalence_1"])
    #if (index+1) % 50 == 0: print(index)

print('index=%d' % index)
print("success: %d out of %d" % (success_count, transform_count))

