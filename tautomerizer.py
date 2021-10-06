#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
import pandas
from PIL import ImageDraw
from PIL import ImageFont
import os

def get_rdkit_tautomers(mol):
    """ based on: https://gist.github.com/iwatobipen/ca1999b6e4637daf88f315b412220737
    """
    tenum = rdMolStandardize.TautomerEnumerator()
    tenum.Canonicalize(mol)
    res = tenum.Enumerate(mol)
    return list(res)


class Tautomerizer:
    """CACTVS Tautomer rules using RDKit"""

    def __init__(self, smirks_filename):
        self.reactions = self.load_smirks(smirks_filename)

    def load_smirks(self, filename):
        rxns = {}
        n_fail = 0
        with open(filename) as f:
            for line in f:
                if line.startswith("#"): continue
                smirks, rule_id = line.split()
                try:
                    rxn = rdChemReactions.ReactionFromSmarts(smirks)
                    rxns[rule_id] = rxn
                except:
                    print("Failed loading rule %s from %s" % (rule_id, filename))
                    n_fail += 1
        print("Loaded %d rules from %s. Failed reading %d smirks." % (len(rxns), filename, n_fail))
        return rxns

    def apply_rules(self, mol):
        #mol = Chem.AddHs(mol)
        tautomers = {}
        for rule_id in self.reactions:
            rxn = self.reactions[rule_id]
            products = rxn.RunReactants((mol,))
            n = len(products)
            if n == 0: continue
            uniq = set()
            for p in products:
                if len(p) != 1:
                    print("Found product tuple in tuple of tuples with %d molecules" % len(p))
                    print(rule_id, Chem.MolToSmiles(mol), products)
                    #raise RuntimeError
                try:
                    Chem.SanitizeMol(p[0])
                except:
                    print("Sanitization failed")
                    continue
                uniq.add(Chem.MolToSmiles(p[0]))
            uniq = [Chem.MolFromSmiles(s) for s in uniq]
            #uniq = [Chem.AddHs(mol) for mol in uniq if mol is not None]
            uniq = [mol for mol in uniq if mol is not None]
            tautomers[rule_id] = uniq
        return tautomers

    def get_tautomers(self, mol):
        t = self.apply_rules(mol)
        input_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
        uniq = set()
        for rule_id in t:
            for mol in t[rule_id]:
                smiles = Chem.MolToSmiles(mol)
                if smiles == input_smiles: continue
                uniq.add(smiles)
        # now let's do it again
        mols = list([Chem.MolFromSmiles(s) for s in uniq])
        for mol in mols:
            t = self.apply_rules(mol)
            for rule_id in t:
                for mol in t[rule_id]:
                    smiles = Chem.MolToSmiles(mol)
                    if smiles == input_smiles: continue
                    uniq.add(smiles)
        #print(input_smiles)
        #print(uniq)
        mols = list([Chem.MolFromSmiles(s) for s in uniq])
        return mols

    def show_transforms(self, mol):
        tautomers = self.apply_rules(mol)
        mols = [mol]
        labels = ['input']
        for rule_id in tautomers:
            for tautomer in tautomers[rule_id]:
                labels.append(rule_id)
                mols.append(Chem.RemoveHs(tautomer))
        img = Draw.MolsToGridImage(mols, legends=labels, subImgSize=(300, 300), molsPerRow=len(mols))
        return img, tautomers

    def show_transforms2(self, mol):
        """use MCS to align molecules and highlight matched atoms"""
        tautomers = self.apply_rules(mol)
        mols = [Chem.RemoveHs(mol)]
        labels = ['input']
        for rule_id in tautomers:
            for tautomer in tautomers[rule_id]:
                labels.append(rule_id)
                mols.append(Chem.RemoveHs(tautomer))

        mcs = rdFMCS.FindMCS(mols)
        template = Chem.MolFromSmarts(mcs.smartsString)
        AllChem.Compute2DCoords(template)
        for m in mols:
            AllChem.GenerateDepictionMatching2DStructure(m, template)
        img = Draw.MolsToGridImage(mols, legends=labels, subImgSize=(300, 300), molsPerRow=len(mols))
        return img, tautomers


#kekule_supplier = Chem.ResonanceMolSupplier(mol, Chem.KEKULE_ALL)
#products = tuple()
#products += rxn.RunReactants((mol,))
#for kekule_mol in kekule_supplier:
#    products += rxn.RunReactants((kekule_mol,))

def cmd_lineparser():
    parser = argparse.ArgumentParser(description='Generate tautomers')
    parser.add_argument('-s', '--smiles', required=True)
    parser.add_argument('-r', '--reaction_smarts', default='smirks.txt')
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    import argparse
    args = cmd_lineparser()
    tautomerizer = Tautomerizer(args.reaction_smarts)
    mol = Chem.MolFromSmiles(args.smiles)
    #print('mol: ', Chem.MolToSmiles(mol))

    #p = Chem.MolFromSmarts("[nX2,NX2,S,O,Se,Te:1]=,:[C,c,nX2,NX2:6][C,c:5]=,:[C,c,nX2:2][N,n,S,s,O,o,Se,Te:3][#1:4]")
    #print(mol_h.HasSubstructMatch(p))
    #products = tautomerizer.apply_rules(mol)
    #for p in products:
    #    try:
    #        Chem.SanitizeMol(p[0])
    #        print('yupiii')
    #    except:
    #        continue
    #print('here')
    img, products = tautomerizer.show_transforms(mol)
    print(Chem.MolToSmiles(Chem.RemoveHs(mol)), 'input')
    for rule_id in products:
        for mol in products[rule_id]:
            print(Chem.MolToSmiles(Chem.RemoveHs(mol)), rule_id)
    img.save("tmp.png")
