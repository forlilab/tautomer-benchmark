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
        self.trajectory = [] # transformations that yielded each tautomer
        self.tautomers = []  # populated after each call to get_tautomers()

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

    def _remove_duplicates(self, tautomers, trajectory):
        duplicates = []
        duplicates_map = {}
        for i in range(len(tautomers)):
            if i in duplicates: continue 
            for j in range(i+1, len(tautomers)):
                if tautomers[i] == tautomers[j]:
                    duplicates.append(j)
                    duplicates_map[j] = i
        assert(len(set(duplicates)) == len(duplicates)) # no duplicates in duplicates!
        for j in sorted(duplicates)[::-1]:
            tautomers.pop(j)
            i = duplicates_map[j] # the equivalent tautomer that will be kept
            trajectory[i].extend(trajectory[j])
            trajectory.pop(j)
        return

    def get_tautomers(self, mol, do_second_round=True):
        self.input_mol = mol
        self.trajectory = [] # clean up from previous calls
        trajectory = []
        tautomers = []
        input_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
        products = self.apply_rules(mol)
        for (rule_id, mols) in products.items():
            for mol in mols:
                smiles = Chem.MolToSmiles(mol)
                if smiles == input_smiles: continue
                tautomers.append(smiles)
                trajectory.append([rule_id])
        self._remove_duplicates(tautomers, trajectory)
        # second round of transformations: tautomerize the tautomers
        if do_second_round:
            mols = list([Chem.MolFromSmiles(s) for s in tautomers])
            for (index, mol) in enumerate(mols):
                products = self.apply_rules(mol)
                for (rule_id, mols) in products.items():
                    for mol in mols:
                        smiles = Chem.MolToSmiles(mol)
                        if smiles == input_smiles: continue
                        tautomers.append(smiles)
                        rules = []
                        for rule in trajectory[index]:
                            rules.append('%s->%s' % (rule, rule_id))
                        trajectory.append(rules)
            self._remove_duplicates(tautomers, trajectory)
        #print(input_smiles)
        #print(uniq)
        self.trajectory = trajectory # just in case we ever reassign trajectory
        mols = list([Chem.MolFromSmiles(s) for s in tautomers])
        self.tautomers = mols
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

    def show_transforms3(self, use_mcs=True):
        n = len(self.tautomers)
        mols = [self.input_mol]*n
        labels = ['']*n
        for index in range(n):
            tautomer = self.tautomers[index]
            label = self.trajectory[index]
            mols.append(tautomer)
            labels.append('%s' % label)

        if use_mcs:
            for mol in mols:
                mcs = rdFMCS.FindMCS((self.input_mol, mol))
                template = Chem.MolFromSmarts(mcs.smartsString)
                AllChem.Compute2DCoords(template)
                AllChem.GenerateDepictionMatching2DStructure(mol, template)
        img = Draw.MolsToGridImage(mols, legends=labels, subImgSize=(300, 300), molsPerRow=n)
        #img = Draw.MolsToGridImage(mols, subImgSize=(300, 300), molsPerRow=n)
        return img



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
    input_mol = Chem.MolFromSmiles(args.smiles)
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
    img, products = tautomerizer.show_transforms(input_mol)
    print(Chem.MolToSmiles(input_mol), 'input')
    for rule_id in products:
        for mol in products[rule_id]:
            print(Chem.MolToSmiles(Chem.RemoveHs(mol)), rule_id)
    img.save("tmp.png")

    print('---------------------------')
    products = tautomerizer.get_tautomers(input_mol)
    for mol, rules in zip(tautomerizer.tautomers, tautomerizer.trajectory):
        print(Chem.MolToSmiles(mol), rules)
    img = tautomerizer.show_transforms3(use_mcs=False)
    img.save("show.png")

