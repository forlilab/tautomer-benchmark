#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
import sys


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
                    print("Sanitization failed", file=sys.stderr)
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

    def get_tautomers(self,
            mol,
            do_second_round=True,
            remove_less_aromatic=True,
            remove_fewer_amides=True):
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
        # remove tautomers that disrupt aromaticity 
        if remove_less_aromatic:
            n_aromatic_atoms = []
            for tautomer in tautomers:
                mol = Chem.MolFromSmiles(tautomer)
                n = sum([atom.GetIsAromatic() for atom in mol.GetAtoms()])
                n_aromatic_atoms.append(n)
            n = sum([atom.GetIsAromatic() for atom in self.input_mol.GetAtoms()])
            n_aromatic_atoms.append(n)
            n = len(tautomers)
            for i in range(n):
                j = n - i - 1
                if n_aromatic_atoms[j] < max(n_aromatic_atoms):
                    tautomers.pop(j)
                    trajectory.pop(j)
        # remove tautomers that get rid of amides
        amide = Chem.MolFromSmarts('O=[CX3][NX3]') # does NOT match 2-Pyridone (intentionally)
        if remove_fewer_amides:
            n_amides = []
            for tautomer in tautomers:
                mol = Chem.MolFromSmiles(tautomer)
                n = len(mol.GetSubstructMatches(amide))
                n_amides.append(n)
            n = len(self.input_mol.GetSubstructMatches(amide))
            n_amides.append(n)
            n = len(tautomers)
            for i in range(n):
                j = n - i - 1
                if n_amides[j] < max(n_amides):
                    tautomers.pop(j)
                    trajectory.pop(j)
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

    def show_transforms3(self, use_mcs=False, highlight=True):
        n = len(self.tautomers)
        mols = [self.input_mol]*n
        labels = ['']*n
        highlights = []
        for index in range(n):
            tautomer = self.tautomers[index]
            mols.append(tautomer)
            rules = self.trajectory[index]
            #label = ', '.join(rules)
            label = ', '.join(rules)
            labels.append(label)
            if not highlight: highlights.append([])
        if use_mcs:
            for mol in mols:
                mcs = rdFMCS.FindMCS((self.input_mol, mol))
                template = Chem.MolFromSmarts(mcs.smartsString)
                AllChem.Compute2DCoords(template)
                AllChem.GenerateDepictionMatching2DStructure(mol, template)
        if highlight:
            pass

        img = Draw.MolsToGridImage(mols, legends=labels,
                highlightAtomLists=highlights,
                subImgSize=(300, 300), molsPerRow=n)
        return img

def cmd_lineparser():
    parser = argparse.ArgumentParser(description='Generate tautomers')
    parser.add_argument('--sdf')
    parser.add_argument('-s', '--smiles')
    parser.add_argument('-r', '--reaction_smarts', default='smirks.txt')
    args = parser.parse_args()
    neither = args.sdf is None and args.smiles is None
    both = args.sdf is not None and args.smiles is not None
    if neither or both: 
        print("Need either --smiles or --sdf", file=sys.stderr)
        sys.exit()
    return args

if __name__ == "__main__":
    import argparse
    args = cmd_lineparser()
    tautomerizer = Tautomerizer(args.reaction_smarts)

    if args.smiles is not None:
        input_mol = Chem.MolFromSmiles(args.smiles)
        products = tautomerizer.get_tautomers(input_mol)
        for mol, rules in zip(tautomerizer.tautomers, tautomerizer.trajectory):
            print(Chem.MolToSmiles(mol), rules)
        img = tautomerizer.show_transforms3(use_mcs=False)
        img.save("tmp.png")

    elif args.sdf is not None:
        import numpy as np
        num_tautomers = []
        mol_names = []
        supp = Chem.SDMolSupplier(args.sdf)
        counter = 0
        for mol in supp:
            name = mol.GetProp('_Name')
            mol_names.append(name)
            products = tautomerizer.get_tautomers(mol)
            n = len(products)
            num_tautomers.append(n)
            counter += 1
            if counter % 2000 == 0:
                print(counter)
                print('max: %d' % max(num_tautomers))
                for i in range(max(num_tautomers)+1):
                    print('%3d: %6d' % (i, num_tautomers.count(i)))

        print('------------------')
        for i in range(max(num_tautomers)+1):
            print('%3d: %6d' % (i, num_tautomers.count(i)))
        print("input molecules: %d" % counter)
        print("new tautomers: %d" % sum(num_tautomers))
