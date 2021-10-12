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

    default_smirks = [
        '[#6X{2-3}z{0-1},N,n,S,s,O,o,Se,Te:1]:,=[NX2,nX2,#6X3,c,P,p:2][Nh,nh,Sh,Oh,Seh,Teh:3]>>[#1][CX4z{0-1},#7,S,O,Se,Te:1][#7X2,CX3z{0-1},#6,P,p:2]=[#7,#16,#8,Se,Te:3] P-06',
        '[CX4z{0-1}h,Nh,nh,Sh,Oh,Seh,Teh:1][NX2,nX2,CX3z{0-1},c,P,p:2]:,=[N,n,S,s,O,o,Se,Te:3]>>[CX{2-3}z{0-1},N,n,S,s,O,o,Se,Te:1]=[NX2,nX2,CX3,c,P,p:2][N,n,S,O,Se,Te:3][#1] P-06r',
        '[nX2,NX2,S,O,Se,Te:1]=,:[C,c,nX2,NX2:6][C,c,NX2,nX2:5]=,:[C,c,nX2,NX2:2][Nh,nh,Sh,sh,Oh,oh,Seh,Teh:3]>>[#1][N,S,O,Se,Te:1][C,NX2:6]=,:[C,N:5][C,#7X2:2]=,:[NX2,S,O,Se,Te:3] P-07',
        '[nX2,NX2,S,O,Se,Te,Cz0X3:1]:,=[c,C,NX2,nX2:6][C,c,NX2,nX2:5]:,=[C,c,NX2,nX2:2][C,c,NX2,nX2:7]:,=[C,c,NX2,nX2:8][Nh,nh,Sh,sh,Oh,oh,Seh,Teh,CX4z0h:3]>>[#1][N,n,S,O,Se,Te,Cz0X4:1][C,c,NX2,nX2:6]=[C,c:5][C,c,NX2,nX2:2]=[C,c,NX2,nX2:7][C,c,NX2,nX2:8]=[NX2,S,O,Se,Te,CX3z0:3] P-09'
        ]

    def __init__(self, smirks_filename=None):
        if smirks_filename is not None:
            with open(smirks_filename) as f:
                smirks_lines = f.readlines()
        else:
            smirks_lines = self.default_smirks
        reactions = self.load_smirks(smirks_lines)
        self.reactions = reactions
        # the following are populated by get_tautomers()
        self.trajectory = []    # transformations that yielded each tautomer
        self.tautomers = []     # mol objects

    def load_smirks(self, smirks_lines):
        rxns = {}
        n_fail = 0
        for line in smirks_lines:
            if line.startswith("#"): continue
            smirks, rule_id = line.split()
            try:
                rxn = rdChemReactions.ReactionFromSmarts(smirks)
                rxns[rule_id] = rxn
            except:
                print("Failed loading rule %s from %s" % (rule_id, filename))
                n_fail += 1
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
                    #print("Sanitization failed", file=sys.stderr)
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

    def show_transforms3(self, use_mcs=False, highlight=False):
        n = len(self.tautomers)
        mols = [self.input_mol]*n
        labels = ['input']*n
        highlights = [[]]*(2*n)
        for index in range(n):
            tautomer = self.tautomers[index]
            mols.append(tautomer)
            rules = self.trajectory[index]
            label = rules[0] 
            if len(rules) > 1:
                label += ' (%d more paths)' % (len(rules) - 1)
            labels.append(label)
        if use_mcs:
            for mol in mols:
                mcs = rdFMCS.FindMCS((self.input_mol, mol))
                template = Chem.MolFromSmarts(mcs.smartsString)
                AllChem.Compute2DCoords(template)
                AllChem.GenerateDepictionMatching2DStructure(mol, template)
        if highlight:
            # https://gist.github.com/iwatobipen/6d8708d8c77c615cfffbb89409be730d
            for index in range(n):
                mcs = rdFMCS.FindMCS([mols[index], mols[index+n]])
                mcs_pattern = Chem.MolFromSmarts(mcs.smartsString)
                input_match = mols[index].GetSubstructMatch(mcs_pattern)
                tauto_match = mols[index+n].GetSubstructMatch(mcs_pattern)
                highlights[index]   = [a.GetIdx() for a in mols[index].GetAtoms()   if a.GetIdx() not in input_match]
                highlights[index+n] = [a.GetIdx() for a in mols[index+n].GetAtoms() if a.GetIdx() not in tauto_match]
        if n == 0:
            labels = ['input']
            mols = [self.input_mol]
            n = 1
        img = Draw.MolsToGridImage(mols, legends=labels,
                highlightAtomLists=highlights,
                subImgSize=(300, 300), molsPerRow=n)
        return img

def cmd_lineparser():
    epilog_msg = (
        'Notes:\n'
        '  - image written to "tautomerizer-output.png" by default\n'
        '    when in single molecule mode (without -m/--multi_molecule)\n'
        '  - RDKit warning are expected\n'
        '\nExamples:\n'
        '    ./tautomerizer.py -s "C1NC=NC=1C(=O)NN2C(=O)NC=C2"\n'
        '    ./tautomerizer.py --sdf molecule.sdf"\n'
        '    ./tautomerizer.py --sdf many_molecules.sdf" -m\n'
    )
        
    parser = argparse.ArgumentParser(
        epilog=epilog_msg,
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--sdf', help='filename of input molecule(s)')
    parser.add_argument('-s', '--smiles')
    parser.add_argument('-r', '--reaction_smarts', help='smirks filename defining transformations')
    parser.add_argument('-m', '--multi_molecule', help='process all molecules in SDF, reports histogram', action='store_true')
    parser.add_argument('-p', '--png_filename', help='graphic depiction of tautomers', default='tautomerizer-output.png')
    args = parser.parse_args()
    neither = args.sdf is None and args.smiles is None
    both = args.sdf is not None and args.smiles is not None
    if neither or both: 
        parser.print_help()
        print("\nError:\n    Need either --smiles or --sdf", file=sys.stderr)
        sys.exit()
    if args.smiles is not None and args.multi_molecule:
        parser.print_help()
        print("\nError:\n    --multi_molecule requires --sdf, not --smiles/-s")
        sys.exit()
    return args

if __name__ == "__main__":
    import argparse
    args = cmd_lineparser()
    tautomerizer = Tautomerizer(args.reaction_smarts)

    if args.multi_molecule:
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
            if counter % 100 == 0:
                print("Processed molecules: %d\r" % counter, end='')
                sys.stdout.flush()
        print("Processed molecules: %d" % counter)
        for i in range(max(num_tautomers)+1):
            print('%3d: %6d' % (i, num_tautomers.count(i)))
        print("new tautomers: %d" % sum(num_tautomers))
    else:
        if args.smiles is not None:
            input_mol = Chem.MolFromSmiles(args.smiles)
        else:
            input_mol = Chem.MolFromMolFile(args.sdf)
        products = tautomerizer.get_tautomers(input_mol)
        print(Chem.MolToSmiles(input_mol), 'input')
        for mol, rules in zip(tautomerizer.tautomers, tautomerizer.trajectory):
            print(Chem.MolToSmiles(mol), rules)
        img = tautomerizer.show_transforms3(use_mcs=False)
        img.save(args.png_filename)
