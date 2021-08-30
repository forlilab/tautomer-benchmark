from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdChemReactions
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdFMCS
import pandas
from PIL import ImageDraw
from PIL import ImageFont
import os

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
        print("Loaded %d rules from %s. Failed: %d" % (len(rxns), filename, n_fail))
        return rxns

    def apply_rules(self, mol, max_iter=1):
        mol = Chem.AddHs(mol)
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
                    raise RuntimeError
                uniq.add(Chem.MolToSmiles(p[0]))
            uniq = [Chem.MolFromSmiles(s) for s in uniq]
            uniq = [Chem.AddHs(mol) for mol in uniq]
            tautomers[rule_id] = uniq
        return tautomers

    def show_transforms(self, mol):
        tautomers = self.apply_rules(mol, max_iter=1)
        mols = [Chem.RemoveHs(mol)]
        labels = ['input']
        for rule_id in tautomers:
            for tautomer in tautomers[rule_id]:
                labels.append(rule_id)
                mols.append(Chem.RemoveHs(tautomer))
        img = Draw.MolsToGridImage(mols, legends=labels, subImgSize=(300, 300), molsPerRow=len(mols))
        return img, tautomers




#kekule_supplier = Chem.ResonanceMolSupplier(mol, Chem.KEKULE_ALL)
#products = tuple()
#products += rxn.RunReactants((mol,))
#for kekule_mol in kekule_supplier:
#    products += rxn.RunReactants((kekule_mol,))
