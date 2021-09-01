#!/usr/bin/env python

from rdkit import Chem
from rdkit.Chem import Draw
import sys

mol = Chem.MolFromSmiles(sys.argv[1])
img = Draw.MolToImage(mol)
img.show()
