from rdkit import Chem
from rdkit.Chem import AllChem
import sys
import os

if len(sys.argv) < 2:
    print("usage: pdbgen.py <SMILES string> [output file]")
    exit(1)

mol = Chem.MolFromSmiles(sys.argv[1])
if mol is None:
    print("invalid SMILES string, exiting")
    os._exit(1)

mol = Chem.AddHs(mol)

AllChem.EmbedMolecule(mol, randomSeed=42)
AllChem.UFFOptimizeMolecule(mol)

file = f"{sys.argv[1]}.pdb" if len(sys.argv) < 3 else sys.argv[2]

with open(file, "w") as pdb_file:
    pdb_file.write(Chem.MolToPDBBlock(mol))
