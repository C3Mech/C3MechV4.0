#!/usr/bin/env python3

import sys
from rdkit import Chem
from rdkit import RDLogger


def smiles_to_inchi(smiles):
  mol = Chem.MolFromSmiles(smiles)
  if mol is None:
    return None
  inchi = Chem.MolToInchi(mol)
  return inchi


def main():
  if len(sys.argv) != 2:
    print("Usage: python smiles2inchi.py \"SMILES1,SMILES2,...\"")
    sys.exit(1)

  smiles_list = sys.argv[1].split(',')

  for smiles in smiles_list:
    space_free_smiles = smiles.strip()
    inchi = smiles_to_inchi(space_free_smiles)
    if inchi:
      # "InChI=1S/C3H5NO2/c1-2-3-4(5)6/h2H,1,3H2",C=CCN([O])[O]
      print(f"\"{inchi}\",{space_free_smiles}")
    else:
      print(f"Invalid SMILES: {space_free_smiles}")


if __name__ == "__main__":
  RDLogger.DisableLog('rdApp.*')
  main()
