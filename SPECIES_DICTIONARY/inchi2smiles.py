#!/usr/bin/env python3

import sys
from rdkit import Chem
from rdkit import RDLogger


def inchi_to_smiles(inchi):
  mol = Chem.MolFromInchi(inchi)
  if mol is None:
    return None
  smiles = Chem.MolToSmiles(mol)
  return smiles


def main():
  if len(sys.argv) != 2:
    print("Usage: python inchi2smiles.py \"InChI1;InChI2,...\"")
    sys.exit(1)

  inchi_list = sys.argv[1].split(';')

  for inchi in inchi_list:
    space_free_inchi = inchi.strip()
    smiles = inchi_to_smiles(space_free_inchi)
    if smiles:
      # "InChI=1S/C3H5NO2/c1-2-3-4(5)6/h2H,1,3H2",C=CCN([O])[O]
      print(f"\"{space_free_inchi}\",{smiles}")
    else:
      print(f"Invalid InChI: {smiles}")


if __name__ == "__main__":
  RDLogger.DisableLog('rdApp.*')
  main()
