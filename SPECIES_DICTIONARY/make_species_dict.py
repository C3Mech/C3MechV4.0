#!/usr/bin/env python3

# Copyright (c) Raymond Langer 2024

# Please read the README.md for instructions on how to run the script.

# The routines that identify functional groups might not work as expected in 100%
# of cases. It is a quick and dirty implementation to make the species directory
# more readable. Changes might also be required if you add new species.

import os, sys, importlib, re, git, argparse, cairosvg, abc, yaml, pandas, time
import numpy as np
from shutil import copyfile
from pathlib import Path

sys.path.append("..")
write_species = importlib.import_module("PREPROCESSOR.write-species")

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit import RDLogger
from rdkit.Chem import Draw
from rdkit.Chem.Draw import DrawingOptions


def print_not_found(what, dir_or_file, path):
  print(what + " " + dir_or_file + " '" + path + "' does not exist")


def print_error(msg):
  print("\n#error: " + msg)


class SubModelFiles(yaml.YAMLObject):
  yaml_loader = yaml.SafeLoader
  yaml_tag = u'!SubModelFiles'

  def __init__(self, files, output_directory, species_dictionary):
    self.files = copy.copy(files)
    if (files is None):
      self.files = []

    self.output_directory = copy.copy(output_directory)
    if (output_directory is None):
      self.output_directory = ''

    self.species_dictionary = copy.copy(species_dictionary)
    if (species_dictionary is None):
      self.species_dictionary = ''

  def check(self):
    ok = True
    for filename in self.files:
      if (not os.path.isfile(filename)):
        print_not_found("submodel", "file", filename)
        ok = False

    if (not os.path.isdir(self.output_directory)):
      print_not_found("species dictionary output", "directory",
                      self.output_directory)
      ok = False

    if (not os.path.isfile(self.species_dictionary)):
      print_not_found("species dictionary", "file", self.species_dictionary)
      ok = False
    return ok


def read_yaml_input(filename):
  """ returns a SubModelFiles """
  yaml_input_options = None
  with open(filename) as inp:
    try:
      yaml_input_options = yaml.safe_load(inp)
    except ValueError as e:
      util.print_error("invalid syntax in yaml input '" + filename + "'")
  return yaml_input_options


def make_submodelfiles_from_yaml(filename, cmd_output_directory,
                                 cmd_species_dictionary):
  """ 
      This routine returns a checked SubModelFiles object if successful.
      Note: the routine may quit the script. 
  """
  if (not os.path.isfile(filename)):
    print_not_found("yaml input", "file", filename)
    quit()
  print("reading yaml file \'" + filename + "'")
  submodels = read_yaml_input(filename)

  if (cmd_output_directory != ''):
    if (submodels.output_directory != ''):
      print("using output directory '" + cmd_output_directory +
            "' provided from the command line")
    submodels.output_directory = cmd_output_directory

  if (cmd_species_dictionary != ''):
    if (submodels.species_dictionary != ''):
      print("using species dictionary '" + cmd_species_dictionary +
            "' provided from the command line")
    submodels.species_dictionary = cmd_species_dictionary

  if (not submodels.check()):
    print("error: invalid input")
    quit()

  return submodels


def print_success():
  print("Success!\n")


def parse_composition(composition, elements):
  cur_compo = {}
  count_nothing = 0
  for i in range(4):
    c = i * 5
    if re.match("[\\s0]+", composition[c:c + 5]):
      count_nothing += 1
    elif composition[c + 1] != ' ':
      if composition[c:c + 2] in elements:
        n = re.match("[a-zA-Z]+\\s*(\\d+)", composition[c:c + 5])
        if not n:
          print("ERROR: in ", composition[c:])
          quit()
        if composition[c:c + 2] in cur_compo:
          print("ERROR: element", composition[c:c + 2],
                "cannot be specified twice")
          quit()
        cur_compo[composition[c:c + 2]] = int(n.group(1))
      else:
        print("current characters: '" + composition[c:c + 5] + "'")
        print("composition: '" + composition + "'")
        raise Exception("ERROR: unknown element '" + composition[c:c + 2] +
                        "'")
        quit()
    else:
      if composition[c:c + 1] in elements:
        n = re.match("[a-zA-Z]+\\s*(\\d+)", composition[c:c + 5])
        if not n:
          print("ERROR: in ", composition[c:])
          quit()
        if composition[c:c + 1] in cur_compo:
          print("ERROR: element", composition[c:c + 1],
                "cannot be specified twice")
        cur_compo[composition[c:c + 1]] = int(n.group(1))
      else:
        print("current characters: '" + composition[c:c + 5] + "'")
        print("composition: '" + composition + "'")
        raise Exception("ERROR: unknown element '" + composition[c:c + 1] +
                        "'")
        quit()
  if (count_nothing == 4):
    print("composition: '" + composition + "'")
    print("ERROR: empty composition")
    quit()
  return cur_compo


def get_sum_formula(composition):
  #
  # This is used for a sanitity check for the provided InChI's
  #
  # generates a sum formula that can be compared to the
  # sum formula in the InChI's
  elements = {'C': 1, 'H': 1, 'F': 1, 'N': 1, 'O': 1, 'HE': 1, 'AR': 1}
  if not len(composition) == 20:
    print("composition:", composition)
    raise Exception("ERROR: the composition string has the wrong length")
  cur_compo = parse_composition(composition, elements)
  sum_formula = ""
  for elem in elements:
    if (elem in cur_compo):
      if cur_compo[elem] > 1:
        sum_formula += elem + str(cur_compo[elem])
      elif cur_compo[elem] == 1:
        sum_formula += elem
      else:
        if (cur_compo[elem] != 0):
          print("ERROR: found", cur_compo[elem], elem)
          quit()
  if (sum_formula == 'HF'):
    return 'FH'
  return sum_formula


def get_real_inchi(inchi):
  # Underscores (_) are our way to indicate
  # some species are special (for example, triplet and singlet CH2
  # have the same InChI's)
  if (re.match(".+_.+", inchi)):
    species_inchi = re.match("(InChI=[^!\\s_]+)", inchi)
    return species_inchi.group(1)
  return inchi


def print_python_elements(name, composition):
  """ generates a string in the format
  "CH3CHOCHO": [[], {'H': 5, 'O': 2, 'C': 3}],
  from the element composition given in the thermochemistry file
  """

  if not len(composition) == 20:
    print("composition:", composition)
    raise Exception("ERROR: the composition string has the wrong length")
  elements = {'C': 1, 'H': 1, 'O': 1, 'N': 1, 'HE': 1, 'AR': 1}
  cur_compo = parse_composition(composition, elements)
  py_info = "\"" + name + "\" : " + "[[], {"
  for elem in elements:
    if (elem in cur_compo):
      if cur_compo[elem] >= 1:
        py_info += "'" + elem + "':" + str(cur_compo[elem]) + ","
  py_info = py_info[:-1]
  py_info += "}],"
  return py_info


def check_inchi_composition_consistency(inchi, sum_formula, species):
  with_two_slash = re.match(".+=[^/]+\\/([^/]+)\\/.+", inchi)
  with_one_slash = re.match(".+=[^/]+\\/([^/]+)", inchi)
  if with_two_slash:
    if with_two_slash.group(1).upper() != sum_formula:
      print(inchi)
      print("species:", species)
      print("sum_formula:", sum_formula)
      print("inchi_sum_formula", with_two_slash.group(1))
      print("ERROR: sum formulas deviate")
      quit()
  elif with_one_slash:
    if with_one_slash.group(1).upper() != sum_formula:
      print(inchi)
      print("species:", species)
      print("sum_formula:", sum_formula)
      print("inchi_sum_formula", with_one_slash.group(1))
      print("ERROR: sum formulas deviate")
      quit()
  else:
    print("ERROR: cannot extract sum formula from \"" + inchi + "\"")
    quit()


def add_other_inchis_and_smiles(lines, i, cur_inchi, cur_smiles, i_prev,
                                inchis):
  if (not cur_inchi in lines[i]):
    raise Exception("'" + str(cur_inchi) + "' not in " + lines[i])
  if (i <= i_prev):
    raise Exception("something is wrong: i = " + str(i) + " <= " +
                    str(i_prev) + " = i_prev")

  inchis_for_cur_species = {cur_inchi: i}
  smiles_for_cur_species = {cur_smiles: i - 1}
  lines_wo_inchi = 0

  for j in range(max(i_prev, 0), i):
    another_inchi = re.match("!!\\s*(InChI=[^!\\s]+)", lines[j])
    if (another_inchi):
      this_inchi = another_inchi.group(1)
      this_real_inchi = get_real_inchi(this_inchi)

      if (this_real_inchi in inchis_for_cur_species
          and this_real_inchi != cur_inchi):
        raise Exception(
            "inchi list for a lumped species contains the duplicate '" +
            this_real_inchi + "'")

      if (lines_wo_inchi == 0):
        raise Exception("found another inchi directly before '" +
                        this_real_inchi + "'. This is invalid syntax.")

      inchis_for_cur_species[this_real_inchi] = j
      lines_wo_inchi = 0
    else:
      lines_wo_inchi += 1
      # reset inchi because everything found so far does not belong to
      # the current species
      if (lines_wo_inchi > 1):
        inchis_for_cur_species = {cur_inchi: i}

  for extra_inchi in inchis_for_cur_species:
    if (extra_inchi in inchis and cur_inchi != extra_inchi):
      print("ERROR:", extra_inchi,
            "is a duplicate (for a lumped species with multiple inchis)")
      quit()

  smiles_for_cur_species = {}
  if (len(inchis_for_cur_species) > 1):
    print("\nLUMPED")
  for inchi_for_lumped in inchis_for_cur_species:
    line_number = inchis_for_cur_species[inchi_for_lumped]  # 0-based index
    smiles_for_lumped = check_smiles(lines,
                                     line_number,
                                     inchis,
                                     inchi_for_lumped,
                                     silent=True)
    smiles_for_cur_species[smiles_for_lumped] = line_number - 1
    if (len(inchis_for_cur_species) > 1):
      print(inchi_for_lumped, smiles_for_lumped)
  return inchis_for_cur_species, smiles_for_cur_species


def check_inchi(lines, i, inchis):
  inchi = re.match("!!\\s*(InChI=[^!\\s]+)", lines[i])
  if (inchi):
    if re.match(".+skipped.*", inchi.group(1)):
      print("skipped not expected. fixme")
      print(lines[i])
      quit()
      return -1, "", ""
    this_inchi = inchi.group(1)

    if (not re.match("\\s*!.+", lines[i + 1])):
      if (this_inchi in inchis):
        if (get_real_inchi(this_inchi) == this_inchi):
          print("ERROR:", this_inchi, "is a duplicate")
          quit()
      inchis[this_inchi] = 1
      if (re.match("^[\\s]+$", lines[i + 1])):
        print("ERROR: line", i + 2, "below ", this_inchi, "must not be empty")
        quit()
      if (re.match("^[\\s]+[^\\s]+", lines[i + 1])):
        print("ERROR: line", i + 2, "below ", this_inchi,
              "must not start with white spaces")
        quit()
      re_species_name = re.match("^([^\\s]+)", lines[i + 1])

      #print(print_python_elements(re_species_name.group(1), lines[i + 1][24:44]))

      sum_formula = get_sum_formula(lines[i + 1][24:44])
      check_inchi_composition_consistency(this_inchi, sum_formula,
                                          lines[i + 1][0:24])
      return i, this_inchi, re_species_name.group(1)
    if (not re.match("\\s*!.+", lines[i + 1])
        and (len(lines[i + 1]) < 80 or lines[i + 1][79] != '1')):
      print("search for", this_inchi)
      print(
          "ERROR: There must be a NASA polynomial coefficient set below every InChi"
      )
      quit()
  elif (re.match(".*InChI.*", lines[i], re.IGNORECASE)):
    if (not re.match("^[!\\s]+InChI.*", lines[i], re.IGNORECASE)):
      print("line:", lines[i])
      print("ERROR: line", i + 1, "contains an inchi but the format is wrong")
      quit()
  return -1, "", ""


def check_inchi_smiles_consistency(inchi, smiles, species, simplification=""):
  mol = Chem.MolFromSmiles(smiles)
  if (not mol):
    rdkit_error("SMILES" + simplification, smiles, species)
    return True

  inchi_from_smiles = Chem.MolToInchi(mol)
  if (get_real_inchi(inchi) != inchi_from_smiles):
    #print(inchi)
    print("inchi" + simplification + ":", get_real_inchi(inchi))
    print("smiles" + simplification + ": ", smiles)
    print("from smiles", inchi_from_smiles)
    print("ERROR: inchi/smiles inconsistency")
    return True
  return False


def check_smiles(lines, i, inchis, inchi, silent):
  if (i <= 0):
    print("ERROR: linenumber", i, "<= 0")
    quit()
  smiles = re.match("!!\\s*([^!\\s]+)", lines[i - 1])
  if (not smiles):
    print("smiles line:")
    print(lines[i - 1])
    print("inchi line:")
    print(lines[i])
    print("ERROR: previous line", i - 1, " is not a SMILES line")
    quit()

  this_smiles = smiles.group(1)
  if (this_smiles[-1] == ","):
    this_smiles = this_smiles[:-1]

  if (not silent):
    print(this_smiles)

  if (check_inchi_smiles_consistency(inchi, this_smiles, 'unknown')):
    quit()
  return this_smiles


def tex_fig_end(file):
  file.write("\\end{figure}\n\n\n")


def tex_clear_page(file):
  file.write("\\clearpage\n\n\n")


def get_scale(string):
  scale = 1.0
  if (len(string) > 25):
    scale = 0.9
  if (len(string) > 30):
    scale = 0.86
  if (len(string) > 33):
    scale = 0.84
  if (len(string) > 35):
    scale = 0.82
  if (len(string) > 38):
    scale = 0.77
  if (len(string) > 40):
    scale = 0.72
  if (len(string) > 42):
    scale = 0.70
  return scale


def tex_sub_fig(file, filename, plot_per_line, scale, species_name, smiles,
                inchi, line_break):
  width = 1.0 / plot_per_line - 0.001
  file.write("\\begin{subfigure}{" + str(width) + "\\textwidth}\n")
  file.write("\\centering\n")
  file.write("\\captionsetup{justification=centering}\n")
  file.write("{\\catcode`\\#=12 \\includegraphics[scale=" + str(scale) + "]{" +
             "" + filename + "}}\n")
  # tuned to avoid Overfull \hbox
  scale_smiles = get_scale(smiles)
  scale_name = get_scale(species_name)
  scale_inchi = 0.1

  file.write(
      "{\\catcode`\\#=12 \\catcode`\\&=12 \\catcode`\\%=12 \\catcode`\\$=12 \\catcode`\\_=12 \\caption*{\n{"
      + "\\scalebox{" + str(scale_name) + "}{\\mbox{ " + species_name +
      " }} \\\\ \\hspace{\\textwidth}\n" + "\\url{ " + inchi +
      " } \\\\ \\hspace{\\textwidth}\n" +
      #"\\scalebox{" + str(scale_inchi) +"}{\\mbox{\\url{" + inchi + "}}} \\\\ \\hspace{\\textwidth}\n" +
      "\\scalebox{" + str(scale_smiles) + "}{\\mbox{\\url{" + smiles +
      '}}}\n' + '}}\n')
  file.write("}\n")
  file.write("\\end{subfigure}")
  if (line_break):
    file.write("\n%\n")
  else:
    file.write("%\n")


def tex_escape(line):
  special = {'&', '%', '$', '_', '{',
             '}'}  # '#' is handled with the catcode command
  escaped_line = ''
  for c in line:
    if (c in special):
      escaped_line += "\\" + c
    else:
      escaped_line += c
  return escaped_line


def is_hydro_carbon(mol):
  has_C = False
  for atom in mol.GetAtoms():
    if (atom.GetAtomicNum() != 1 and atom.GetAtomicNum() != 6):
      return False
    elif (atom.GetAtomicNum() == 6):
      has_C = True
  return has_C


def get_n_atom(mol, element):
  count = 0
  for atom in mol.GetAtoms():
    if (atom.GetSymbol() == element):
      count += 1
  return count


def only_single_bonds(mol):
  for bond in mol.GetBonds():
    if (bond.GetBondType() != Chem.BondType.SINGLE):
      return False
  return True


def single_bonds_with_n_doubles(mol, n_double=1):
  n_found_double = 0
  for bond in mol.GetBonds():
    if (bond.GetBondType() != Chem.BondType.SINGLE):
      if ((n_found_double < n_double or n_double < 0)
          and bond.GetBondType() == Chem.BondType.DOUBLE):
        n_found_double += 1
      else:
        return False
  if (n_found_double < n_double):
    return False
  return True


def single_bonds_with_n_triples(mol, n_triple=1):
  n_found_triple = False
  for bond in mol.GetBonds():
    if (bond.GetBondType() != Chem.BondType.SINGLE):
      if ((n_found_triple < n_triple or n_triple < 0)
          and bond.GetBondType() == Chem.BondType.TRIPLE):
        n_found_triple += 1
      else:
        return False
  return True


def is_alkane(mol):
  num_rads = Descriptors.NumRadicalElectrons(mol)
  return is_hydro_carbon(mol) and only_single_bonds(mol) and num_rads == 0


def is_alkyl(mol):
  num_rads = Descriptors.NumRadicalElectrons(mol)
  return is_hydro_carbon(mol) and only_single_bonds(mol) and num_rads


def is_alkene(mol):
  num_rads = Descriptors.NumRadicalElectrons(mol)
  return is_hydro_carbon(mol) and single_bonds_with_n_doubles(
      mol) and num_rads == 0


def is_alkenyl(mol):
  num_rads = Descriptors.NumRadicalElectrons(mol)
  return is_hydro_carbon(mol) and single_bonds_with_n_doubles(mol) and num_rads


def is_alkyne(mol):
  num_rads = Descriptors.NumRadicalElectrons(mol)
  return is_hydro_carbon(mol) and single_bonds_with_n_triples(
      mol) and num_rads == 0


def is_alkynyl(mol):
  num_rads = Descriptors.NumRadicalElectrons(mol)
  return is_hydro_carbon(mol) and single_bonds_with_n_triples(mol) and num_rads


def is_polyene(mol):
  num_rads = Descriptors.NumRadicalElectrons(mol)
  return is_hydro_carbon(mol) and single_bonds_with_n_doubles(
      mol, -1) and num_rads == 0


def is_polyenyl(mol):
  num_rads = Descriptors.NumRadicalElectrons(mol)
  return is_hydro_carbon(mol) and single_bonds_with_n_doubles(mol,
                                                              -1) and num_rads


def is_hydro_carbon_mol(mol):
  num_rads = Descriptors.NumRadicalElectrons(mol)
  return is_hydro_carbon(mol) and num_rads == 0


def is_hydro_carbon_rad(mol):
  num_rads = Descriptors.NumRadicalElectrons(mol)
  return is_hydro_carbon(mol) and num_rads


def is_aldehyde(mol):
  pat = Chem.MolFromSmarts("[CX3H1](=[O])")
  pat2 = Chem.MolFromSmarts("[CX3H1](=[O]).[CX3H1](=[O])")
  if (mol.GetSubstructMatch(pat) and get_n_atom(mol, "O") == 1):
    return True
  if (mol.GetSubstructMatch(pat2) and get_n_atom(mol, "O") == 2):
    return True
  else:
    return False


def is_RO(mol):
  pat = Chem.MolFromSmarts("[#6][OX1]")
  if (mol.GetSubstructMatch(pat) and get_n_atom(mol, "O") == 1):
    return True
  else:
    return False


def is_RO2(mol):
  pat = Chem.MolFromSmarts("[OX2][OX1]")
  if (mol.GetSubstructMatch(pat) and get_n_atom(mol, "O") == 2):
    if (Descriptors.NumRadicalElectrons(mol) != 1):
      print("something went wrong in is_RO2")
      quit()
    return True
  else:
    return False


def is_RO2H_mol(mol):
  pat = Chem.MolFromSmarts("[OX2][OH]")
  if (mol.GetSubstructMatch(pat) and get_n_atom(mol, "O") == 2
      and Descriptors.NumRadicalElectrons(mol) == 0):
    return True
  else:
    return False


def is_RO2H_rad(mol):
  pat = Chem.MolFromSmarts("[OX2][OH]")
  if (mol.GetSubstructMatch(pat) and get_n_atom(mol, "O") == 2
      and Descriptors.NumRadicalElectrons(mol) != 0):
    return True
  else:
    return False


def is_ether(mol):
  pat = Chem.MolFromSmarts("[#6][OX2][#6]")
  if (mol.GetSubstructMatch(pat) and get_n_atom(mol, "O") == 1):
    return True
  else:
    return False


def is_ether_ooh(mol):
  pat = Chem.MolFromSmarts("[#6][OX2][#6].[OX2][OH]")
  if (mol.GetSubstructMatch(pat)):
    return True
  else:
    return False


def is_O2QOOH(mol):
  pat_o2 = Chem.MolFromSmarts("[OX2][OX1]")
  if (mol.GetSubstructMatch(pat_o2)):
    pat_ooh = Chem.MolFromSmarts("[OX2][OH]")
    if (mol.GetSubstructMatch(pat_ooh)):
      return True
  return False


def is_keto_hydroperoxide(mol):
  # warning: this also fits aldehydes
  pat_keto = Chem.MolFromSmarts("[#6]=O")
  if (mol.GetSubstructMatch(pat_keto)):
    pat_ooh = Chem.MolFromSmarts("[OX2][OH]")
    if (mol.GetSubstructMatch(pat_ooh)):
      return True
  return False


def is_keto_alcohol(mol):
  pat = Chem.MolFromSmarts("O=[#6].[#6][OH]")
  if (mol.GetSubstructMatch(pat)):
    return True
  return False


def is_ketone(mol):
  pat_keto = Chem.MolFromSmarts("[CX3]=O")
  pat_keto2 = Chem.MolFromSmarts("O=[#6].[#6]=O")
  if (mol.GetSubstructMatch(pat_keto) and not is_aldehyde(mol)
      and get_n_atom(mol, "O") == 1):
    return True
  elif (mol.GetSubstructMatch(pat_keto2) and not is_aldehyde(mol)
        and get_n_atom(mol, "O") == 2):
    return True
  return False


def is_POOH(mol):
  pat_ooh = Chem.MolFromSmarts("[OH][OX2].[OX2][OH]")
  if (mol.GetSubstructMatch(pat_ooh)):
    return True
  return False


def is_alcohol(mol):
  pat_oh = Chem.MolFromSmarts("[OH][#6]")
  if (mol.GetSubstructMatch(pat_oh) and get_n_atom(mol, "O") == 1):
    return True
  return False


def is_alcohol_peroxyl(mol):
  pat = Chem.MolFromSmarts("[OH][#6].[OX2][OX1]")
  if (mol.GetSubstructMatch(pat)):
    return True
  return False


def is_keto_peroxyl(mol):
  pat = Chem.MolFromSmarts("[#6]=O.[OX2][OX1]")
  if (mol.GetSubstructMatch(pat)):
    return True
  return False


def is_keto_oxyl(mol):
  pat = Chem.MolFromSmarts("[#6]=O.[#6][OX1]")
  if (mol.GetSubstructMatch(pat)):
    return True
  return False


def is_hydroperoxide_oxyl(mol):
  pat = Chem.MolFromSmarts("[OH][OX2][#6].[#6][OX1]")
  if (mol.GetSubstructMatch(pat)):
    return True
  return False


def is_alcohol_oxyl(mol):
  pat = Chem.MolFromSmarts("[OH][#6].[#6][OX1]")
  if (mol.GetSubstructMatch(pat)):
    return True
  return False


def is_alcohol_hydroperoxide(mol):
  pat = Chem.MolFromSmarts("[OH][#6].[#6][OX2][OH]")
  if (mol.GetSubstructMatch(pat)):
    return True
  return False


def is_peroxide(mol):
  pat = Chem.MolFromSmarts("[#6][O][O][#6]")
  if (mol.GetSubstructMatch(pat)):
    return True
  return False


def is_aromatic(mol):
  for atom in mol.GetAtoms():
    if (atom.GetIsAromatic()):
      return True
  return False


def draw_species(tex_file, mol, fig_dir, timestamp_species_dict, name,
                 tex_opts, smiles, inchi, species_in_cur_line,
                 rows_in_cur_fig):
  fig_filename = os.path.join(fig_dir, name + ".pdf")
  if (not os.path.isfile(fig_filename)
      or os.path.getmtime(fig_filename) <= timestamp_species_dict):
    Draw.MolToFile(mol, os.path.join(fig_dir, "temp.svg"))
    cairosvg.svg2pdf(url=os.path.join(fig_dir, "temp.svg"),
                     write_to=fig_filename)
  if (species_in_cur_line == tex_opts.species_per_line - 1):
    tex_sub_fig(tex_file, name + ".pdf", tex_opts.species_per_line,
                tex_opts.scale, name, smiles, inchi, True)
    species_in_cur_line = 0
    rows_in_cur_fig += 1
    if (rows_in_cur_fig == tex_opts.rows):
      tex_fig_end(tex_file)
      #tex_clear_page(tex_file)
      tex_fig_begin(tex_file)
      rows_in_cur_fig = 0
  else:
    tex_sub_fig(tex_file, name + ".pdf", tex_opts.species_per_line,
                tex_opts.scale, name, smiles, inchi, False)
    species_in_cur_line += 1
  return species_in_cur_line, rows_in_cur_fig


class SpeciesType:
  """ identifiers and RDKit molecules """

  def __init__(self, name, inchi, smiles, mol, fallback):
    self.name = name
    self.inchi = inchi
    self.smiles = smiles
    self.mol = mol
    self.fallback = fallback


class SpeciesDict:

  def __init__(self, species_name2inchi, species_name2smiles, canonical2data,
               fallback):
    self.species_name2inchi = species_name2inchi
    self.species_name2smiles = species_name2smiles
    self.canonical2data = canonical2data
    self.fallback = fallback

  def get_valid_inchis(self):
    """ returns a dictionary; keys are model names, values are InChIs 
        lists of length 1 unless a species is lumped """
    best_inchi = self.species_name2inchi
    for name in self.species_name2inchi:
      if (self.fallback and name in self.fallback.invalid_stereo_inchi):
        if (len(self.species_name2inchi[name]) != 1):
          raise Exception(
              "At least one of the InChIs in " +
              str(self.species_name2inchi[name]) +
              " is supposed to invalid. This should never hapen for lumped species."
          )
        best_inchi[name] = [self.fallback.species_name2inchi_wo_stereo[name]]
    return best_inchi

  def get_valid_primary_keys(self):
    valid_keys = self.canonical2data.copy()
    # filter invalid inchis
    for primary_key in self.canonical2data:
      name = self.canonical2data[primary_key]['model_name']
      if (self.fallback and name in self.fallback.invalid_stereo_inchi):
        print("REMOVING", primary_key)
        del valid_keys[primary_key]
    return valid_keys

  def get_species_name2primary_keys(self):
    species_name2primary_keys = {}
    primary_keys = self.get_valid_primary_keys()
    for primary_key in primary_keys:
      name = primary_keys[primary_key]['model_name']
      if (not name in species_name2primary_keys):
        species_name2primary_keys[name] = [primary_key]
      else:
        species_name2primary_keys[name].append(primary_key)
    return species_name2primary_keys

  def get_dataframe(self, set_index):
    model_names = []
    inchis = []
    smiles = []
    excited = []
    multiplicity = []
    lumped = []
    stereochemistry = []

    for primary_key in self.canonical2data:
      cur_dict = self.canonical2data[primary_key]
      model_names.append(cur_dict['model_name'])
      inchis.append(primary_key[0])
      smiles.append(cur_dict['smiles'])
      excited.append(primary_key[1])
      multiplicity.append(primary_key[2])
      lumped.append(cur_dict['lumped'])
      stereochemistry.append(cur_dict['stereochemistry'])

    data = {
        'model_name': model_names,
        'inchi': inchis,
        'smiles': smiles,
        'excited': excited,
        'multiplicity': multiplicity,
        'lumped': lumped,
        'stereochemistry': stereochemistry
    }
    df = pandas.DataFrame.from_dict(data)
    if (set_index):
      df.set_index(['inchi', 'excited', 'multiplicity'], inplace=True)
    return df

  def write_to_csv(self, filename):
    df = self.get_dataframe(set_index=False)
    df.to_csv(filename,
              index=False,
              columns=[
                  "model_name", "inchi", "smiles", "excited", "multiplicity",
                  "lumped", "stereochemistry"
              ])
    return df


class FallBack:

  def __init__(self, invalid_stereo_inchi, invalid_stereo_smiles,
               species_name2inchi_wo_stereo, species_name2smiles_wo_stereo):
    self.invalid_stereo_inchi = invalid_stereo_inchi
    self.invalid_stereo_smiles = invalid_stereo_smiles
    self.species_name2inchi_wo_stereo = species_name2inchi_wo_stereo
    self.species_name2smiles_wo_stereo = species_name2smiles_wo_stereo


class Input:
  csv_filename = ""
  yaml_filename = ""
  mech_filename = ""
  thermo_filename = ""
  style_file_dir = ""

  def __init__(self, csv_filename, yaml_filename, mech_filename,
               thermo_filename, style_file_dir):
    self.csv_filename = csv_filename
    self.yaml_filename = yaml_filename
    self.mech_filename = mech_filename
    self.thermo_filename = thermo_filename
    self.style_file_dir = style_file_dir


class Output:
  output_dir = ""
  title = ""
  authors = ""
  sha = ""

  def __init__(self, output_dir, title, authors, sha):
    self.output_dir = output_dir
    self.title = title
    self.authors = authors
    self.sha = sha


class TexOptions:
  rows = 7
  species_per_line = 5
  scale = 0.25
  # distance between rows
  baselinestretch = 0.64
  side_margin = 0.60


class SpeciesClassification(object, metaclass=abc.ABCMeta):

  @abc.abstractmethod
  def write_all(self, tex_file, fig_dir, tex_opts):
    raise NotImplementedError(
        'users must define write_all() to use this base class')

  @abc.abstractmethod
  def add_species(self, dry_run, name, inchis, smiles, fallback_data):
    raise NotImplementedError(
        'users must define add_species() to use this base class')


class DefaultClassification(SpeciesClassification):

  def write_all(self, tex_file, fig_dir, tex_opts, timestamp_species_dict):
    print_species_group(self.c0, "$\\textrm{C}_0$ species", tex_file, fig_dir,
                        tex_opts, timestamp_species_dict, False)
    print_species_group(self.c1, "$\\textrm{C}_1$ species", tex_file, fig_dir,
                        tex_opts, timestamp_species_dict)
    print_species_group(self.c2, "$\\textrm{C}_2$ species", tex_file, fig_dir,
                        tex_opts, timestamp_species_dict)

    print_species_group(self.alkane, "Alkanes", tex_file, fig_dir, tex_opts,
                        timestamp_species_dict)
    print_species_group(self.alkyl, "Alkyl radicals", tex_file, fig_dir,
                        tex_opts, timestamp_species_dict)

    print_species_group(self.alkene, "Alkenes", tex_file, fig_dir, tex_opts,
                        timestamp_species_dict)
    print_species_group(self.alkenyl, "Alkenyl radicals", tex_file, fig_dir,
                        tex_opts, timestamp_species_dict)

    print_species_group(self.alkyne, "Alkynes", tex_file, fig_dir, tex_opts,
                        timestamp_species_dict)
    print_species_group(self.alkynyl, "Alkynyl radicals", tex_file, fig_dir,
                        tex_opts, timestamp_species_dict)

    print_species_group(self.polyene, "Polyenes", tex_file, fig_dir, tex_opts,
                        timestamp_species_dict)
    print_species_group(self.polyenyl, "Polyenyl radicals", tex_file, fig_dir,
                        tex_opts, timestamp_species_dict)

    print_species_group(self.aromatic_mol, "Aromatics (molecules)", tex_file,
                        fig_dir, tex_opts, timestamp_species_dict)
    print_species_group(self.aromatic_rad, "Aromatics (radicals)", tex_file,
                        fig_dir, tex_opts, timestamp_species_dict)

    print_species_group(self.hydro_carbon_mol, "Other hydro-carbon molecules",
                        tex_file, fig_dir, tex_opts, timestamp_species_dict)
    print_species_group(self.hydro_carbon_rad, "Other hydro-carbon radicals",
                        tex_file, fig_dir, tex_opts, timestamp_species_dict)

    print_species_group(self.RO2, "Peroxy radicals", tex_file, fig_dir,
                        tex_opts, timestamp_species_dict)
    print_species_group(self.RO, "Oxy radicals", tex_file, fig_dir, tex_opts,
                        timestamp_species_dict)
    print_species_group(self.RO2H_mol, "Hydroperoxides (molecules)", tex_file,
                        fig_dir, tex_opts, timestamp_species_dict)
    print_species_group(self.RO2H_rad, "Hydroperoxides (radicals)", tex_file,
                        fig_dir, tex_opts, timestamp_species_dict)
    print_species_group(self.O2QOOH, "Peroxy hydroperoxides", tex_file,
                        fig_dir, tex_opts, timestamp_species_dict)
    print_species_group(self.aldehyde, "Aldehydes (molecules and radicals)",
                        tex_file, fig_dir, tex_opts, timestamp_species_dict)
    print_species_group(self.ketone, "Ketones (molecules and radicals)",
                        tex_file, fig_dir, tex_opts, timestamp_species_dict)
    print_species_group(self.ether, "Ethers (molecules and radicals)",
                        tex_file, fig_dir, tex_opts, timestamp_species_dict)
    print_species_group(self.ether_ooh, "Ether-hydroperoxides", tex_file,
                        fig_dir, tex_opts, timestamp_species_dict)
    print_species_group(self.keto_hydroperoxide, "Carbonyl-hydroperoxides",
                        tex_file, fig_dir, tex_opts, timestamp_species_dict)
    print_species_group(self.POOH, "Di-hydroperoxy radicals", tex_file,
                        fig_dir, tex_opts, timestamp_species_dict)
    print_species_group(self.alcohol, "Alcohols (molecules and radicals)",
                        tex_file, fig_dir, tex_opts, timestamp_species_dict)

    print_species_group(self.keto_alcohol, "Miscellaneous", tex_file, fig_dir,
                        tex_opts, timestamp_species_dict)
    print_species_group(self.alcohol_peroxyl, "", tex_file, fig_dir, tex_opts,
                        timestamp_species_dict)
    print_species_group(self.keto_peroxyl, "", tex_file, fig_dir, tex_opts,
                        timestamp_species_dict)
    print_species_group(self.keto_oxyl, "", tex_file, fig_dir, tex_opts,
                        timestamp_species_dict)
    print_species_group(self.hydroperoxide_oxyl, "", tex_file, fig_dir,
                        tex_opts, timestamp_species_dict)
    print_species_group(self.alcohol_oxyl, "", tex_file, fig_dir, tex_opts,
                        timestamp_species_dict)
    print_species_group(self.alcohol_hydroperoxide, "", tex_file, fig_dir,
                        tex_opts, timestamp_species_dict)
    print_species_group(self.peroxide, "", tex_file, fig_dir, tex_opts,
                        timestamp_species_dict)
    print_species_group(self.misc, "", tex_file, fig_dir, tex_opts,
                        timestamp_species_dict)

  def make_rdkit_mol(self, name, inchi, smiles, fallback_data):
    mole = None
    fallback = False
    if (not name in fallback_data.invalid_stereo_inchi):
      mol = Chem.MolFromInchi(inchi)
    elif (not name in fallback_data.invalid_stereo_smiles):
      mol = Chem.MolFromSmiles(smiles)
    elif (name in fallback_data.species_name2inchi_wo_stereo):
      mol = Chem.MolFromInchi(fallback_data.species_name2inchi_wo_stereo[name])
      fallback = False
    elif (name in fallback_data.species_name2smiles_wo_stereo):
      raise Exception(
          "This should never happen. InChIs without stereochemistry must be valid."
      )
      mol = Chem.MolFromSmiles(
          fallback_data.species_name2smiles_wo_stereo[name])
      fallback = False
    else:
      raise Exception("constructing rdkit molecule is impossible for '" +
                      name + "'")

    if (not mol):
      raise Exception(name + ", " + inchi + ", " + smiles)
    return mol, fallback

  def add_species(self, dry_run, name, inchis, smiles_list, fallback_data):
    count = 0
    for inchi, smiles in zip(inchis, smiles_list):
      count += 1
      mol, fallback = self.make_rdkit_mol(name, inchi, smiles, fallback_data)
      if (len(inchis) > 1):
        spec = SpeciesType(name + "_LMPD_" + str(count), inchi, smiles, mol,
                           fallback)
      else:
        spec = SpeciesType(name, inchi, smiles, mol, fallback)
      if (dry_run or get_n_atom(mol, "C") == 0):
        self.c0.append(spec)
      elif (get_n_atom(mol, "C") == 1):
        self.c1.append(spec)
      elif (get_n_atom(mol, "C") == 2):
        self.c2.append(spec)
      elif (is_aromatic(mol)):
        if (Descriptors.NumRadicalElectrons(mol) != 0):
          self.aromatic_rad.append(spec)
        else:
          self.aromatic_mol.append(spec)
      elif (is_alkane(mol)):
        self.alkane.append(spec)
      elif (is_alkyl(mol)):
        self.alkyl.append(spec)
      elif (is_alkene(mol)):
        self.alkene.append(spec)
      elif (is_alkenyl(mol)):
        self.alkenyl.append(spec)
      elif (is_alkyne(mol)):
        self.alkyne.append(spec)
      elif (is_alkynyl(mol)):
        self.alkynyl.append(spec)
      elif (is_polyene(mol)):
        self.polyene.append(spec)
      elif (is_polyenyl(mol)):
        self.polyenyl.append(spec)
      elif (is_hydro_carbon_mol(mol)):
        self.hydro_carbon_mol.append(spec)
      elif (is_hydro_carbon_rad(mol)):
        self.hydro_carbon_rad.append(spec)
      elif (is_RO2(mol)):
        self.RO2.append(spec)
      elif (is_aldehyde(mol)):
        self.aldehyde.append(spec)
      elif (is_RO(mol)):
        self.RO.append(spec)
      elif (is_RO2H_mol(mol)):
        self.RO2H_mol.append(spec)
      elif (is_RO2H_rad(mol)):
        self.RO2H_rad.append(spec)
      elif (is_ether(mol)):
        self.ether.append(spec)
      elif (is_ether_ooh(mol)):
        self.ether_ooh.append(spec)
      elif (is_O2QOOH(mol)):
        self.O2QOOH.append(spec)
      elif (is_keto_hydroperoxide(mol)):
        self.keto_hydroperoxide.append(spec)
      elif (is_keto_alcohol(mol)):
        self.keto_alcohol.append(spec)
      elif (is_ketone(mol)):
        self.ketone.append(spec)
      elif (is_POOH(mol)):
        self.POOH.append(spec)
      elif (is_alcohol(mol)):
        self.alcohol.append(spec)
      elif (is_alcohol_peroxyl(mol)):
        self.alcohol_peroxyl.append(spec)
      elif (is_keto_peroxyl(mol)):
        self.keto_peroxyl.append(spec)
      elif (is_keto_oxyl(mol)):
        self.keto_oxyl.append(spec)
      elif (is_hydroperoxide_oxyl(mol)):
        self.hydroperoxide_oxyl.append(spec)
      elif (is_alcohol_oxyl(mol)):
        self.alcohol_oxyl.append(spec)
      elif (is_alcohol_hydroperoxide(mol)):
        self.alcohol_hydroperoxide.append(spec)
      elif (is_peroxide(mol)):
        self.peroxide.append(spec)
      else:
        self.misc.append(spec)

  c0 = []
  c1 = []
  c2 = []
  aromatic_mol = []
  aromatic_rad = []
  alkane = []
  alkyl = []
  alkene = []
  alkenyl = []
  alkyne = []
  alkynyl = []
  polyene = []
  polyenyl = []
  hydro_carbon_mol = []
  hydro_carbon_rad = []
  RO2 = []
  aldehyde = []
  RO = []
  RO2H_mol = []
  RO2H_rad = []
  ether = []
  ether_ooh = []
  O2QOOH = []
  keto_hydroperoxide = []
  keto_alcohol = []
  ketone = []
  POOH = []
  alcohol = []
  alcohol_peroxyl = []
  keto_peroxyl = []
  keto_oxyl = []
  hydroperoxide_oxyl = []
  alcohol_oxyl = []
  alcohol_hydroperoxide = []
  peroxide = []
  misc = []


def tex_fig_begin(file):
  file.write("\n\n\\begin{figure}[!ht]\n")
  file.write("\\centering\n")


def print_species_group(species,
                        title,
                        tex_file,
                        fig_dir,
                        tex_opts,
                        timestamp_species_dict,
                        newpage=True):
  if (not len(species)):
    return
  if (title != ""):
    if (newpage):
      tex_file.write("\\clearpage\n")
    tex_file.write("\\section{" + title + "}")

  tex_fig_begin(tex_file)
  species_in_cur_line = 0
  rows_in_cur_fig = 0
  for i in range(len(species)):
    species_in_cur_line, rows_in_cur_fig = draw_species(
        tex_file, species[i].mol, fig_dir, timestamp_species_dict,
        species[i].name, tex_opts, species[i].smiles, species[i].inchi,
        species_in_cur_line, rows_in_cur_fig)
  tex_fig_end(tex_file)


def parse_species_line(species_dict, line):
  without_comment = re.match("\\s*([^!]+)\\s*", line)
  if (without_comment):
    species = without_comment.group(1).strip().split()
    for s in species:
      s = s.upper()
      if (s != "SPECIES"):
        species_dict[s] = 0
  return species_dict


def make_species_from_kinetics_file(lines):
  species_dict = {}
  is_species_section = False
  for i in range(len(lines)):
    if (re.match("\\s*SPECIES", lines[i])):
      is_species_section = True
      print("start species")
    if (is_species_section and re.match("\\s*END", lines[i])):
      break
    if (is_species_section):
      species_dict = parse_species_line(species_dict, lines[i])
  return species_dict


def get_thermo_lines(thermo_filename):
  check_file = Path(thermo_filename)
  lines_thermo = []
  if check_file.is_file():
    thermo_file = open(thermo_filename)
    lines_thermo = thermo_file.readlines()
    thermo_file.close()
  else:
    print_not_found("chemkin thermochemistry", "file", thermo_filename)
    quit()
  return lines_thermo


def itv_mode_check(dry_run, species_list, lines, sc):
  inchis = {}
  count_species = 0
  i_prev = -1
  i_cur = -1
  species_name2inchi = {}
  species_name2smiles = {}
  for i in range(len(lines)):
    inchi_line, inchi, name = check_inchi(lines, i, inchis)
    inchi = get_real_inchi(inchi)
    if (inchi_line >= 0):
      if (name.upper() not in species_list):
        print("'" + name + "' is not in species_list")
      else:
        species_list[name] = 1

      count_species += 1
      smiles = check_smiles(lines, i, inchis, inchi, True)
      #print(name, inchi, smiles)

      i_prev = i_cur
      i_cur = i
      all_inchis, all_smiles = add_other_inchis_and_smiles(
          lines, i, inchi, smiles, i_prev + 1, inchis)
      if (len(all_inchis) > 1):
        print("'" + name + "' is lumped")

      species_name2inchi[name] = all_inchis
      species_name2smiles[name] = all_smiles

      sc.add_species(dry_run, name, [inchi], [smiles], FallBack({}, {}, {},
                                                                {}))

  return sc, count_species, species_name2inchi, species_name2smiles


def strip_nonascii(s):
  return s.encode('ascii', 'ignore').decode()


def check_composition_line(lines, i, inchis):
  tmp_line = strip_nonascii(lines[i])
  clean_line = tmp_line.split('!')[0].rstrip()
  if len(clean_line) >= 80 and clean_line[79] in ['1', '2', '3', '4']:
    if line[79] == '4':
      try:
        label, thermo, comp = read_NASA7_entry(thermo, TintDefault, comments)
      except Exception as e:
        error_line_number = self.line_number - len(current) + 1
        error_entry = ''.join(current).rstrip()
        logger.info('Error while reading thermo entry starting on line {0}:\n'
                    '"""\n{1}\n"""'.format(error_line_number, error_entry))
  inchi = re.match("!!\\s*(InChI=[^!\\s]+)", lines[i])


def get_species_name2sum_formula(lines, species_list_or_dict, thermo_filename,
                                 silent):
  """ returns a dict """
  clean_lines = []
  for i in range(len(lines)):
    tmp_line = strip_nonascii(lines[i])
    clean_lines.append(tmp_line.split('!')[0].rstrip())

  identifier2sum_formula = {}
  for i in range(3, len(clean_lines)):
    if (len(clean_lines[i]) >= 80 and clean_lines[i][79] == '4'
        and len(clean_lines[i - 1]) >= 80 and clean_lines[i - 1][79] == '3'
        and len(clean_lines[i - 2]) >= 80 and clean_lines[i - 2][79] == '2'
        and len(clean_lines[i - 3]) >= 80 and clean_lines[i - 3][79] == '1'):
      the_line = clean_lines[i - 3]
      identifier = the_line[0:24].split(' ')[0].rstrip().upper()
      composition = the_line[24:44]
      sum_formula = get_sum_formula(composition)
      identifier2sum_formula[identifier] = sum_formula

  missing_nasa = {}
  found_nasa = {}
  for identifier in species_list_or_dict:
    if (identifier not in identifier2sum_formula):
      missing_nasa[identifier] = 1
    else:
      found_nasa[identifier] = 1

  if (not silent and len(species_list_or_dict)):
    print("found NASA polynomial coefficients for " + str(len(found_nasa)) +
          "/" + str(len(species_list_or_dict)) + " species")
  if (len(missing_nasa) > 0 and len(species_list_or_dict)):
    print(
        "\n#warning: Could not find thermochemistry data for the following " +
        str(len(missing_nasa)) + " species:")
    for identifier in missing_nasa:
      print(identifier)
    print(
        "Thermochemistry data are expected to be available for each of the species. Is the thermochemistry file '"
        + thermo_filename +
        "' intended for the provided kinetic (sub)model(s)?\n")
  return identifier2sum_formula, found_nasa


def get_species_name2sum_formula_from_file(thermo_filename):
  print("reading the thermochemistry file '" + thermo_filename + "'")
  lines = get_thermo_lines(thermo_filename)
  species_name2sum_formula, found_nasa = get_species_name2sum_formula(
      lines, {}, thermo_filename, silent=True)
  if (len(species_name2sum_formula) == 0):
    print_warning("could not find thermochemistry data for any species")
  return species_name2sum_formula


def print_df_duplicates(df, col):
  df = df[df.duplicated(subset=col, keep=False)]
  print(df.to_string(index=False))
  #v = df[col].value_counts()
  #print(df[df[col].isin(v.index[v.gt(1)])].to_string(index=False))


def check_required_column(df, col):
  if (not col in df):
    print("read the following dataframe:")
    print(df.to_string(index=False))
    print("Could not find the required column '" + col + "'")
    quit()


def check_numeric(df, col):
  try:
    pandas.to_numeric(df[col])
  except ValueError:
    mask = pandas.to_numeric(df[col], errors='coerce').isna()
    print(df.loc[mask].to_string(index=False))
    print_error("could not interpret the above value(s) in the column '" +
                col + "'  as number(s)")
    quit()


def check_valid(df, col, valid, reason):
  df_invalid = df.loc[~df[col].isin(valid)]
  if (len(df_invalid.index)):
    print(df_invalid.to_string(index=False))
    print_error("above value(s) in the column '" + col +
                "' are invalid because " + reason)
    quit()


def remove_inchi_layers(inchi, delimiter_prefixes):
  for prefix in delimiter_prefixes:
    inchi = re.sub("/" + prefix + "[^/]*", "", inchi)
  return inchi


def remove_stereolayer(inchi):
  return remove_inchi_layers(inchi, ["b", "t", "m", "s"])


def check_stereochemistry(df):
  df_stereo = df.loc[df['stereochemistry'] == 0].copy()
  table = str.maketrans(dict.fromkeys('/\\-+'))
  df_stereo['smiles'] = df_stereo['smiles'].apply(
      lambda smiles: smiles.translate(table))
  df_stereo['inchi'] = df_stereo['inchi'].apply(remove_stereolayer)
  error_fragment1 = "the above printed "
  error_fragment2 = " remain unique after removing stereochemistry information."
  question = "\n\nIs it essential in these cases to specify RDKit incompatible"
  question += " information about stereochemistry? "
  question += "If so, the script needs to be revised."
  if (len(df_stereo.index) and df_stereo['inchi'].is_unique):
    print("\n")
    uni_loc = df_stereo.drop_duplicates(subset='inchi', keep=False).index
    print(df.loc[uni_loc].to_string(index=False))
    print_error(error_fragment1 + "InChI(s)" + error_fragment2 + question)
    quit()
  if (len(df_stereo.index) and df_stereo['smiles'].is_unique):
    print("\n")
    uni_loc = df_stereo.drop_duplicates(subset='smiles', keep=False).index
    print(df.loc[uni_loc].to_string(index=False))
    print_error(error_fragment1 + "smiles" + error_fragment2 + question)
    quit()

  species_name2inchi = df_stereo.set_index('model_name')['inchi'].to_dict()
  species_name2smiles = df_stereo.set_index('model_name')['smiles'].to_dict()
  return species_name2inchi, species_name2smiles


def check_nonzero_zero_multiplicity(df, valid_multiplicity):
  valid_multiplicity_wo_zero = [m for m in valid_multiplicity if m != 0]
  for index, row in df.iterrows():
    df_cur = df.loc[(df['inchi'] == row['inchi'])
                    & (df['excited'] == row['excited'])]
    if (len(df_cur.index) > 1):
      df_invalid = df_cur.loc[~df_cur['multiplicity'].
                              isin(valid_multiplicity_wo_zero)]
      if (len(df_invalid.index)):
        print(df_cur.to_string(index=False))
        reason = "species multiplicities must not be 0 and non-zero for otherwise identical species"
        print_error("above value(s) in the column '" + 'multiplicity' +
                    "' are invalid because " + reason)
        quit()


def read_csv_species_dict(csv_filename, species_list, fatal_error_only,
                          show_lumped):
  """ returns two dicts dict mapping species names to identifiers """
  check_file = Path(csv_filename)
  required_columns = [
      "model_name", "inchi", "smiles", "excited", "multiplicity", "lumped",
      "stereochemistry"
  ]
  if check_file.is_file():
    print("reading '" + csv_filename + "'")
    df = pandas.read_csv(csv_filename)
    df.insert(0, 'line', range(2, 2 + len(df)))
    for col in required_columns:
      check_required_column(df, col)

    if (df.isnull().values.any()):
      print(df[df.isnull().any(axis=1)].to_string(index=False))
      print_error("the above lines of the species dictionary '" +
                  csv_filename + "' must not contain nulls")
      quit()

    check_numeric(df, 'excited')
    check_numeric(df, 'multiplicity')
    check_numeric(df, 'lumped')
    check_numeric(df, 'stereochemistry')

    valid_excited = [0, 1]
    check_valid(df, 'excited', valid_excited,
                "a species is either excited or not")
    valid_multiplicity = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    check_valid(
        df, 'multiplicity', valid_multiplicity,
        "the multiplicity must be 0 or 2S+1 (S is the total spin quantum number)"
    )
    valid_lumped = [0, 1]
    check_valid(df, 'lumped', valid_lumped,
                "species are either lumped (1) or not (0)")
    valid_stereochemistry = [0, 1]
    check_valid(
        df, 'stereochemistry', valid_stereochemistry,
        "the stereochemistry information is either RDKit-compatible (1) or not (0)"
    )
    if (fatal_error_only):
      print("#warning: only fatal errors are considered")
      df = df[df['model_name'].isin(species_list)]

    any_error = False
    uniqueness_required = ['inchi', 'excited', 'multiplicity']
    if (not df.set_index(uniqueness_required).index.is_unique):
      print_error("duplicate species. The composite identifier " +
                  str(uniqueness_required) + " is non-unique for:")
      print_df_duplicates(df, uniqueness_required)
      any_error = True

    if (not fatal_error_only):
      check_nonzero_zero_multiplicity(df, valid_multiplicity)

    if (any_error):
      print("\nYou need to fix the above errors in '" + csv_filename + "'")
      quit()

    count = 0
    # example element:
    # canonical2data[('InChI=1S/He', 0, 0)]['smiles'] = '[He]'
    canonical2data = {}
    # .items() iterates over (column name, Series) pairs.
    for column, val in df.set_index(uniqueness_required).to_dict().items():
      for primary_key in val:
        if not primary_key in canonical2data:
          canonical2data[primary_key] = {}
        canonical2data[primary_key][column] = val[primary_key]

    species_name2inchi_wo_stereo, species_name2smiles_wo_stereo = check_stereochemistry(
        df)

    species_name2inchi_lumped = {}
    species_name2smiles_lumped = {}
    if (not df['model_name'].is_unique):
      df_duplicates = df[df.duplicated(subset='model_name', keep=False)]
      check_valid(
          df_duplicates, 'lumped', [1],
          "species names that appear multiple times must be marked as lumped")
      check_valid(df_duplicates, 'excited', [0],
                  "excited species are not expected to be lumped")
      check_valid(
          df_duplicates, 'stereochemistry', [1],
          "stereochemistry must be RDKit-compatible for lumped species")
      if (show_lumped):
        print("\n\nThe following species are lumped:")
        print_df_duplicates(df, ['model_name'])
      species_name2inchi_lumped = df_duplicates.set_index(
          'model_name')['inchi'].to_dict()
      species_name2smiles_lumped = df_duplicates.set_index(
          'model_name')['smiles'].to_dict()
      for name in species_name2inchi_lumped:
        df_lumped_species = df_duplicates.loc[df['model_name'] == name]
        species_name2inchi_lumped[name] = df_lumped_species['inchi'].to_list()
        species_name2smiles_lumped[name] = df_lumped_species['smiles'].to_list(
        )
        #print(species_name2smiles_lumped[name])

    species_name2inchi = df.set_index('model_name')['inchi'].to_dict()
    species_name2smiles = df.set_index('model_name')['smiles'].to_dict()
    for name in species_name2inchi:
      if (name in species_name2inchi_lumped):
        species_name2inchi[name] = species_name2inchi_lumped[name]
      else:
        species_name2inchi[name] = [species_name2inchi[name]]

    for name in species_name2smiles:
      if (name in species_name2smiles_lumped):
        species_name2smiles[name] = species_name2smiles_lumped[name]
      else:
        species_name2smiles[name] = [species_name2smiles[name]]

    return species_name2inchi, species_name2smiles, species_name2inchi_wo_stereo, species_name2smiles_wo_stereo, canonical2data
  else:
    print("csv species dictionary file '" + csv_filename + "' not found")
    quit()


def rdkit_error(what, identifier, species):
  print("RDKit cannot parse " + what + " '" + identifier + "' for '" +
        species + "'")


def rdkit_check(species_name2inchi, species_name2smiles,
                species_name2inchi_wo_stereo, species_name2smiles_wo_stereo):
  if (len(species_name2inchi) != len(species_name2smiles)):
    print("len(species_name2inchi) =", len(species_name2inchi))
    print("len(species_name2smiles) =", len(species_name2smiles))
    raise Exception("len(species_name2inchi) != len(species_name2smiles)")
  if (len(species_name2inchi_wo_stereo) != len(species_name2smiles_wo_stereo)):
    print("len(species_name2inchi_wo_stereo) =",
          len(species_name2inchi_wo_stereo))
    print("len(species_name2smiles_wo_stereo) =",
          len(species_name2smiles_wo_stereo))
    raise Exception(
        "len(species_name2inchi_wo_stereo) != len(species_name2smiles_wo_stereo)"
    )

  error = False

  for species in species_name2inchi_wo_stereo:
    if (not Chem.MolFromInchi(species_name2inchi_wo_stereo[species])):
      rdkit_error("InChI (w/o stereochemistry)", species_name2inchi[species],
                  species)
      error = True

  for species in species_name2smiles_wo_stereo:
    if (not Chem.MolFromSmiles(species_name2smiles_wo_stereo[species])):
      rdkit_error("SMILES (w/o stereochemistry)", species_name2smiles[species],
                  species)
      error = True

  if (error):
    print("\nFix the above InChI/SMILES error(s).")
    quit()

  for species in species_name2inchi_wo_stereo:
    if (check_inchi_smiles_consistency(
        species_name2inchi_wo_stereo[species],
        species_name2smiles_wo_stereo[species],
        species,
        simplification=" (w/o stereochemistry)")):
      error = True

  if (error):
    print("\nFix the above consistency error(s).")
    quit()
  else:
    print("InChI/SMILES check without stereochemistry successful")

  invalid_stereo_inchi = {}
  invalid_stereo_smiles = {}
  for species in species_name2inchi:
    for inchi in species_name2inchi[species]:
      if (not Chem.MolFromInchi(inchi)):
        if (species in species_name2inchi_wo_stereo):
          print("ignoring stereochemistry for '" + inchi + "'")
          invalid_stereo_inchi[species] = inchi
        else:
          rdkit_error("InChI", inchi, species)
          error = True

  for species in species_name2smiles:
    for smiles in species_name2smiles[species]:
      if (not Chem.MolFromSmiles(smiles)):
        if (species in species_name2smiles_wo_stereo):
          print("ignoring stereochemistry for '" + smiles + "'")
          invalid_stereo_smiles[species] = smiles
        else:
          rdkit_error("SMILES", smiles, species)
          error = True

  if (error):
    print("\nFix the above InChI/SMILES error(s).")
    quit()

  for species in species_name2inchi:
    if (len(species_name2inchi[species]) != len(species_name2smiles[species])):
      print("InChIs:", species_name2inchi[species])
      print("SMILES:", species_name2smiles[species])
      raise Exception("len(species_name2inchi[" + species +
                      "]) == len(species_name2smiles[" + species +
                      "]) is required. Something is wrong in the above lists.")

    for i in range(len(species_name2inchi[species])):
      if (not species in invalid_stereo_inchi
          and not species in invalid_stereo_smiles
          and check_inchi_smiles_consistency(species_name2inchi[species][i],
                                             species_name2smiles[species][i],
                                             species)):
        error = True
  if (error):
    print("\nFix the above consistency error(s).")
    quit()
  else:
    print("InChI/SMILES check with stereochemistry successful")

  return invalid_stereo_inchi, invalid_stereo_smiles


def make_species_dictionary(csv_filename,
                            species_list,
                            fatal_error_only=False,
                            show_lumped=False):
  species_name2inchi, species_name2smiles, species_name2inchi_wo_stereo, species_name2smiles_wo_stereo, canonical2data = read_csv_species_dict(
      csv_filename, species_list, fatal_error_only, show_lumped)

  invalid_stereo_inchi, invalid_stereo_smiles = rdkit_check(
      species_name2inchi, species_name2smiles, species_name2inchi_wo_stereo,
      species_name2smiles_wo_stereo)

  fallback_data = FallBack(invalid_stereo_inchi, invalid_stereo_smiles,
                           species_name2inchi_wo_stereo,
                           species_name2smiles_wo_stereo)
  return SpeciesDict(species_name2inchi, species_name2smiles, canonical2data,
                     fallback_data)


def make_species_dictionary_itv(thermochemistry_filename,
                                species_list,
                                fatal_error_only=False,
                                show_lumped=False,
                                c3_species_dict=None):
  """
  This routine reads the species dictionary from the thermochemistry 
  file and matches it against the c3_species_dict. Reading from the 
  species dictionary from the thermochemistry file is deprecated as 
  of December 2024. This routine was only written to convert the old 
  format the more detailed new one.
  """
  dry_run = False
  sc = DefaultClassification()
  lines = get_thermo_lines(thermochemistry_filename)
  sc, count_species, species_name2inchi, species_name2smiles = itv_mode_check(
      dry_run, species_list, lines, sc)

  itv_canonical2data = {}
  itv_not_found = {}
  for name in species_name2inchi:
    itv_not_found[name] = 1

  if (c3_species_dict):
    c3_unique_inchis = {}
    c3_primary_keys = c3_species_dict.get_valid_primary_keys()
    c3_species_name2primary_keys = c3_species_dict.get_species_name2primary_keys(
    )
    c3_unique_inchi_aux = {}
    for primary_key in c3_primary_keys:
      if (not primary_key[0] in c3_unique_inchi_aux):
        c3_unique_inchi_aux[primary_key[0]] = 1
      else:
        c3_unique_inchi_aux[primary_key[0]] += 1
    for inchi in c3_unique_inchi_aux:
      if (c3_unique_inchi_aux[inchi] == 1):
        c3_unique_inchis[inchi] = 1

    for itv_name in species_name2inchi:
      for itv_inchi in species_name2inchi[itv_name]:
        for c3_name in c3_species_name2primary_keys:
          for c3_primary_key in c3_species_name2primary_keys[c3_name]:
            c3_inchi = c3_primary_key[0]
            c3_excited = c3_primary_key[1]
            c3_multiplicity = c3_primary_key[2]

            if (c3_inchi == itv_inchi
                and (c3_multiplicity == 0 or itv_inchi in c3_unique_inchis)
                and c3_excited == 0 and len(species_name2inchi[itv_name]) == 1
                and len(c3_species_name2primary_keys[c3_name]) == 1):
              itv_canonical2data[(c3_inchi, c3_excited, c3_multiplicity)] = {}
              itv_canonical2data[(
                  c3_inchi, c3_excited,
                  c3_multiplicity)] = c3_primary_keys[c3_primary_key]
              itv_canonical2data[(c3_inchi, c3_excited,
                                  c3_multiplicity)]['model_name'] = itv_name
              itv_not_found[itv_name] = 0
            elif (c3_inchi == itv_inchi
                  and len(species_name2inchi[itv_name]) == 1
                  and len(c3_species_name2primary_keys[c3_name]) > 1):
              c3_multiplicities = {}
              c3_excited_vals = {}
              c3_multiplicity = 0
              c3_excited = 0
              for c3_primary_key in c3_species_name2primary_keys[c3_name]:
                c3_multiplicity = c3_primary_key[2]
                c3_excited = c3_primary_key[1]
                c3_multiplicities[c3_primary_key[2]] = 1
                c3_excited_vals[c3_primary_key[1]] = 1

              if (len(c3_excited_vals) == 1 and len(c3_multiplicities) == 1):
                primary_key = (itv_inchi, c3_multiplicity, c3_excited)
                itv_canonical2data[primary_key] = {}
                itv_canonical2data[primary_key]['inchi'] = itv_inchi
                itv_canonical2data[primary_key]['excited'] = primary_key[1]
                itv_canonical2data[primary_key]['multiplicty'] = primary_key[2]
                itv_canonical2data[primary_key]['model_name'] = itv_name
                if (len(species_name2inchi[itv_name]) == 1):
                  itv_canonical2data[primary_key]['lumped'] = 0
                else:
                  itv_canonical2data[primary_key]['lumped'] = 1
                itv_canonical2data[primary_key]['stereochemistry'] = 1
                smiles_lst = list(species_name2smiles[itv_name].keys())
                itv_canonical2data[primary_key]['smiles'] = smiles_lst[0]
                itv_not_found[itv_name] = 0
              else:
                print(itv_name, c3_name, "(c3 lumped)")
            elif (c3_inchi == itv_inchi
                  and len(species_name2inchi[itv_name]) > 1):
              c3_multiplicities = {}
              c3_excited_vals = {}
              c3_multiplicity = 0
              c3_excited = 0
              for c3_primary_key in c3_species_name2primary_keys[c3_name]:
                c3_multiplicity = c3_primary_key[2]
                c3_excited = c3_primary_key[1]
                c3_multiplicities[c3_primary_key[2]] = 1
                c3_excited_vals[c3_primary_key[1]] = 1

              if (len(c3_excited_vals) == 1 and len(c3_multiplicities) == 1):
                primary_key = (itv_inchi, c3_multiplicity, c3_excited)
                itv_canonical2data[primary_key] = {}
                itv_canonical2data[primary_key]['inchi'] = itv_inchi
                itv_canonical2data[primary_key]['excited'] = primary_key[1]
                itv_canonical2data[primary_key]['multiplicty'] = primary_key[2]
                itv_canonical2data[primary_key]['model_name'] = itv_name
                if (len(species_name2inchi[itv_name]) == 1):
                  itv_canonical2data[primary_key]['lumped'] = 0
                else:
                  itv_canonical2data[primary_key]['lumped'] = 1
                itv_canonical2data[primary_key]['stereochemistry'] = 1
                smiles_lst = list(species_name2smiles[itv_name].keys())
                smiles_count = 0
                for tmp_inchi in species_name2inchi[itv_name]:
                  if (tmp_inchi == itv_inchi):
                    break
                  smiles_count += 1
                itv_canonical2data[primary_key]['smiles'] = smiles_lst[
                    smiles_count]
                itv_not_found[itv_name] = 0
            elif (c3_inchi == itv_inchi
                  and len(species_name2inchi[itv_name]) == 1
                  and len(c3_species_name2primary_keys[c3_name]) == 1):
              itv2c3name = {
                  'S-CH2': 'CH2(S)',
                  'T-CH2': 'CH2',
                  'S-C3H2': 'C3H2(S)',
              }
              # ignored because Raymond knows we don't have this in the ITV model
              c3_ignore_rl = {
                  'CHV': 1,
                  'OHV': 1,
              }
              # ignore c3name for which Raymond did the matching manually
              for itv_tmp_name in itv2c3name:
                c3_ignore_rl[itv2c3name[itv_tmp_name]] = 1

              if (itv_name == c3_name
                  or (itv_name in itv2c3name
                      and itv2c3name[itv_name] == c3_name)):
                itv_canonical2data[(c3_inchi, c3_excited,
                                    c3_multiplicity)] = {}
                itv_canonical2data[(
                    c3_inchi, c3_excited,
                    c3_multiplicity)] = c3_primary_keys[c3_primary_key]
                itv_canonical2data[(c3_inchi, c3_excited,
                                    c3_multiplicity)]['model_name'] = itv_name
                itv_not_found[itv_name] = 0
              elif (c3_name in c3_ignore_rl or itv_name in itv2c3name):
                pass
              else:
                print(itv_name, c3_name, "excited or multiplicity specified")

    c3_inchi2name = {}
    for c3_name in c3_species_name2primary_keys:
      for c3_inchi in c3_species_name2primary_keys[c3_name]:
        c3_inchi2name[c3_inchi[0]] = c3_name

    print("c3 species dictionary covers " + str(len(c3_inchi2name)) +
          " inchis")
    count_no_overlap = {}
    for itv_name in species_name2inchi:
      overlap = False
      for itv_inchi in species_name2inchi[itv_name]:
        #print("checking " + itv_inchi)
        if (itv_inchi in c3_inchi2name):
          #print(itv_name + " overlaps with C3")
          overlap = True
          break

      if (not overlap):
        #print("adopted ITV definition for " + itv_name)
        itv_not_found[itv_name] = 0

        inchi_count = 0
        for itv_inchi in species_name2inchi[itv_name]:
          primary_key = (itv_inchi, 0, 0)
          itv_canonical2data[primary_key] = {}
          itv_canonical2data[primary_key]['inchi'] = itv_inchi
          itv_canonical2data[primary_key]['excited'] = primary_key[1]
          itv_canonical2data[primary_key]['multiplicty'] = primary_key[2]
          itv_canonical2data[primary_key]['model_name'] = itv_name
          if (len(species_name2inchi[itv_name]) == 1):
            itv_canonical2data[primary_key]['lumped'] = 0
          else:
            itv_canonical2data[primary_key]['lumped'] = 1
          itv_canonical2data[primary_key]['stereochemistry'] = 1
          smiles_lst = list(species_name2smiles[itv_name].keys())
          itv_canonical2data[primary_key]['smiles'] = smiles_lst[inchi_count]
          count_no_overlap[itv_name] = 1
          inchi_count += 1
    print(
        str(len(count_no_overlap)) +
        " ITV species have no overlap with C3Mech")

  print("Could not find the following species")
  for name in itv_not_found:
    if (itv_not_found[name]):
      print(name)
  itv_names_assigned = {}
  for primary_key in itv_canonical2data:
    itv_names_assigned[itv_canonical2data[primary_key]['model_name']] = 1
  print("Found " + str(len(itv_names_assigned)) + "/" +
        str(len(species_name2inchi)) +
        " entries for the new ITV species dictionary")
  itv_species_dict = SpeciesDict(species_name2inchi, species_name2smiles,
                                 itv_canonical2data, FallBack({}, {}, {}, {}))
  df = itv_species_dict.write_to_csv('ITV_species_dict.csv')
  return itv_species_dict


def make_species_classification(inp, sc, dry_run, species_list, itv_mode,
                                fatal_error_only, show_lumped, silent):
  print("process InChIs and SMILES formulas...")
  count_species = 0
  timestamp_species_dict = 0.0
  if (itv_mode):
    if (not silent):
      print("reading the thermochemistry file '" + inp.thermo_filename + "'")
    lines = get_thermo_lines(inp.thermo_filename)
    sc, count_species, species_name2inchi, species_name2smiles = itv_mode_check(
        dry_run, species_list, lines, sc)
    timestamp_species_dict = os.path.getmtime(inp.thermo_filename)
  else:
    timestamp_species_dict = os.path.getmtime(inp.csv_filename)
    spec_dict = make_species_dictionary(inp.csv_filename, species_list,
                                        fatal_error_only, show_lumped)
    # the thermochemistry file is optional in itv_mode
    if (inp.thermo_filename != ''):
      count = 0
      count_missing_inchi = 0
      print(
          "perform optional crosscheck with thermochemistry element compositions"
      )
      print("reading the thermochemistry file '" + inp.thermo_filename + "'")
      lines = get_thermo_lines(inp.thermo_filename)
      species_name2sum_formula, found_nasa = get_species_name2sum_formula(
          lines, species_list, inp.thermo_filename, False)
      for species in species_list:
        if (species in spec_dict.species_name2inchi
            and species in species_name2sum_formula):
          #print(species)
          for inchi in spec_dict.species_name2inchi[species]:
            check_inchi_composition_consistency(
                inchi, species_name2sum_formula[species], species)
          count += 1
        if (not species in spec_dict.species_name2inchi
            and species in species_name2sum_formula):
          count_missing_inchi += 1
      print(
          "successfully crosschecked " + str(count) + "/" +
          str(len(species_list)) +
          " InChI sum formulas against thermochemistry element compositions (skipped "
          + str(count_missing_inchi) + "/" + str(len(species_list) - count) +
          " due to missing InChIs)")
      #if(len(identifier2sum_formula) > 0):

      missing_inchi = {}
      for species in species_list:
        if (species in spec_dict.species_name2inchi):
          sc.add_species(dry_run, species,
                         spec_dict.species_name2inchi[species],
                         spec_dict.species_name2smiles[species],
                         spec_dict.fallback)
          count_species += 1
        else:
          missing_inchi[species] = 1

      if (len(missing_inchi) > 0):
        print("InChIs are missing for the following " +
              str(len(missing_inchi)) + " species:")
        for species in missing_inchi:
          print(species)
  return sc, count_species, timestamp_species_dict


def check_species_list(species_list, thermo_filename):
  all_found = True
  for s in species_list:
    if (species_list[s] == 0):
      if (all_found):
        print("Could not find the following species:")
      print(s)
      all_found = False

  if (not all_found):
    print(
        "WARNING: COULD NOT FIND AN INCHI FOR EVERY SPECIES. CHECK THE THERMOCHEMISTRY INPUT: '"
        + thermo_filename + "'\n")
  else:
    print_success()


def get_git_sha(path):
  try:
    repo = git.Repo(path, search_parent_directories=True)
    return repo.head.object.hexsha
  except git.exc.InvalidGitRepositoryError:
    return "no git repository found"


def copy_style_files(src_dir, dest_dir):
  copyfile(os.path.join(src_dir, "placeins.sty"),
           os.path.join(dest_dir, "placeins.sty"))
  copyfile(os.path.join(src_dir, "xurl.sty"),
           os.path.join(dest_dir, "xurl.sty"))
  copyfile(os.path.join(src_dir, "url.sty"), os.path.join(dest_dir, "url.sty"))
  copyfile(os.path.join(src_dir, "hyphenat.sty"),
           os.path.join(dest_dir, "hyphenat.sty"))
  copyfile(os.path.join(src_dir, "miscdoc.sty"),
           os.path.join(dest_dir, "miscdoc.sty"))


def set_RDKit_drawing_option():
  # configuration
  DrawingOptions.atomLabelFontSize = 75
  DrawingOptions.dotsPerAngstrom = 100
  DrawingOptions.bondLineWidth = 3.0


def get_tex_begin_and_end(title, authors, sha, tex_opts):
  species_dict_begin = [
      "\\let\\mypdfximage\\pdfximage\n",
      "\\def\\pdfximage{\\immediate\\mypdfximage}\n", "\n",
      "\\documentclass[a4paper,fontsize=7pt]{scrartcl}\n", "\n",
      "\\usepackage[english]{babel}\n", "\\usepackage{xurl}\n",
      "\\usepackage{lmodern}\n", "\\usepackage[none]{hyphenat}\n",
      "\\usepackage{graphicx}\n", "\\renewcommand{\\baselinestretch}{" +
      str(tex_opts.baselinestretch) + "}\n",
      "\\usepackage[font=small]{caption}\n", "\\usepackage{subcaption}\n",
      "\\usepackage[top=0.6in,bottom=0.6in,left=" + str(tex_opts.side_margin) +
      "in,right=" + str(tex_opts.side_margin) + "in]{geometry}\n",
      "\\usepackage[figurename=Fig.]{caption}\n",
      "\\usepackage[T1]{fontenc} % allows for the correct printing of _\n",
      "\\usepackage[section]{placeins}\n", "\n", "\\title{" + title + "}\n",
      "\\subtitle{model hash: " + sha + "}\n", "\\author{" + authors + "}\n",
      "\\date{\\vspace{-5ex}}\n",
      "%\\pdfsuppresswarningpagegroup=1 % works only with newer pdflatex versions\n",
      "\n", "\\pdfminorversion=7\n", "\\begin{document}\n", "\\maketitle\n"
      "Each figure shows a species with its name in the kinetic model,\n",
      "its International Chemical Identifier (InChI),\n",
      "and its simplified molecular-input line-entry system (SMILES) formula.\n"
  ]

  species_dict_end = ["\\end{document}\n"]

  return species_dict_begin, species_dict_end


def write_dict(title, authors, sha, classification, tex_opts,
               timestamp_species_dict, output_dir):
  species_dict_begin, species_dict_end = get_tex_begin_and_end(
      title, authors, sha, tex_opts)
  tex_file = open(os.path.join(output_dir, "species_dict.tex"), "w")
  for line in species_dict_begin:
    tex_file.write(line)

  classification.write_all(tex_file, output_dir, tex_opts,
                           timestamp_species_dict)

  for line in species_dict_end:
    tex_file.write(line)
  tex_file.close()


def get_species_list_submodel(yaml_filename, cmd_output_directory,
                              cmd_species_dictionary):
  submodel_files = make_submodelfiles_from_yaml(yaml_filename,
                                                cmd_output_directory,
                                                cmd_species_dictionary)
  species_list = write_species.make_species_list('', submodel_files.files)

  species_dict = {k.upper(): 0 for k in species_list}
  for s in species_dict:
    species_dict[s] = 0
  return species_dict, submodel_files


def get_species_list(kinetics_filename, silent=True):
  species_dict = {}
  print("Reading the chemkin kinetics file '" + kinetics_filename + "'")
  check_file = Path(kinetics_filename)
  lines_mech = []
  if check_file.is_file():
    kinetics_file = open(kinetics_filename, encoding='cp1252')
    lines_mech = kinetics_file.readlines()
    kinetics_file.close()
  else:
    print_not_found("chemkin kinetics", "file", kinetics_filename)
    quit()
  species_dict = make_species_from_kinetics_file(lines_mech)

  species_dict = {k.upper(): v for k, v in species_dict.items()}
  for s in species_dict:
    species_dict[s] = 0
  return species_dict


def write_species_dict(inp,
                       out,
                       dry_run,
                       itv_mode,
                       fatal_error_only,
                       show_lumped,
                       sc=DefaultClassification(),
                       silent=False):
  if (not silent):
    print("Constructing the species list")

  start = time.time()

  species_list = {}
  if (inp.yaml_filename != ''):
    # construct species list from submodule(s) specified in the
    # yaml file
    species_list, submodel_files = get_species_list_submodel(
        inp.yaml_filename, out.output_dir, inp.csv_filename)
    output_dir = submodel_files.output_directory
    inp.csv_filename = submodel_files.species_dictionary
  else:
    species_list = get_species_list(inp.mech_filename, silent)
    check_output_dir = Path(out.output_dir)
    if (not check_output_dir.is_dir()):
      print("output directory '" + out.output_dir + "' does not exist.")
      quit()
    output_dir = out.output_dir

  print("found " + str(len(species_list)) +
        " species from the selected kinetic (sub)model(s)")
  if (not silent):
    print_success()
  set_RDKit_drawing_option()

  if (not silent):
    print("Copying style files to the output directory '" + output_dir + "'")

  if (not dry_run):
    copy_style_files(inp.style_file_dir, output_dir)
  if (not silent):
    print_success()

  end = time.time()
  ## print("Time until make_species_classification:", end - start)

  start2 = time.time()
  sc, n_species_found, timestamp_species_dict = make_species_classification(
      inp, sc, dry_run, species_list, itv_mode, fatal_error_only, show_lumped,
      silent)
  end2 = time.time()
  ## print("Time required for make_species_classification:", end2 - start2)
  if (itv_mode):
    check_species_list(species_list, inp.thermo_filename)

  if (dry_run):
    print("skipping output generation")
    quit()
  print("generating species images and tex input...")
  start_write_dict = time.time()
  write_dict(out.title, out.authors, out.sha, sc, TexOptions(),
             timestamp_species_dict, output_dir)
  end_write_dict = time.time()
  ## print("Time required for output generation:",
  ##       end_write_dict - start_write_dict)
  print_success()

  print("compile species_dict.pdf with: 'cd \"" + output_dir +
        "\" && pdflatex species_dict.tex'")


if __name__ == "__main__":
  RDLogger.DisableLog('rdApp.*')

  preprocessor = os.path.join("..", "PREPROCESSOR")
  submechanisms = os.path.join("..", "SUBMECHANISMS")
  parser = argparse.ArgumentParser(
      description='Generates figures and latex inputs')
  parser.add_argument(
      '-t',
      '--thermo',
      help='chemkin thermochemistry file (may contain SMILES and InChIs)',
      default=os.path.join(submechanisms, "SOURCE-C3Mech.THERM"))
  parser.add_argument(
      '-k',
      '--kinetics',
      help=
      'chemkin kinetics file that can be used to read in the species. requires -y \'\' and will be ignored otherwise',
      default=os.path.join(preprocessor,
                           os.path.join("Input_Abstractions", "C3.dic")))
  parser.add_argument(
      '-y',
      '--yaml',
      help=
      'yaml file listing submodels considered in the species dictionary generation',
      default="submodels.yaml")
  parser.add_argument(
      '-o',
      '--output_dir',
      help=
      'output directory which will be used to generate the species dictionary',
      default="")
  parser.add_argument(
      '-c',
      '--csv',
      help=
      'csv file with the species defintions used to generate the species dictionary',
      default="")
  parser.add_argument('-i',
                      '--itv',
                      action='store_true',
                      help='the thermochemistry file is in ITV format')
  parser.add_argument('-l',
                      '--lumped',
                      help='show lumped species',
                      action='store_true')
  parser.add_argument(
      '-f',
      '--fatal_error_only',
      help=
      'only species currently used in the (sub)model(s) are considered in the checks',
      action='store_true')
  parser.add_argument('-n',
                      '--name',
                      help='name (title) of species dictionary',
                      default="C3Mech species dictionary")
  parser.add_argument(
      '-a',
      '--authors',
      help='list of authors',
      default=
      "Yuki Murakami, Luna Pratali Maffei, Raymond Langer, Sanket Girhe,\\\\\nScott W. Wagnon, Goutham Kukkadapu, Matteo Pelucchi, Kuiwen Zhang,\\\\\nMandhapati Raju, Tanusree Chatterjee, William J. Pitz, Tiziano Faravelli,\\\\\nHeinz Pitsch, Peter Kelly Senecal, Henry J. Curran"
  )
  parser.add_argument('-g',
                      '--git',
                      help='directory that contains the git repository',
                      default=".")
  parser.add_argument('-v',
                      '--verbose',
                      action='store_true',
                      help="enable verbose execution")
  parser.add_argument(
      '-d',
      '--dry_run',
      action='store_true',
      help="perform only checks without output generation (runs faster)")
  args = vars(parser.parse_args())

  inp = Input(args['csv'], args['yaml'], args['kinetics'], args['thermo'],
              "style_files_spec_dict")

  check_git = Path(args['git'])
  if (not check_git.is_dir()):
    raise Exception("The directory '" + args['git'] +
                    "' containing the git repository does not exist.")
  out = Output(args['output_dir'], args['name'], args['authors'],
               get_git_sha(args['git']))
  silent = not args['verbose']

  write_species_dict(inp,
                     out,
                     dry_run=args["dry_run"],
                     itv_mode=args['itv'],
                     fatal_error_only=args['fatal_error_only'],
                     show_lumped=args['lumped'],
                     silent=silent)
