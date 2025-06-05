#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import datetime
import yaml

from shutil import copyfile

VERSION = "4.0"


def make_species_list(directory, files_list):
  # start reading from first "kinetic_module" you find
  species_list = []
  for file in files_list:
    print("reading '" + os.path.basename(file) + "'")
    with open(os.path.join(directory, file), encoding='cp1252') as ckiblock:
      check_rxns = 0
      for line in ckiblock:
        if '!\KINETICS_MODULE' in line:
          check_rxns = 1
        lineunc = line.split('!')[0].strip()
        if len(lineunc) > 0:
          if check_rxns == 1:
            # check different types of separators hierarchically
            sep = ''
            if '<=>' in lineunc:
              sep = '<=>'
            elif '=>' in lineunc:
              sep = '=>'
            elif '=' in lineunc:
              sep = '='

            if sep != '':
              rcts_rough, prds_rough = lineunc.split(sep)
              rcts = rcts_rough.strip().split('(+M)')[0].split('+')
              # hopefully reactants is enough
              prds = prds_rough.strip().split('(+M)')[0].split('+')
              spcs = rcts + prds
              #print(spcs)
              # species may contain spaces for how they were derived: remove
              for spc in spcs:
                spc_clean = spc.strip()
                if spc_clean[0] == '2':  # few reactants written as 2A = B
                  spc_clean = spc_clean[1:]
                # check it is actually a species
                if len(spc_clean.split()) > 1:
                  spc_clean = spc_clean.split(
                  )[0]  # only possible species is the first one - in case it is a float, it won't be included

                if spc_clean[0] in [
                    '0', '1', '3', '4', '5', '6', '7', '8', '9', '.'
                ]:
                  continue
                try:
                  float(spc_clean)  # try to convert to float after splitting
                except ValueError:
                  # if unclosed bracket: residue from a (+M) written with spaces perhaps
                  if spc_clean[-1] == '(':
                    spc_clean = spc_clean.split('(')[0]
                  elif spc_clean[-1] == ')' and '(' not in spc_clean:
                    spc_clean = spc_clean.split(')')[0]

                  if spc_clean not in species_list and spc_clean != 'M':
                    # print(spc_clean)
                    species_list.append(spc_clean)

  mandatory = ['HE', 'N2', 'AR', 'C2H6', 'CO2', 'CH4']
  for m in mandatory:
    if m not in species_list:
      species_list.insert(0, m)

  return species_list


def write_species_list(species_list, filename):
  # write spc str
  commentstr = '\n!+++++++++++++++++++++ '
  speciesallstr = '\n!\SPECIES_MODULE: ALL '
  speciesendstr = '\n!\END_SPECIES_MODULE: ALL '
  allspcstr = '\n'
  count = 0
  spcsperline = 3
  for spc in species_list:
    count += 1
    allspcstr += spc
    allspcstr += ' ' * (30 - len(spc))
    if count % spcsperline == 0:
      allspcstr += '\n'

  total_str = commentstr + speciesallstr + commentstr + allspcstr + commentstr + speciesendstr + commentstr
  if filename:
    spcfile = open(filename, "w")
    spcfile.writelines(total_str)
    spcfile.close()
  return total_str
  #print(species_list)


def clean_therm(therm_file, therm_output, species_list, datetime):
  therm = open(therm_file, 'r', encoding="utf-8")
  thermlines = therm.readlines()
  therm.close()

  print("writing '" + therm_output + "'")
  with open(therm_output, 'w') as thermfile:
    thermfile.write('! Thermodynamic data for C3MechV' + VERSION +
                    ', generated on {0:%d/%m/%Y %H:%M:%S}.\n'.format(datetime))
    thermfile.write('THERMO ALL\n')
    thermfile.write('300.   1000.   5000.\n')

    species_dict = {s.upper(): 0 for s in species_list}

    for i in range(len(thermlines)):
      wrk_line = thermlines[i].split("!")[0]
      wrk_line = wrk_line.rstrip()
      if len(wrk_line) >= 80 and wrk_line[79] == '1':
        sp_name = wrk_line[0:25].split()[0].upper()
        if sp_name in species_dict and not species_dict[
            sp_name] and i + 3 < len(thermlines):
          thermfile.write(thermlines[i])
          thermfile.write(thermlines[i + 1])
          thermfile.write(thermlines[i + 2])
          thermfile.write(thermlines[i + 3])
          species_dict[sp_name] = 1

    for s, v in species_dict.items():
      if v == 0:
        print('{0} : thermochemistry data not found...'.format(s))

    thermfile.write('END\n')


def clean_tran(tran_file, tran_output, species_list, datetime):
  tran = open(tran_file, 'r', encoding="utf-8")
  tranlines = tran.readlines()
  tran.close()

  species_dict = {s.upper(): 0 for s in species_list}
  print("writing '" + tran_output + "'")
  with open(tran_output, 'w') as tranfile:
    tranfile.write('! Transport data for C3MechV' + VERSION +
                   ', generated on {0:%d/%m/%Y %H:%M:%S}.\n'.format(datetime))
    for i in range(len(tranlines)):
      wrk_line = tranlines[i].split("!")[0]
      wrk_line = wrk_line.rstrip()
      if wrk_line == '':
        continue
      sp_name = wrk_line.split()[0].upper()
      if sp_name in species_dict and not species_dict[sp_name]:
        species_dict[sp_name] = 1
        tran_component = tranlines[i].split()
        if ('!' in tranlines[i]):
          tranfile.write(
              '{0:18}    {1}    {2:8.2f}    {3:6.2f}    {4:6.2f}    {5:6.2f}    {6:6.2f}    {7}\n'
              .format(tran_component[0], tran_component[1],
                      float(tran_component[2]), float(tran_component[3]),
                      float(tran_component[4]), float(tran_component[5]),
                      float(tran_component[6]), ' '.join(tran_component[7:])))
        else:
          tranfile.write(
              '{0:18}    {1}    {2:8.2f}    {3:6.2f}    {4:6.2f}    {5:6.2f}    {6:6.2f}\n'
              .format(tran_component[0], tran_component[1],
                      float(tran_component[2]), float(tran_component[3]),
                      float(tran_component[4]), float(tran_component[5]),
                      float(tran_component[6])))
        #tranfile.write('{0}\n'.format(tranlines[line].rstrip()))

    for s, v in species_dict.items():
      if v == 0:
        print('{0} : transport data not found...'.format(s))


def get_submech_dir():
  return get_absolute_path(os.path.join("..", "SUBMECHANISMS"))


def print_not_found(what, dir_or_file, path):
  print(what + " " + dir_or_file + " '" + path + "' does not exist")


def read_yaml_input(filename):
  """ returns a SubModelFiles """
  yaml_input_options = None
  with open(filename) as inp:
    try:
      yaml_input_options = yaml.safe_load(inp)
    except ValueError as e:
      util.print_error("invalid syntax in yaml input '" + filename + "'")
  return yaml_input_options


def make_submodelfiles_from_yaml(filename, directory):
  """ 
      This routine returns a checked SubModelFiles object if successful.
      Note: the routine may quit the script. 
  """
  if (not os.path.isfile(filename)):
    print_not_found("yaml input", "file", filename)
    quit()
  print("reading yaml file \'" + os.path.basename(filename) + "'")
  submodels = read_yaml_input(filename)
  submodels.insert_model_path(directory)

  if (not submodels.check()):
    print("error: invalid input")
    quit()

  return submodels


def get_species_list_submodel(yaml_filename):
  submodel_files = make_submodelfiles_from_yaml(yaml_filename,
                                                get_submech_dir())

  species_list = make_species_list('', submodel_files.get_files())

  species_dict = {k.upper(): 0 for k in species_list}
  for s in species_dict:
    species_dict[s] = 0
  return species_dict, submodel_files


class SubModelFiles(yaml.YAMLObject):
  yaml_loader = yaml.SafeLoader
  yaml_tag = u'!SubModelFiles'

  def __init__(self, core, submechs):
    self.core = copy.copy(core)
    self.submechs = copy.copy(submechs)
    if (submechs is None):
      self.submechs = []

  def check(self):
    ok = True
    if (not os.path.isfile(self.core)):
      print_not_found("coremech", "file", self.core)
      ok = False
    for filename in self.submechs:
      if (not os.path.isfile(filename)):
        print_not_found("submech", "file", filename)
        ok = False
    return ok

  def insert_model_path(self, directory):
    self.core = os.path.join(directory, self.core)
    self.submechs = [os.path.join(directory, f) for f in self.submechs]

  def get_files(self):
    return [self.core] + self.submechs


def get_absolute_path(relative_path):
  # Get the directory of the currently executed script
  script_dir = os.path.dirname(os.path.abspath(__file__))

  # Merge the script directory with the relative path
  absolute_path = os.path.join(script_dir, relative_path)

  return absolute_path


def insert_species_list(species_list, file_path, new_lines):
  with open(file_path, 'r') as file:
    lines = file.readlines()

  species_found = False
  end_found = False

  for line in lines:
    # Check for SPECIES keyword (case-insensitive)
    if line.strip().lower().startswith('species'):
      species_found = True
      new_lines.append(line)  # Copy the SPECIES line
      continue

    # If we have found SPECIES but not END yet, handle lines to replace
    if species_found and not end_found:
      if line.strip().lower().startswith('end'):
        end_found = True
        new_lines.append(species_list + "\n")  # Insert species_list before END
        new_lines.append(line)  # Copy the END line
        continue

      # If we're in between SPECIES and END, skip these lines (to be replaced)
      continue

    # If we haven't found SPECIES or after END, copy other lines as is
    new_lines.append(line)

  return new_lines


def merge_models(species_list, submodel_files, output_filename, datetime):
  species_section_str = write_species_list(species_list, '')
  new_lines = []
  new_lines.append('! Kinetics data for C3MechV' + VERSION +
                   ', generated on {0:%d/%m/%Y %H:%M:%S}.\n'.format(datetime) +
                   '\n')
  new_lines.append(
      '! The following sub-models were considered in this file:\n')
  new_lines.append('! - ' + os.path.basename(submodel_files.core) + '\n')
  for file_path in submodel_files.submechs:
    new_lines.append('! - ' + os.path.basename(file_path) + '\n')
  new_lines = insert_species_list(species_section_str, submodel_files.core,
                                  new_lines)
  for file_path in submodel_files.submechs:
    with open(file_path, 'r') as file:
      lines = file.readlines()
      new_lines += lines
  new_lines.append("END\n\n")
  print("writing '" + output_filename + "'")
  with open(output_filename, 'w', encoding="utf-8") as file:
    file.writelines(new_lines)


if __name__ == "__main__":

  print("read input...")
  DATETIME = datetime.datetime.now()
  species_list, submodel_files = get_species_list_submodel("submodels.yaml")
  print("\ngenerate output...")
  merge_models(species_list, submodel_files,
               os.path.join("output", "C3Mech.CKI"), DATETIME)
  therm_input = os.path.join(get_submech_dir(), "SOURCE-C3Mech.THERM")
  clean_therm(therm_input, os.path.join("output", "C3Mech.THERM"),
              species_list, DATETIME)
  trans_input = os.path.join(get_submech_dir(), "SOURCE-C3Mech.TRAN")
  clean_tran(trans_input, os.path.join("output", "C3Mech.TRAN"), species_list,
             DATETIME)
