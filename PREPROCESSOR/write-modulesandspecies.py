#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, argparse, datetime 
from shutil import copyfile

def copy_llnl(dest_dir):
  src_dir = os.path.join("..", "SUBMECHANISMS", "LLNL")
  copyfile(os.path.join(src_dir, "LLNL_BLOCK.CKI"), os.path.join(dest_dir, "LLNL_BLOCK.CKI"))
  copyfile(os.path.join(src_dir, "LLNL_BLOCK_HT.CKI"), os.path.join(dest_dir, "LLNL_BLOCK_HT.CKI"))
  copyfile(os.path.join(src_dir, "LLNL_BLOCK_DMC_EC.CKI"), os.path.join(dest_dir, "LLNL_BLOCK_DMC_EC.CKI"))

def copy_pah(dest_dir):
  src_dir = os.path.join("..", "SUBMECHANISMS", "ITV_POLIMI")
  copyfile(os.path.join(src_dir, "PAH_BLOCK.CKI"), os.path.join(dest_dir, "PAH_BLOCK.CKI"))

def copy_nuig(dest_dir):
  src_dir = os.path.join("..", "SUBMECHANISMS", "NUIG", "LT-HT")
  copyfile(os.path.join(src_dir, "NUIG_C-N_LT-HT.MECH"), os.path.join(dest_dir,   "NUIG_C-N_LT-HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C0_LT-HT.MECH"), os.path.join(dest_dir,    "NUIG_C0_LT-HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C1-C2_LT-HT.MECH"), os.path.join(dest_dir, "NUIG_C1-C2_LT-HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C3-C4_LT-HT.MECH"), os.path.join(dest_dir, "NUIG_C3-C4_LT-HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C5_LT-HT.MECH"), os.path.join(dest_dir,    "NUIG_C5_LT-HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C5cy_LT-HT.MECH"), os.path.join(dest_dir,  "NUIG_C5cy_LT-HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C6_LT-HT.MECH"), os.path.join(dest_dir,    "NUIG_C6_LT-HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C6cy_LT-HT.MECH"), os.path.join(dest_dir,  "NUIG_C6cy_LT-HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C7_LT-HT.MECH"), os.path.join(dest_dir,    "NUIG_C7_LT-HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_MTBE_LT-HT.MECH"), os.path.join(dest_dir,  "NUIG_MTBE_LT-HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_N_LT-HT.MECH"), os.path.join(dest_dir,     "NUIG_N_LT-HT.MECH"))

def copy_nuig_ht(dest_dir):
  src_dir = os.path.join("..", "SUBMECHANISMS", "NUIG", "HT")
  copyfile(os.path.join(src_dir, "NUIG_C-N_HT.MECH"), os.path.join(dest_dir,   "NUIG_C-N_HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C3-C4_HT.MECH"), os.path.join(dest_dir, "NUIG_C3-C4_HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C5_HT.MECH"), os.path.join(dest_dir,    "NUIG_C5_HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C5cy_HT.MECH"), os.path.join(dest_dir,  "NUIG_C5cy_HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C6_HT.MECH"), os.path.join(dest_dir,    "NUIG_C6_HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C6cy_HT.MECH"), os.path.join(dest_dir,  "NUIG_C6cy_HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_C7_HT.MECH"), os.path.join(dest_dir,    "NUIG_C7_HT.MECH"))
  copyfile(os.path.join(src_dir, "NUIG_N_HT.MECH"), os.path.join(dest_dir,     "NUIG_N_HT.MECH"))

def get_files_list():
    files_list = []
    
    with open('input_tot.dic') as inputtot:
        checksubmechs = 0
        for line in inputtot:
            # del commented line
            lineunc = line.split('//')[0]
            if len(lineunc) > 0:
                # search for coremech
                if 'CoreMechanism' in lineunc:
                    file = line.split('CoreMechanism')[1].strip().split(';')[0]
                    files_list.append(file)
                if checksubmechs == 1:
                    files = lineunc.split('}')[0].split(';')[0].strip().split()
                    if '' in files:
                        files.remove('')
                    if 'AAA_allspecies.txt' in files:
                        files.remove('AAA_allspecies.txt')
                        
                    if len(files) > 0:
                        files_list += files
                        
                if 'SubMechanisms' in lineunc:
                    checksubmechs = 1
                    file = lineunc.split('@SubMechanisms')[1].strip().split(';')[0].strip()
                    if len(file) > 0:
                        files_list.append(file)
    print(files_list)
    return files_list
    inputtot.close()                

def make_species_list(directory, files_list):
    # start reading from first "kinetic_module" you find
    species_list = []
    for file in files_list:
        print("reading '" + file + "'")
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
                                if spc_clean[0] == '2': # few reactants written as 2A = B
                                    spc_clean = spc_clean[1:]
                                # check it is actually a species
                                if len(spc_clean.split()) > 1:
                                    spc_clean = spc_clean.split()[0] # only possible species is the first one - in case it is a float, it won't be included
                                    
                                if spc_clean[0] in ['0', '1', '3', '4', '5', '6', '7', '8', '9', '.']:
                                    continue
                                try:
                                    float(spc_clean) # try to convert to float after splitting
                                except ValueError:
                                    # if unclosed bracket: residue from a (+M) written with spaces perhaps
                                    if spc_clean[-1] == '(':
                                        spc_clean = spc_clean.split('(')[0]
                                    elif spc_clean[-1] == ')' and '(' not in spc_clean: 
                                        spc_clean = spc_clean.split(')')[0]
                                                                       
                                    if spc_clean not in species_list and spc_clean != 'M':
                                        # print(spc_clean)
                                        species_list.append(spc_clean)
    
    
    print("done")
    mandatory = ['HE', 'N2', 'AR', 'C2H6', 'CO2', 'CH4']
    for m in mandatory:
        if m not in species_list:
            species_list.insert(0, m)

    return species_list

def write_species_list(species_list):
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
        allspcstr += ' '*(30 - len(spc))
        if count % spcsperline == 0:
            allspcstr += '\n'

    total_str = commentstr + speciesallstr + commentstr + allspcstr + commentstr + speciesendstr + commentstr
    spcfile = open(os.path.join('AAA_Input_Kinetics', 'AAA_allspecies.txt'), "w")
    spcfile.writelines(total_str)
    spcfile.close()
    #print(species_list)

def clean_therm(therm_file, therm_output, species_list, datetime):
    therm = open(therm_file, 'r', encoding="utf-8")
    thermlines = therm.readlines()
    therm.close()
    #print(thermlines)

    with open(therm_output, 'w') as thermfile:
        thermfile.write('! Thermodynamic data for C3Mech, generated on {0:%d/%m/%Y %H:%M:%S}.\n'.format(datetime))
        thermfile.write('THERMO ALL\n')
        thermfile.write('300.   1000.   5000.\n')
        for sp_name in species_list:
            #print(sp_name)
            counter = 0
            for line in range(len(thermlines)):
                #print(thermlines[line].split()[0])
                #print(line)
                if not thermlines[line].strip(' ') == '\n' and thermlines[line].startswith('!') == False:
                    if sp_name == thermlines[line].split()[0]:
                        #print(sp_name)
                        #print(thermlines[line])
                        thermfile.write(thermlines[line])
                        thermfile.write(thermlines[line+1])
                        thermfile.write(thermlines[line+2])
                        thermfile.write(thermlines[line+3])
                        counter+=1
                        break
                    else:
                        continue

            if counter == 0:
                #print(sp_name)
                print('{0} : thermochemistry data not found...'.format(sp_name))

        thermfile.write('END\n')

def clean_tran(tran_file, tran_output, species_list, datetime):
    tran = open(tran_file, 'r', encoding="utf-8")
    tranlines = tran.readlines()
    tran.close()

    with open(tran_output, 'w') as tranfile:
        tranfile.write('! Transport data for C3Mech, generated on {0:%d/%m/%Y %H:%M:%S}.\n'.format(datetime))
        for sp_name in species_list:
            counter = 0
            for line in range(len(tranlines)):
                if (tranlines[line].strip(' ') != '\n') and (tranlines[line].startswith('!') == False):
                    if sp_name == tranlines[line].split()[0]:
                        tran_component = tranlines[line].split()
                        if ('!' in tranlines[line]):
                            tranfile.write('{0:18}    {1}    {2:8.2f}    {3:6.2f}    {4:6.2f}    {5:6.2f}    {6:6.2f}    {7}\n'.format(tran_component[0],tran_component[1],float(tran_component[2]),float(tran_component[3]),float(tran_component[4]),float(tran_component[5]),float(tran_component[6]), ' '.join(tran_component[7:])))
                        else:
                            tranfile.write('{0:18}    {1}    {2:8.2f}    {3:6.2f}    {4:6.2f}    {5:6.2f}    {6:6.2f}\n'.format(tran_component[0],tran_component[1],float(tran_component[2]),float(tran_component[3]),float(tran_component[4]),float(tran_component[5]),float(tran_component[6])))
                        #tranfile.write('{0}\n'.format(tranlines[line].rstrip()))
                        counter+=1
                        break
                    else:
                        continue
            if counter == 0:
                #print(sp_name)
                print('{0} : transport data not found...'.format(sp_name))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    gen_name = "C3Mech"
    parser.add_argument("-m", "--minimal", help="generate mininal " + gen_name + ".THERM and " + gen_name + ".TRANS without duplicates",
                    action="store_false")
    target_directory = "AAA_Input_Kinetics"
    parser.add_argument("-c", "--copy", help="copy submodels to " + target_directory,
                    action="store_false")
    args = parser.parse_args()
    
    if(args.copy):
        copy_nuig_ht(target_directory)
        copy_nuig(target_directory)
        copy_pah(target_directory)
        copy_llnl(target_directory)
    files_list = get_files_list()
    species_list = make_species_list('AAA_Input_Kinetics', files_list)
    write_species_list(species_list)

    if args.minimal:
        DATETIME  = datetime.datetime.now()
        clean_therm('../SUBMECHANISMS/SOURCE-C3Mech.THERM','C3Mech.THERM',species_list, DATETIME)
        clean_tran('../SUBMECHANISMS/SOURCE-C3Mech.TRAN','C3Mech.TRAN',species_list, DATETIME)
