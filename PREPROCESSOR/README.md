# steps for mechanism compilation

Reminder: the mechanism compilation requires a couple of additional steps compared to standard CHEMKIN-like mechanism generally provided in the literature.
This is justified by the fact that the current full model contains a huge number of species and reactions that are unnecessary when e.g., high-temperature flames or small fuels such as hydrogen and ammonia are simulated.
Such compilation assembles only the necessary module blocks to facilitate and speed up kinetic simulations, while keeping the model as one. 

## 1. make sure to be able to execute `TailorSMOKEpp`

Modules are combined together with TailorSMOKE.

* WINDOWS: make sure that the `TailorSMOKEpp.exe` and `Run_Tailor_assemble.bat` are executable.  
* MAC, LINUX: TailorSMOKE bin and libraries are found in PREPROCESSOR/tailorsmokepp (bin, lib). Make sure these paths are added to environment variables (e.g., PATH, LD_LIBRARY_PATH).
(note: only works with gcc version >= 9)  

## 2. set up the module list in input_tot.dic

According to the desired mechanism type (see Appendix below), a different set of modules needs to be indicated.  
File to edit is `input_tot.dic`. Note: comments are set as //.

`@InputFolder`: input folder with modules (do not edit- it will be filled automatically, see step 3)  
`@OutputFolder`: output directory where assembled modules are stored.  
`@MechanismName`: name of the output mechanism in the output directory.  
`@CoreMechanism`: core mechanism: edit only for cantera users to: `NUIG_C0_LT-HT_Cantera.MECH`.  
`@SubMechanisms`: list of modules to be included. For suggested module lists: see Appendix below. Warning: `AAA_allspecies.txt` should not be commented.  

Note: for cantera users, 

## 3. Fill in the `AAA_Input_Kinetics` running the python script `write-modulesandspecies.py` 

to run the preprocessor:  
`python write-modulesandspecies.py` will

- fill the `AAA_Input_Kinetics` folder according to the modules listed in `input_tot.dic`.  
- write clean thermo and transport files as `C3Mech.THERM` and `C3Mech.TRAN`.  

## 4. Assemble the blocks with TailorSMOKE

- WINDOWS: run the `Run_Tailor_assemble.dat` batch script.  
- MAC, LINUX: run `TailorSMOKEpp.sh`  --input input_tot.txt
(for tailorsmoke installation: add bin and lib of tailorsmokepp.tar.gz to your paths; 
the Input_Abstractions folder will then contain the merged mechanisms in chemkin format (comments are included)  

Output: folder as indicated in `@OutputFolder` specification containing the combined blocks.  

Then, clean up the output (comments will be kept, but unnecessary TailorSMOKE output lines are removed):
`python clean-up-merge.py`  
Note: if name of the output folder is changed, it should be changed also in this script.
The output mech will then be rewritten.

The final mechanism can then be compiled with the preferred simulation code (CHEMKIN, CANTERA, OpenSMOKE++, FlameMaster, etc.)  
Example: with `OpenSMOKE++` on LINUX/MAC, run

`OpenSMOKEpp_CHEMKIN_PreProcessor.sh --input input_final.dic`.



# Appendix: Recommended blocks for mechanism compilation for different subsets

To facilitate compilation, here below are listed recommended subsets for different types of fuels. These subsets correspond to those used in the publication.
These can be directly pasted (or uncommented) in the  `@SubMechanisms` list in `input_tot.dic`.
Disclaimer: the authors are not responsible for incorrect use of such subsets for simulation cases incompatible with the intended purpose.


Note: all the following subsets require the `NUIG_C0_LT-HT.MECH` block set as `@CoreMechanism` in `input_tot.dic`, hence this module is omitted in the lists.


Reminder: `AAA_allspecies.txt` should never be removed from the `@SubMechanisms` list.

LT and HT below refer to low-temperature and high-temperature mechanisms, respectively.

## HT MAH/PAH formation (Premixed flames and counterflow laminar diffusion flames + LBVs + HT flow reactors and ST)

NUIG_C1-C2_LT-HT.MECH  
NUIG_C3-C4_HT.MECH   
NUIG_C5cy_HT.MECH  
NUIG_C5_HT.MECH  
NUIG_C6cy_HT.MECH  
PAH_BLOCK.CKI 
    
    
## LT MAH/PAH oxidation (LT JSR/Flow reactor oxidation + RCM)

NUIG_C1-C2_LT-HT.MECH  
NUIG_C3-C4_LT-HT.MECH  
NUIG_C5cy_LT-HT.MECH  
NUIG_C5_LT-HT.MECH  
NUIG_C6cy_LT-HT.MECH  
PAH_BLOCK.CKI  

## LT TPRF mixtures (Flow reactors and RCM)

NUIG_C1-C2_LT-HT.MECH  
NUIG_C3-C4_LT-HT.MECH  
NUIG_C5cy_LT-HT.MECH  
NUIG_C5_LT-HT.MECH  
NUIG_C6cy_LT-HT.MECH  
NUIG_C6_LT-HT.MECH  
NUIG_C7_LT-HT.MECH  
PAH_BLOCK.CKI  
LLNL_BLOCK.CKI  
   
  



