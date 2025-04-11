# C3MechV4.0
C3MechV4.0 release associated with the publication: XXXXXX

# Highlights
•	Unified combustion kinetic model for a wide range of fuel surrogate mixtures

•	Updated version of the widely adopted C3MechV3.3 combustion kinetic model

•	Simulates conventional and renewable fuels (alkanes, H2, NH3, DME, DMC/EC)

•	Extensive model testing for both pure fuels and fuel mixture experiments



# Repository structure


### SPECIES_DICTIONARY:
Species dictionary of C3Mech can be generated in this directory.  

### SUBMECHANISMS:
This directory includes all the latest (work-in-progress) submechanisms of each group. 
SOURCE-C3Mech.THERM/.TRAN are the only THERM/TRAN files to edit. write-species.py in PREPROCESSOR will read these files and generate "cleaned" THERM/TRAN files based on species included in given compiled mechanism.

### PREPROCESSOR: 
Mechanism compile can be done in this directory. Submechanisms and thermochemistry/transport propaties will be extracted from SUBMECHANISMS directory  
