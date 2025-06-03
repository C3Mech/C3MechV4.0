# C3MechV4.0

C3MechV4.0 release associated with the publication: **XXXXXX**.

# Highlights

- **Unified combustion kinetic model** covering a wide range of fuel surrogate mixtures  
- **Updated version** of the widely adopted C3MechV3.3 combustion kinetic model  
- **Supports both conventional and renewable fuels** (alkanes, H₂, NH₃, DME, DMC/EC)  
- **Extensive testing** for pure fuels and fuel mixture experiments  

# Repository structure

### SPECIES_DICTIONARY:
Contains the species dictionary for C3Mech.  
You can run the provided script to check constraints and generate a PDF overview of species.

### SUBMECHANISMS:
Holds the latest submechanisms for each group.  
- **SOURCE-C3Mech.THERM** and **SOURCE-C3Mech.TRAN** are the only THERM and TRAN files that require editing. The script `compile_c3mech.py` in the `PREPROCESSOR` directory reads these files and generates cleaned versions based on the species included in the final compiled mechanism.  
- Recommended modules for various fuel types are noted in the `PREPROCESSOR/submodels.yaml`.

### PREPROCESSOR:
Used to compile mechanisms. Submechanisms, thermochemistry, and transport properties are extracted from the `SUBMECHANISMS` directory. Please refer to the `PREPROCESSOR/README.md` for detailed instructions on how to compile mechanisms. Within `PREPROCESSOR`, the **`preselection/`** subdirectory provides a set of commonly used sub-model combinations. These files are offered for convenience so that they can be downloaded and used directly.

