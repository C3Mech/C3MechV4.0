# C3MechV4.0

**C3MechV4.0** is a universal chemical kinetic model developed for both conventional and renewable fuels (e.g., alkanes, hydrogen, ammonia, dimethyl ether, dimethyl carbonate/ethylene carbonate). It has been validated for a wide range of conditions and uses a hierarchical structure that allows the targeted compilation of sub-models. This repository is linked to the publication ["Modeling Combustion Chemistry using C3MechV4.0: an extension to mixtures of hydrogen, ammonia, alkanes, and cycloalkanes"]() and provides the following data and resources:

- All sub-models maintained by the C3 consortium  
- A species directory (CSV & PDF), along with a Python script for checking species data and converting them into a PDF  
- A compiler script to generate CHEMKIN mechanisms from user-selected sub-models  
- Several precompiled sub-model sets (in CHEMKIN and Cantera YAML) for direct use  

## Repository contents

1. **SUBMECHANISMS**  
This directory contains sub-models ("mechanism" files with kinetic data) organized by the research groups that developed them and the jointly used files `SOURCE-C3Mech.THERM` and `SOURCE-C3Mech.TRAN` that provide all thermochemistry and transport data for C3MechV4.0. These files are recommended as a starting point for model development based on C3MechV4.0.

2. **SPECIES_DICTIONARY**  
This directory provides a CSV file with plain species data and a PDF generated from these, which serves as an overview. The Python script `make_species_dict.py` implements the C3Mech species model and checks constraints on the data in the CSV file. It can optionally generate LaTeX code and plots of species structures, which can then be compiled into a PDF. The script can also be adapted to manage species data in other chemical kinetic models. Refer to the `README.md` for more details on the species model, the data in the CSV file, and the script's usage. 

3. **COMPILER**  
This directory includes the script `compile_c3mech.py` and instructions on creating a single CHEMKIN-format mechanism from user-selected sub-models. For details on usage, refer to the `README.md`.

4. **PRESELECTION**  
This directory contains frequently used sub-model combinations that have already been compiled in CHEMKIN and Cantera YAML formats. They can be used directly without the additional compilation step and are provided for convenience.

For questions or suggestions, please open an issue or [contact us](mailto:r.langer@itv.rwth-aachen.de).
