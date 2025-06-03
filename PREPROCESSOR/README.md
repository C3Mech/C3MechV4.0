# Mechanism compilation

This directory provides the workflow to compile sub-models into a single file in CHEMKIN format.

## How to Use

1. **Select sub-models**  
 Open and edit `submodels.yaml` to specify the sub-models you want to include in the final mechanism.

2. **Compile the mechanism**  
 Run the following command in the terminal:
 ```sh
 ./compile_c3mech.py
   ```
The command writes the files `C3Mech.CKI,` `C3Mech.THERM,` and `C3Mech.TRAN` to the `output/` directory.

## Notes and Recommendations

- **Selecting necessary blocks:** Selecting only the necessary module blocks facilitates and speeds up kinetic simulations. 
- **User responsibility:** The user must select an appropriate model subset for a simulation case. 
- **Preselected combinations:** Frequently used sub-model combinations are provided in the `preselection/` directory. These subsets correspond to those used in associated publications where applicable.

For further questions or if you need additional help, please feel free to open an issue or reach out.