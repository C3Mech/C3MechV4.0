# Mechanism compilation

This directory provides the workflow to compile sub-models into a single file in CHEMKIN format.

## How to use

1. **Select sub-models**  
 Open and edit `submodels.yaml` to specify the sub-models you want to include in the final mechanism.

2. **Compile the mechanism**  
 Run the following command in the terminal:
 ```sh
 ./compile_c3mech.py
   ```
The command writes the files `C3Mech.CKI,` `C3Mech.THERM,` and `C3Mech.TRAN` to the `output/` directory.

## Installing dependencies

You can install the dependencies, for example, with:

```sh
pip3 install pyyaml --user
```

## Notes and recommendations

- **Selecting necessary blocks:** Selecting only the necessary module blocks facilitates and speeds up kinetic simulations. 
- **User responsibility:** The user must select an appropriate model subset for a simulation case. 

If you have further questions or need additional help, please feel free to open an issue or [contact us](r.langer@itv.rwth-aachen.de).
