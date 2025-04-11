# Species dictionary

All species in the C3MechV4.0 are defined in `species_dict.csv`. This CSV file has
seven columns with the following meanings:
- `model_name`: name used in the CHEMKIN files (upper case ASCII string)
- `inchi`:  International Chemical Identifier (ASCII string)
- `smiles`: a string following the simplified molecular-input line-entry system 
  (ASCII string)
- `excited`: `1` if excited; otherwise, `0` (boolean)
- `multiplicity`: multiplicity of a species' energy level (`0` if unspecified,
  `1` for singlets, `2` for doublets, `3` for triplets, etc.) (non-negative 
  integer)
- `lumped`: `1` if the `model_name` represents multiple InChIs that should be 
  tracked in the dictionary; otherwise, `0` (boolean)
- `rdkit_stereochemistry`: `1` for RDKit-compatible stereochemistry information; 
  otherwise `0` (boolean)

The combination `inchi`-`excited`-`multiplicity` must be unique for each row.
Simplified identifiers neglecting stereochemistry information derived from the
values in the `inchi` and `smiles` columns must still be RDKit compatible and 
consistent if the value in the `stereochemistry` column is `0`.

# Installing dependencies

You can install the dependencies, for example, with:

```sh
pip3 install rdkit CairoSVG gitpython pyyaml pandas numpy yapf --user
```

Alternatively, you can use conda:

```sh
conda env create -f environment.yaml
```

The scripts were tested with Python versions 3.8.5 and 3.9.6.

# Generating the species dictionary

The following commands generate the output file `species_dict.pdf`:

```sh
# Select the relevant submodels for the species dictionary the yaml input (default is submodels.yaml).
# Species images will only be updated if they are older than the species dictionary CSV file. 
# Erase the content of the output directory to force a regeneration of the species images. 
./make_species_dict.py
# use ./make_species_dict.py -d to perform checks without output generation (faster)
cd output
pdflatex species_dict.tex
```
