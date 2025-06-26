# BuildCstruct - Armchair Carbon Nanotube Generator

## Overview

BuildCstruct is a Python 3 script for generating armchair carbon nanotube (CNT) structures with optional functionalization. This is a modified version of the original buildCstruct 1.2 by Andrea Minoia and Martin Voegele, streamlined to focus specifically on armchair CNTs and made compatible with Python 3.

**Modified by:** mdanh  
**Date:** June 26th 2025  
**Original Authors:** Andrea Minoia, Martin Voegele  
**Original Date:** January 27th 2018

## Features

- **Armchair CNT Generation**: Creates single-wall armchair carbon nanotubes with customizable diameter and length
- **Functionalization Support**: Adds functional groups (OH, COOH, COO⁻) to CNT edges
- **Multiple Output Formats**: Supports XYZ, GRO, and MOL2 file formats
- **Charge Assignment**: Automatically assigns partial charges based on functionalization type
- **Connectivity Information**: Generates bond connectivity for MOL2 format files

## Dependencies

- **Python 3.x**
- **NumPy**: For mathematical operations and array handling
- **Standard Python libraries**: sys, optparse, os

## Installation

No installation required. Simply ensure you have Python 3 and NumPy installed:

```bash
pip install numpy
```

## Usage

### Basic Syntax

```bash
python buildCstruct_modified.py [options] output_file
```

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-c, --credits` | Display credits and exit | - |
| `-g, --geometry INDEX LENGTH` | Specify CNT geometry (chiral index and length in Å) | Required |
| `-f, --funct TYPE` | Functionalization type: none, oh, cooh, coo | none |
| `--xyz` | Output in XYZ format | - |
| `--gro` | Output in GROMACS GRO format | - |
| `--mol2` | Output in MOL2 format | - |

### Examples

#### 1. Basic Armchair CNT (6,6) with 40 Å length
```bash
python buildCstruct_modified.py -g 6 40 --xyz cnt_6_6.xyz
```

#### 2. Hydroxyl-functionalized CNT
```bash
python buildCstruct_modified.py -g 8 50 -f oh --mol2 cnt_8_8_oh.mol2
```

#### 3. Carboxyl-functionalized CNT (COOH)
```bash
python buildCstruct_modified.py -g 10 60 -f cooh --mol2 cnt_10_10_cooh.mol2
```

#### 4. Carboxylate-functionalized CNT (COO⁻)
```bash
python buildCstruct_modified.py -g 6 40 -f coo --mol2 cnt_6_6_coo.mol2
```

## Technical Details

### Armchair CNT Structure

Armchair carbon nanotubes are characterized by their chiral vector (n,n), where n is the chiral index. The script generates CNTs with:

- **Diameter**: Calculated based on the chiral index and C-C bond length (1.3874 Å)
- **Length**: User-specified in Angstroms
- **Carbon arrangement**: Armchair edge configuration

### Functionalization Chemistry

#### 1. None (Default)
- **Description**: Pristine CNT with hydrogen termination at edges
- **Edge atoms**: Saturated with hydrogen atoms
- **Charge**: Neutral carbons (-0.16 e), hydrogens (+0.13 e)

#### 2. Hydroxyl (-OH)
- **Description**: Hydroxyl groups attached to every second carbon at CNT edges
- **Chemistry**: C-OH bonds formed at edge carbons
- **Charges**: 
  - Functionalized carbons: varied charges
  - Oxygen: -0.53 e
  - OH hydrogen: +0.37 e
  - Remaining edge carbons: hydrogen-terminated

#### 3. Carboxyl (-COOH)
- **Description**: Carboxylic acid groups at every second edge carbon
- **Chemistry**: C-COOH bonds with protonated carboxyl groups
- **Charges**:
  - Carboxyl carbon: +0.7 e
  - Oxygen (C=O): -0.55 e
  - Oxygen (C-OH): -0.6 e
  - OH hydrogen: +0.44 e

#### 4. Carboxylate (-COO⁻)
- **Description**: Deprotonated carboxylic acid groups (carboxylate anions)
- **Chemistry**: C-COO⁻ bonds without protons
- **Charges**:
  - Carboxyl carbon: +0.83 e
  - Both oxygens: -0.84 e each

### Coordinate System

- **X-axis**: Radial direction (CNT cross-section)
- **Y-axis**: CNT axis (length direction)
- **Z-axis**: Radial direction (CNT cross-section)
- **Origin**: CNT center

## Code Structure

### Main Functions

#### `main()`
- **Purpose**: Main execution function
- **Process**: 
  1. Parse command line arguments
  2. Generate CNT structure
  3. Apply functionalization
  4. Write output file

#### `armcnt(n, l, ccbond, funct)`
- **Purpose**: Generate armchair CNT coordinates
- **Parameters**:
  - `n`: Chiral index
  - `l`: CNT length (Å)
  - `ccbond`: C-C bond length (1.3874 Å)
  - `funct`: Functionalization type
- **Returns**: Coordinates, number of atoms per ring, periodic length, total atoms

#### Functionalization Functions

##### `add_H(coords, natx, funct_OH)`
- **Purpose**: Add hydrogen atoms or hydroxyl groups to CNT edges
- **Parameters**:
  - `coords`: Atomic coordinates list
  - `natx`: Number of atoms per ring
  - `funct_OH`: Boolean for OH functionalization

##### `add_COO(coords, natx, is_protonated)`
- **Purpose**: Add carboxyl/carboxylate groups to CNT edges
- **Parameters**:
  - `coords`: Atomic coordinates list
  - `natx`: Number of atoms per ring
  - `is_protonated`: Boolean for COOH vs COO⁻

#### Output Functions

##### `write_xyz(file, data)`
- **Purpose**: Write XYZ format file
- **Format**: Simple atomic coordinates

##### `write_gro(file, data, pbc1)`
- **Purpose**: Write GROMACS GRO format file
- **Format**: GROMACS structure with box information

##### `write_mol2(file, data, bondlist)`
- **Purpose**: Write MOL2 format file
- **Format**: Complete molecular structure with bonds and charges

#### Utility Functions

##### `connect(coords, natx, nohcoords)`
- **Purpose**: Calculate bond connectivity
- **Method**: Distance-based bond detection
- **Returns**: Connectivity matrix and bond list

##### `getdist(at1, at2)`
- **Purpose**: Calculate distance between two atoms
- **Returns**: Euclidean distance in Angstroms

## File Formats

### XYZ Format
```
<number_of_atoms>
Generated by buildCstruct - Armchair CNT
<atom_symbol> <x> <y> <z>
...
```

### GRO Format
```
Generated by buildCstruct - Armchair CNT
<number_of_atoms>
<residue_id><residue_name><atom_name><atom_id> <x> <y> <z>
...
<box_x> <box_y> <box_z>
```

### MOL2 Format
```
@<TRIPOS>MOLECULE
CNT
<atoms> <bonds> 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
<atom_id> <atom_name> <x> <y> <z> <atom_type> <subst_id> <subst_name> <charge>
...

@<TRIPOS>BOND
<bond_id> <atom1> <atom2> <bond_type>
...
```

## Charge Assignment Details

The script uses GAFF (General Amber Force Field) compatible charges:

### Carbon Atoms
- **Bulk carbons**: 0.0 e (neutral aromatic carbons)
- **Edge carbons**: Varies by functionalization and position
- **Carboxyl carbons**: +0.7 to +0.83 e depending on protonation

### Hydrogen Atoms
- **Alkyl hydrogens**: +0.13 to +0.22 e
- **Hydroxyl hydrogens**: +0.37 to +0.44 e

### Oxygen Atoms
- **Hydroxyl oxygen**: -0.53 e
- **Carbonyl oxygen**: -0.55 to -0.58 e
- **Carboxylate oxygen**: -0.84 to -0.86 e

## Geometric Parameters

- **C-C bond length**: 1.3874 Å (aromatic)
- **C-H bond length**: 1.087 Å
- **C-O bond length**: 1.362 Å
- **O-H bond length**: 0.974 Å
- **Bond angles**: 120° (sp² hybridization)

## Limitations

1. **Single-wall only**: Does not generate multi-wall CNTs
2. **Armchair only**: No zigzag or chiral CNTs
3. **Edge functionalization only**: No sidewall functionalization
4. **Non-periodic**: Generates finite-length CNTs only

## Troubleshooting

### Common Issues

1. **"Numpy not installed"**
   - Solution: Install NumPy with `pip install numpy`

2. **"You have given me more than one argument"**
   - Solution: Ensure output filename is the last argument

3. **Missing geometry parameters**
   - Solution: Always specify `-g index length`

### Error Messages

- **Missing arguments**: Use `--help` to see required parameters
- **Invalid functionalization**: Use one of: none, oh, cooh, coo
- **File write errors**: Check file permissions and disk space

## References

#  Literature:

 - M. Vögele, J. Köfinger, G. Hummer: 
   Simulations of Carbon Nanotube Porins in Lipid Bilayers.
   Faraday Discuss., 2018, Accepted Manuscript, DOI: 10.1039/C8FD00011E  
   http://pubs.rsc.org/en/content/articlelanding/2018/fd/c8fd00011e

(buildCstruct)
 - A. Minoia, L. Chen, D. Beljonne, L. Lazzaroni:
   Molecular Modeling Study of the Structure andStability of Polymer/Carbon Nanotube Interfaces.
   Polymer 53 (2012) 5480-5490
