# CellSimGMX

## Explanation of repository

This is a repo for the CellSimGMX project. This branch does not use PACKMOL anymore to construct the cells but builds arbitrary cells without any external tools. 
Only GROMACS is required to run simulations. 

## Cloning the repo and contributing

1. 
2. First clone the repo in a local directory 
```sh
git clone --branch master https://github.com/matspunt/CELL_MODEL.git
```
 The ```master``` branch will be kept stable. Development will happen in the ```development``` branch. 
 
3. If you are contributing to the development branch, I recommend you to clone it in a separate directory since dynamic branches can get tricky with version control. To clone the ```development``` branch, simply repeat the command:
```sh
git clone --branch development https://github.com/matspunt/CELL_MODEL.git
```
4. Within the ```development``` branch, you can push commits to the branch, see here for instructions: https://linuxhint.com/push-to-specific-branch-in-git/
Or, alternatively, you can use an IDE like VS Code which has integrated git functionality

## Installing and using the code
1. To install, ```conda``` or ```pip``` virtual environment is recommended. Install from source ```pip install -e .``` in source directory (you need ```setuptools``` installed). 
2. Add pip console scripts to PATH if not already there: ```export PATH="~/.local/bin:$PATH"``` in ```~/.bashrc```
2. Package then runs with ```cellsimgmx```. Check options with ```cellsimgmx -h```, or you can use ```python3 -m cellsimgx``` if you do not want to use the console script.  
3. 

## Logic of the code and basic usage tutorial

The code has three objects "Cell", "Tissue" and "System", where "System" is the Tissue object with GMX functionality added. 

## GitHub todo
  
- Ensure ```settings_parser.py``` works with the new .JSON categories and settings (without PACKMOL!)
- Remove ```packmol.py```. Rewrite to ```cell.py```. This should have the functionality of packing arbitrary number of beads on certain shape surfaces (square, cube, ellipsoid, sphere) etc.
    Properties of cell shape should be configurable by ```input.json```. In addition, it should reference the ```forcefield.itp``` and distribute particle types accordingly. 
- Rewrite ```tissue.py``` to directly use a Cell object (instead of reading in the coordinate file from disk), and use the Cell object to create Tissues
- Include documentation using sphinx (or a similar framework)
- 

## Future goals

- [ ] Add quality control metrics that can give insight in the behaviour of the simulations. These include things such as the 
    * I) averaged density of particles per a given area
    * II) the Radial Distribution Functions of a given particle type
    * III) the residence time of two neighbouring particles
    * IV) the geometric shape of a given combination of parameters (obtained from clustering)
    * V) the number of crashes reported for duplicates of that simulation. 
- [ ] Add rerun capabilites to the ```tools``` that is capable of decomposing the LJ interactions in interactions per bead type. 
- [ ] Regarding LAMMPS functionality, figure out the force field definition/parsing by looking at e.g. Moltemplate. 
- Extend simulation logic to multiple CELLs (tissues)

## Basic parameter description

Currently all simulations in the GitHub were run at 310 K (NVT) and 1 bar (NPT). 

