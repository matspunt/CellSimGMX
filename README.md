# CELL_MODEL

## Explanation of repository

This is a repo for the CELL_MODEL project. The file structure of the repo is the following:

```
├── coord
├── log
├── mdps
├── params
├── README.md
├── src
│   ├── forcefield.itp
│   ├── gromacs.py
│   ├── input.json
│   ├── lammps.py
│   ├── main.py
│   ├── packmol.py
│   └── settings_parser.py
├── tools
│   ├── CELL.gro
│   ├── cell_layer_tools.py
│   ├── tools_README
│   └── using_tools.py
└── traj [UNSYNCED]
├── CELL_sim_database.csv
```

 ```src ``` contains the not-so modular code that hopefully will be packaged into a proper Python package at a later stage.  Within 'src', ```main.py ``` contains the current functionality of the CELL code.
 
  ```tools ``` contains modular code that is designed to help with setting up different GROMACS simulations. A description of the functions will follow later. 
  
  Then, all other directories (i.e.  ```coord ```,  ```log ```,  ```mdps ``` and  ```params ```) contain matching subdirectories whose name is defined by the main parameters in their simulations. A central database is utilized (```CELL_sim_database.csv```) in which the simulation is described, and some information about the simulations is written. A description of how this is formatted will follow later. 

## Future goals

- [ ] Add quality control metrics that can give insight in the behaviour of the simulations. These include things such as the I) averaged density of particles per a given area, II) the Radial Distribution Functions of a given particle type, III) the residence time of two neighbouring particles, IV) the geometric shape of a given combination of parameters, and V) the number of crashes reported for duplicates of that simulation. 
- [ ] Include a basic toolkit that can create multiple CELL packing (both crystal-like, and localized densities). Shreyas working on this. 
- [ ] Add rerun capabilites to the ```tools``` that is capable of decomposing the LJ interactions in interactions per bead type. 
- [ ] Regarding LAMMPS functionality, figure out the force field definition/parsing by looking at e.g. Moltemplate. 
- [ ] Add Shreyas his Python toy model to the ```tools``` code. 

## Cloning the repo and contributing

0. To be a collaborator, you need to be added to the repo by me (=Mats). If you haven't been added yet, let me know and I will do so. 
1. We use [PACKMOL](https://github.com/m3g/packmol) to build single cells in predefined compositions. One can install it locally or use a Python packaged version from ```conda``` or ```pip```. These may not be completely up to date with the PACKMOL codebase but we don't use any new features so this should be fine. Currently, there are no other dependencies. I suggest making a special conda environment for this project using ```conda create --name CELL``` and work within that. 
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

## Using the code or tools. 
1. The ```src``` code is self-explanatory, and classes are well-documented. Just check ```main.py```
2. The tools are separated in different files 
MORE WILL BE ADDED LATER WHEN THE TOOLS ARE REWRITTEN. 

