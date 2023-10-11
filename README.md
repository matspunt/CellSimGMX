# CELL_MODEL

## Explanation of repository

This is a repo for the CELL_MODEL project. The file structure of the repo is the following:

```
├── coord
│   └── M1_mass72.0_r1.8_fk250.0_sigma0.47_eps1.0_md_Berendsen_no
│       ├── 1-min.gro
│       ├── 2-run.gro
│       ├── CELL.gro
│       └── CELL.xyz
├── log
│   └── M1_mass72.0_r1.8_fk250.0_sigma0.47_eps1.0_md_Berendsen_no
│       ├── 1-min.log
│       ├── 2-run.log
│       ├── gmx_run.log
│       ├── PACKMOL_build-06-06-10-41-24.log
├── params
│   └── M1_mass72.0_r1.8_fk250.0_sigma0.47_eps1.0_md_Berendsen_no
│       ├── CELL.top
│       └── forcefield.itp
├── README.md
├── settings
│   └── M1_mass72.0_r1.8_fk250.0_sigma0.47_eps1.0_md_Berendsen_no
│       ├── 1-min.tpr
│       ├── 2-run.tpr
│       ├── min.mdp
│       └── run.mdp
├── simulations_CELL.csv
├── src
│   ├── gromacs.py
│   ├── input.json
│   ├── lammps.py
│   ├── main.py
│   ├── packmol.py
│   └── settings_parser.py
├── tools
│   ├── CELL_analysis_library.py
│   ├── CELL_tools_library.py
│   ├── dict_examples.py
│   ├── template_param_range.py
│   └── template_using_analysis.py
└── traj
    ├── UNSYNCED, BUT CONTAINS TRAJECTORIES


```

 ```src ``` contains code that hopefully will be packaged into a proper Python package at a later stage.
 
  ```tools ``` contains code that is designed to help with setting up different GROMACS simulations.
  
  Then, all other directories (i.e.  ```log ```,  ```mdps ``` and  ```params ```) contain matching subdirectories whose name is defined by the main parameters in their simulations. The coordinates (```log ```) are hosted on ```/BIOBACKUP/biofys4/common/Mats/CELL_MODEL```. A central database is utilized (```simulations_CELL.csv```) in which the simulations are described, and some information about the simulations is written. A description of how this is formatted will follow later. 

## NOTE
  
- THIS BRANCH IS OBSOLETE AND WILL BE SLOWLY INTEGRATED INTO ```no_packmol``` branch!