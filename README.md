# CELL_MODEL

This is a repo for the CELL_MODEL project. 

TODO: 

Proper Python package code:
- Gromacs_Run class, generalize simulation setup independent of folder structure (add .mdp parser!)
- LAMMPS parser, figure out force field definitions in Martini --> use Moltemplate to prepare random Martini system. Or is Ksenia going to take care of this?

Cell layer testing:
- Packing in a lattice or random orientation using either PACKMOL or gmx genconf. 
- Automated sanity checking, MSD analysis, energy tracking etc. --> talk about with Shreyas. 
- Add existing scripts to "cell_layer.py" as modules (Classes or functions) and with argparse so setup is easier. 

