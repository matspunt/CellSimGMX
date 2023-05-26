# CELL_MODEL

This is a repo for the CELL_MODEL project. 

TODO: 
- Add tolerance to be read from .JSON into PackmolReader class.
- Make the JSON parser a general Class, and not a part of PACKMOL.
- In Gromacs_IO class, add "number_of_cells" from JSON into there.
=======
- Add tolerance to be read from .JSON into PackmolReader class. 
- Make the JSON parser a general Class, and not a part of PACKMOL. 
- In Gromacs_IO class, add "number_of_cells" from JSON into there. 
- Gromacs_Run class, generalize simulation setup independent of folder structure (add .mdp parser!)
- LAMMPS parser, figure out force field definitions in Martini --> use Moltemplate to prepare random Martini system. Or is Ksenia going to take care of this?

Multiple cells:
- Packing in a lattice or random orientation using either PACKMOL or gmx genconf. 

MSD analysis:
- Automated sanity checking --> talk about with Shreyas. 
