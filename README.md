# CellSimGMX

<p align="left">
  <img src="docs/logo.png" width="550" title="CellSimGMX logo">
</p>

## CellSimGMX

This is a repo for the CellSimGMX project. Only basic Python modules (e.g. numpy, matplotlib) are necessary to use the toolkit. Gromacs is solely used to run simulations, and not for input file generation. This would allow compatability with other MD engines through coordinate and topology conversion.

**THIS IS THE DEVELOPMENT BRANCH! Changes will be merged into the master branch if stable/functional**
**Do not blindly merge development with master!**

## Documentation/installation

Note: documentation is currently built from the master branch. 

Documentation: https://cellsimgmx.readthedocs.io

Use gromacs-2021.6 --> VS instabilities. 

## Notes/todos
- Note: stress test tissue packing logic for large number of cells!!
- Apply weak position restraint to matrix particles at the .MDP level. 
- Although relative paths work fine, absolute paths are probably best in feeding to the programme. 

## Project to do's applied to the code

Essential (High priority):

- [ ] Implement virtual site logic --> might need gromacs-2021.6 for this. 
- [ ] Generate randomly rotated cells within a Layer 
- [ ] Write  ```gromacs.py``` which takes all input files and actually runs and sets up the simulations with error handling --> add to logger!!!
- [ ] Allow the selection of different tissue layers (like in epidermis). Cell type in a layer would be identical. At the cell level, multiple cells have to be generated (with also allowed structural differences in nr of particles, radius and shape). Layer could be a new object between Cell and Tissue?
- [ ] Fix disordered packing bug!
- [ ] Expand on the matrix options. Distance between particles (eliminate vacuum), and expose the matrix scaling factor to the .JSON. 

Medium priority:
- [ ] Include: https://github.com/shirtsgroup/physical_validation to assess physical reasonableness systems. Perhaps this is overkill, and it would make more sense to look at validation of force field parameters. 
- [ ] Introduce analysis modules which either wrap GMX tools or our own implementation which calculate some properties of the systems. 
- [ ] Introduce a nucleus within the Cell (a mini cell with smaller sigma?) at random positions in the tisuse. Keep the center particle independent of the nucleus. 
- [ ] Allow coordinate files to be written as PDBs with CONECT records (to visualize surface bonds). 

## Analysis modules

- [ ] Add quality control metrics that can give insight in the behaviour of the simulations. These include things such as the 
    * I) averaged density of particles per a given area
    * II) the Radial Distribution Functions of a given particle type
    * III) the residence time of two neighbouring particles
    * IV) the geometric shape of a given combination of parameters (obtained from clustering)
    * V) the number of crashes reported for duplicates of that simulation. 
- [ ] A metric which can somehow quantify the shape within bulk tissue (hard to inspect visually)\
- [ ] A metric which assesses holes occuring on the surface of cells. 

## Basic model description

- Running simulations at 310 K (NVT) and 1 bar (NPT), using the v-rescale thermostat and c-rescale barostat. Nose-hoover is not a good idea with leapfrog due to sometimes simulating few DOF. 

- A cell radius of 1.85 nm is currently used. In combination with this, a force constant of 250 needs to be used at minimum at <40 fs time steps to prevent energy leakage in NVE. 

- Sigma = 0.47 nm (Martini 3 R-bead) with epsilon varying from 2 - 10 currently. (5.5 is the most attractive in M3). No 'electrostatic' interactions. 

- 180 surface particles connected to the center through harmonic springs. Perhaps we will add surface springs in the model (though have to be careful about substantiating this from biology POV). 

- Tissue can be packed in various initial shapes, and various packing arrangements. 