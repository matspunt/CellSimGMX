# CellSimGMX

<p align="left">
  <img src="docs/logo.png" width="550" title="CellSimGMX logo">
</p>

## CellSimGMX

**Deprecated repo. For latest code and documentation, see: https://github.com/Vattulainen-Group/gmx_CellModel/wiki**

This is a repo for the CellSimGMX project. Only basic Python modules (e.g. numpy, matplotlib) are necessary to use the toolkit. Gromacs is solely used to run simulations, and not for input file generation. This would allow compatability with other MD engines through coordinate and topology conversion.
## Documentation/installation

Note: documentation is currently bcuilt from the master branch. 

Documentation: https://cellsimgmx.readthedocs.io

Use gromacs-2021.6 --> when running with VS / Go like potentials. 

## Notes/todos
- Note: stress test tissue packing logic for large number of cells!!
- Apply weak position restraint to matrix particles at the .MDP level. 

## Project to do's applied to the code

Essential (High priority):

- [ ] For automated simulation running, find a way to deal with freezing simulations (and report in cellsimgmx.log file?) --> look at updates of file on disk
- [ ] Write a wrapper for running multiple simulations. 
- [ ] Repopulate repo. 
- [ ] Implement virtual site logic. 
- [ ] Delete master branch. 
- [ ] Update documentation
- [ ] Tabulated potential as an option for itp parsing. 
- [ ] Center bead with large sigma (sigma is 10% volume of whole cell) (as a particle type), in first iteration at geometrical center [NUCLEUS option in JSON file]
- [ ] matrix --> grid but also hexagonally packed (+spacing between particles + number of layers + spacing between layers + distance to first cell layer on matrix)
- [ ] Allow volume to be filled with cytoplasmic beads. 
- [ ] Generate randomly rotated cells within a Layer 
- [ ] Allow the selection of different tissue layers (like in epidermis). Cell type in a layer would be identical. At the cell level, multiple cells have to be generated (with also allowed structural differences in nr of particles, radius and shape). Layer could be a new object between Cell and Tissue?
- [ ] replace particles randomly with different radii!!!
- [ ] Expand on the matrix options. Distance between particles (eliminate vacuum), and expose the matrix scaling factor to the .JSON. 

Medium priority:
- [ ] Include: https://github.com/shirtsgroup/physical_validation to assess physical reasonableness systems. Perhaps this is overkill, and it would make more sense to look at validation of force field parameters. 
- [ ] Allow coordinate files to be written as PDBs with CONECT records (to visualize surface bonds). 
- [ ] Add Saana's code to GitHub as analysis module. 

## Basic model description

- Running simulations at 310 K (NVT) and 1 bar (NPT), using the v-rescale thermostat and c-rescale barostat. Nose-hoover is not a good idea with leapfrog due to sometimes simulating few DOF. 

- A cell radius of 1.85 nm is currently used. In combination with this, a force constant of 250 needs to be used at minimum at <40 fs time steps to prevent energy leakage in NVE. 

- Sigma = 0.47 nm (Martini 3 R-bead) with epsilon varying from 2 - 10 currently. (5.5 is the most attractive in M3). No 'electrostatic' interactions. 

- 180 surface particles connected to the center through harmonic springs. Perhaps we will add surface springs in the model (though have to be careful about substantiating this from biology POV). 

- Tissue can be packed in various initial shapes, and various packing arrangements. 
