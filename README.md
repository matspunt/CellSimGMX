# CellSimGMX

<p align="left">
  <img src="docs/logo.png" width="550" title="CellSimGMX logo">
</p>

## Explanation of repository

This is a repo for the CellSimGMX project. This branch does not use PACKMOL anymore to construct the cells but builds arbitrary cells & tissues without any external tools. 
Only basic Python modules (e.g. numpy, matplotlib) are necessary to use the toolkit. Gromacs is solely used
to run simulations, and not for input file generation. This would allow compatability with other MD engines through
coordinate and topology conversion. 

## Accessing and installing the code

1. First clone the repo in a local directory 
```sh
git clone https://github.com/matspunt/CellSimGMX.git
```
2. By default you should be on the ```no_packmol``` branch, but just to make sure run ```git checkout no_packmol```
3. To install, ```conda``` or ```pip``` virtual environment is recommended: ```conda create -n cellsimgmx```
* If NOT using ```conda``` environments, run ```pip install .``` in the source directory (dir with ```setup.py```). 
* If using ```conda``` env, you need to make sure the system ```pip``` is NOT used by running ```conda info --envs``` and explicitly linking to the pip in that environment:
```/path/to/miniconda3/envs/cellsimgmx/bin/pip pip install .``` 
4. Add pip console scripts to PATH if not already there: ```export PATH="~/.local/bin:$PATH"``` in ```~/.bashrc```
5. You can now run the programme with ```cellsimgmx``` or ```cellsimgmx -h```, or you can use ```python3 -m cellsimgx``` if you do not want to use the console script.  

## Inputs of the package

The code comes with a force field (```forcefield.itp```) for all the VdW particles in the system. An input file ```input.json``` is used to determine the specifications of the system
used for the simulation (amount of junction particles, distribution of them etc). The parameter descriptions should speak for themselves. 

## Bugs / to dos
- Disordered tissue packing sometimes breaks cells due to wrong distance treshold logic (easy to fix)
- Matrix input for now only allows a single atom name ('MX1') and a single layer. 
- 

## Notes for Mats

- Note: stress test tissue packing logic for large number of cells!!
- Include documentation. 
- Apply weak position restraint to matrix particles. 
- Explain equilibration routine in code documentation!!!!
- Explain naming conventions documentation!!!!
- Optimize matrix distance based on simulation logic.
- Although relative paths work fine, absolute paths are probably best in feeding to the programme. 

## Project to do's applied to the code

Essential features (**HIGH PRIORITY**)
- [x] Ensure ```settings_parser.py``` works with the new .JSON categories and settings. Tthe JSON parsing simply reads all variables and saves them as varnames. This way it can be extended easily without having to rewrite the parsing. 
- [x] Remove PACKMOL functionality ```cell.py```. This should have the functionality of packing arbitrary number of beads on certain shape surfaces (square, cube, ellipsoid, sphere) etc.
- [x] In addition, it distributes particle types in a systematic way and builds the GROMACS topology from the forcefield.
- [x] Within this logic, also include generation of the surface bonds (can be configured for arbitrary number of n-neighbours). 
- [x] Create various Tissue packings (```tissue.py```) based on ```input.json``` settings.
- [x] Change how the matrix is defined (as individual particles) and create logic that works for matrix setting turned on or off. 
- [ ] Write ```simulation.py``` or ```gromacs.py``` which takes all input files and actually runs and sets up the simulations with error handling. 
- [x] Come up with files and simulation directory layout. How is everything run and named? How should the files be organized for the user?

In the future (**MEDIUM PRIORITY**)
- [ ] Allow inhomogeneous system creation from different cell types (i.e. with a different distribution of junction beads, VdW properties etc) instead of merely duplicating a single type. 
- [ ] Introduce analysis modules which either wrap GMX tools or our own implementation which calculate some properties of the systems. 
- [ ] Introduce a nucleus within the Cell (a mini cell with smaller sigma?) at random positions in the tisuse. Keep the center particle independent of the nucleus. 

Nice-to-have's (**LOW PRIORITY**)
- [ ] Include documentation of the code using sphinx / doxygen. 
- [ ] Make installation compatible with conda (currently only pip is supported)

## Analysis modules

- [ ] Add quality control metrics that can give insight in the behaviour of the simulations. These include things such as the 
    * I) averaged density of particles per a given area
    * II) the Radial Distribution Functions of a given particle type
    * III) the residence time of two neighbouring particles
    * IV) the geometric shape of a given combination of parameters (obtained from clustering)
    * V) the number of crashes reported for duplicates of that simulation. 
- [ ] A metric which can somehow quantify the shape within bulk tissue (hard to inspect visually)\
- [ ] A metric which assesses holes occuring on the surface of cells. 
- [ ] ????

## Basic model description

- Currently all simulations in the GitHub were run at 310 K (NVT) and 1 bar (NPT), using the v-rescale thermostat and Parrinello-Rahman barostat. 

- A cell radius of 1.85 nm is currently used. In combination with this, a force constant of 250 needs to be used at minimum at <40 fs time steps to prevent energy leakage in NVE. 

- Sigma = 0.47 nm (Martini 3 R-bead) with epsilon varying from 2 - 10 currently. (5.5 is the most attractive in M3). No 'electrostatic' interactions. 

- 180 surface particles connected to the center through harmonic springs. Perhaps we will add surface springs in the model (though have to be careful about substantiating this from biology POV). 

- Tissue can be packed in various initial shapes, and various packing arrangements. 