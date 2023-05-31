# CELL_MODEL

This is a repo for the CELL_MODEL project. The repo is organized in the following way: 

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
└── traj
```

NOTE:

- "src" contains actual python package files and functionality
- "tools" contains useful tools in setting up simulations and testing the project model

TODO: 

Proper Python package code:
- Gromacs_Run class, generalize simulation setup independent of folder structure
- Add multiple CELL parsing in the code. 
- LAMMPS parser, figure out force field definitions in Martini --> use Moltemplate to prepare random Martini system. Or is Ksenia going to take care of this?

Project tools:

- Create PACKMOL wrapper for building ordered systems of multiple cells, and random systems of multiple cells
- Add on to tools in "cell_layer.py"
- Expand on GitHub by saving project progress here, also so other people can access them. 

## Getting Started


### Requirements


### 

1. 
2. Clone the repo
   ```sh
   git clone https://github.com/your_username_/Project-Name.git
   ```
3. Install NPM packages
   ```sh
   npm install
   ```
4. Enter your API in `config.js`
   ```js
   const API_KEY = 'ENTER YOUR API';
   ```



## Future development goals

- [] Add back to top links
- [ ] Add Additional Templates w/ Examples
- [ ] Add "components" document to easily copy & paste sections of the readme

