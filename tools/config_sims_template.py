'''Template script for setting up simulations with the CELL tools code in GROMACS.

This assumes you have installed PACKMOL (the script will print the PACKMOL binary that is used)
You can clone the GitHub in a suitable directory in '/wrk'. Then, you can configure four aspects of the simulation:

1) Choose a jobname, all simulation information will be saved under this jobname. This hopefully makes it easier to locate
specific simulations when the .csv starts filling up.
2) Choose your forcefield parameters. You can keep this constant, or loop over the dictionary and configure different settings
for each simulation. 
3) Choose PACKMOL settings. If you want to use more membrane beads, don't forget to change the 'packmol_dict' to reflect this.
The final GMX topology is constructed based on the bead names in the "CELL.gro' that PACKMOL creates. 
4) Choose the mdp run settings of your simulations. By default, the standard Martini MDP settings are used. You can change
any setting that you like by changing it in a mdp dictionary. 

You can check "dict_examples.py" for basic examples on how to modify these configurations for multiple bead types. This current
script restricts itself to single cells, with only a single membrane bead type
'''

from CELL_tools import DatabaseMaintenance
from CELL_tools import ForceFieldParser
from CELL_tools import InputSingleCell
from CELL_tools import PackmolExecuterSingleCell
from CELL_tools import Gromacs_IO
from CELL_tools import Gromacs_MDP
from CELL_tools import GromacsRun

##CONFIGURATION 1: Choose a descriptive job name
jobname = "Figuring out force constant between 2000-5000" 

##CONFIGURATION 2: choosing the force field parameters
ff_single_bead_type = {
    'bead_types': {
        'M1': {'name': 'M1', 'mass': 72},
        'C': {'name': 'C', 'mass': 72}
    },
    'bond_types': {
        'C_M1': {'r0': 1.85, 'fk': 1000}
    },
    'nonbond_params': {
        'M1_M1': {'sigma': 0.47, 'epsilon': 2}
    }
}

##CONFIGURATION 3: configuring PACKMOL settings
packmol_dict = {
        'settings': {
            'tolerance' : '4', #tolerance in angstrom
            'radius' : '18', #radius in angstrom
            'shape' : 'sphere' #cube, ellipsoid are also supported, but untested
    },
        'beads': {
        'M1': {'name': 'M1', 'number': 180},
        'C': {'name': 'C', 'number': 1},
    }
}

##CONFIGURATION 4: let's use the Langevin integrator for all simulations
mdp_langevin = {
    'integrator': 'sd',
    'tcoupl': 'no',
    'tau_t': '1.0',
    'ref_t': '298',
    'ld-seed' : '-1',            
}

PackmolExecuterSingleCell.check_packmol_path()

parser = ForceFieldParser(ff_single_bead_type)
parser.generate_forcefield_itp()

packmol = PackmolExecuterSingleCell(packmol_dict)
packmol.run_packmol_single_CELL()

GMX = Gromacs_IO()
GMX.convert_xyz_to_gro("CELL.xyz", "CELL.gro", box_size = 5)
GMX.build_GMX_top_single_CELL('CELL.gro', 'forcefield.itp')

Gromacs_MDP.write_min_mdp_file('min.mdp', 'cg', '5000', '5', '0.01') 
Gromacs_MDP.write_run_mdp_file('run.mdp', '5000', mdp_langevin) #50000000 for 1 us

GromacsRun.run_GMX_basic(filename='CELL.gro', timeout=60, nr_of_threads=12) #2400 is default for 1 us

DatabaseMaintenance.move_files_to_parentdir('.', '../log', ['.log'])
DatabaseMaintenance.move_files_to_parentdir('.', '../params', ['.top', '.itp'])
DatabaseMaintenance.move_files_to_parentdir('.', '../traj', ['.trr', '.edr', '.xtc'])
DatabaseMaintenance.move_files_to_parentdir('.', '../settings', ['.mdp', '.tpr'])
DatabaseMaintenance.move_files_to_parentdir('.', '../coord', ['.gro', '.xyz'])
DatabaseMaintenance.del_files_with_extension('.', '.cpt') #we do not need the checkpoints, feel free to uncomment though