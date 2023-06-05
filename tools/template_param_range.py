'''Template script for setting up simulations with the CELL tools code in GROMACS.

You can configure four aspects of the simulation:

1) Choose a jobname, all simulation information will be saved under this jobname. This hopefully makes it easier to locate
specific simulations when the .csv starts filling up.
2) Choose your force field parameters. You can keep this constant, or loop over the dictionary and configure different settings
for each simulation. In this script, we are looping over different force field parameters
3) Choose PACKMOL settings. If you want to use more membrane beads, don't forget to change the 'packmol_dict' to reflect this.
The final GMX topology is constructed based on the bead names in the "CELL.gro' that PACKMOL creates. 
4) Choose the mdp run settings of your simulations. By default, the standard Martini MDP settings are used. You can change
any setting that you like by changing it in a mdp dictionary. 

You can check "dict_examples.py" for basic examples on how to modify these configurations for multiple bead types. 
The current script restricts itself to single cells, with only a single membrane bead type
'''

from CELL_tools_library import DatabaseMaintenance
from CELL_tools_library import ForceFieldParser
from CELL_tools_library import InputSingleCell
from CELL_tools_library import PackmolExecuterSingleCell
from CELL_tools_library import Gromacs_IO
from CELL_tools_library import Gromacs_MDP
from CELL_tools_library import GromacsRun

import os
import pandas as pd
import numpy as np
import shutil
from itertools import product

##CONFIGURATION 1: Choose a descriptive job name
jobname = "This is an example jobname" 

##CONFIGURATION 2: Setting the range of search parameters for the force field definition

ff_dict = {
    'bead_types': {
        'M1': {'name': 'M1', 'mass': 72},
        'C': {'name': 'C', 'mass': 72}
    },
    'bond_types': {
        'C_M1': {'r0': 1.6, 'fk': 1000}
    },
    'nonbond_params': {
        'M1_M1': {'sigma': 3.5, 'epsilon': 0.5}
    }
}

r0_range = np.arange(1.2, 1.8, 0.2).round(2)
fk_range = np.arange(250, 2000, 250).round(2)
mass_range = np.arange(50, 100, 10).round(2)
sigma_range = [0.47] #let's keep sigma constant for now. 
epsilon_range = np.arange(1, 3, 1).round(2)

##CONFIGURATION 3: configuring PACKMOL settings
packmol_dict = {
        'settings': {
            'tolerance' : '4', #tolerance in angstrom
            'radius' : '18', #radius in angstrom
            'shape' : 'sphere' #cube, ellipsoid are supported, but untested
    },
        'beads': {
        'M1': {'name': 'M1', 'number': 180},
        'C': {'name': 'C', 'number': 1},
    }
}

##CONFIGURATION 4: using leapfrog and Berendsen thermostat for all simulations
mdp_settings= {
    'integrator': 'md',
    'tcoupl': 'Berendsen',
    'tc-grps': 'System',
    'tau_t': '1.0',
    'ref_t': '310',     
    'pcoupl': 'no'
}

# Before running, let's create a backup of our database in case we accidentally start running simulations
# we don't want to run. This can be restored in case of mistakes
shutil.copy2('simulations_CELL.csv', 'simulations_CELL.csv.old')

PackmolExecuterSingleCell.check_packmol_path()
Gromacs_IO.check_GMX_path()

#Create all the parameter combinations
param_comb = pd.DataFrame(list(product(r0_range, fk_range, mass_range, sigma_range, epsilon_range)),
                      columns=['r0_bond', 'fk_bond', 'mass', 'sigma', 'epsilon'])

# In the following part, we check the overlap of directory names
dir_names = []

for _, row in param_comb.iterrows():
    dir_name = "_".join([f"{key}" for key in ff_dict['bead_types'] if key.startswith('M')] +
                        [f"mass{row['mass']}_r{row['r0_bond']}_fk{row['fk_bond']}_sigma{row['sigma']}_eps{row['epsilon']}",
                         *[str(value) for key, value in mdp_settings.items() if key in ['integrator', 'tcoupl', 'pcoupl']]])
    
    dir_names.append(dir_name)

#ask the user to confirm the parameters of the runs
proceed = DatabaseMaintenance.check_repeat_sim('simulations_CELL.csv', dir_names)

if proceed:
    os.makedirs('tmp_run', exist_ok=True)
    os.chdir('tmp_run')
    #when choosing to proceed, we need to loop over both the parameters and the possible new directory names
    for (_, row), dir_name in zip(param_comb.iterrows(), dir_names):
        #Update the user on the simulation progress
        sim_counter = int(_,) + 1
        print("\nCurrently at sim " + str(sim_counter) + ". Dirname: '" + dir_name + "'")
    
        #we create a new dictionary for our parameter combinations
        ff_dict_param_comb = {
            'bead_types': {
                'M1': {'name': 'M1', 'mass': row['mass']},
                'C': {'name': 'C', 'mass': row['mass']}
            },
            'bond_types': {
                'C_M1': {'r0': row['r0_bond'], 'fk': row['fk_bond']}
            },
            'nonbond_params': {
                'M1_M1': {'sigma': row['sigma'], 'epsilon': row['epsilon']}
            }
        }
    
        parser = ForceFieldParser(ff_dict_param_comb)
        parser.generate_forcefield_itp()

        packmol = PackmolExecuterSingleCell(packmol_dict)
        packmol.run_packmol_single_CELL()

        Gromacs_IO.convert_xyz_to_gro("CELL.xyz", "CELL.gro", box_size = 5)
        Gromacs_IO.build_GMX_top_single_CELL('CELL.gro', 'forcefield.itp')

        Gromacs_MDP.write_min_mdp_file('min.mdp', 'cg', '5000', '5', '0.01') #do not need to change
        Gromacs_MDP.write_run_mdp_file('run.mdp', '5000', mdp_settings) #50000000 steps for 1 us

        GromacsRun.run_GMX_basic(filename='CELL.gro', timeout=60, nr_of_threads=12) #2400 is default for 1 us

        #Copy the files to their respective locations
        DatabaseMaintenance.move_files_to_parentdir('.', f'../../log/{dir_name}', ['.log'])
        DatabaseMaintenance.move_files_to_parentdir('.', f'../../params/{dir_name}', ['.top', '.itp'])
        DatabaseMaintenance.move_files_to_parentdir('.', f'../../traj/{dir_name}', ['.trr', '.edr', '.xtc'])
        DatabaseMaintenance.move_files_to_parentdir('.', f'../../settings/{dir_name}', ['.mdp', '.tpr'])
        DatabaseMaintenance.move_files_to_parentdir('.', f'../../coord/{dir_name}', ['.gro', '.xyz'])
    
        #Now we update the database with our simulation settings
        bead_names = [f"{key}" for key in ff_dict['bead_types'] if key.startswith('M')]
    
        entries_db = [{
        'job_name': jobname,
        'dir_name': dir_name,
        'crash_score': 0.0,
        'scoring_dict': {},
        'beads': bead_names,
        'mass': row['mass'],
        'LJ_params': {
            'sigma': row['sigma'],
            'epsilon': row['epsilon']
        },
        'integrator': mdp_settings.get('integrator', ''),
        'tcoupl': mdp_settings.get('tcoupl', ''),
        'pcoupl': mdp_settings.get('pcoupl', '')
        }]
    
        DatabaseMaintenance.save_new_entries_db('../simulations_CELL.csv', entries_db)
    os.chdir('../')     
    #do cleanup and get rid of any leftover files
    shutil.rmtree('tmp_run')
else:
    print("The simulation setup is stopped upon your request.")

### Todo:
### Run the simulations in a temporary directory, and ask the user to merge the requested simulations with the database. 
### Create a simulation and crash score in case the simulation already exists. Append '_repeatX' at the end of the dir name
### if it already exists and the user wants to proceed. 
### 