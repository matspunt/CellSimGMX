'''06/06/23. Template script for setting up simulations with the CELL tools code in GROMACS.

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

DO NOT FORGET SETTING THE TIMEOUT, AND NSTEPS!!!!!!!!!!!!!!!!!!!!!! IN THE MAIN GROMACS.RUN FUNCTION!!!!!! DESCRIBE THIS ON THE GITHUB
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
jobname = "Populating CSV. Random parameter search" 

##CONFIGURATION 2: Setting the force field and the range of parameters (can keep constant if you want to test .mdp settings)

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
#fk_range = np.arange(250, 2000, 250).round(2)
#mass_range = np.arange(50, 100, 25).round(2)
#mass_range = [72]
#sigma_range = [0.47] #let's keep sigma constant for now. 
#epsilon_range = np.arange(1, 3, 1).round(2)

#Example for demonstrating script during presentation
fk_range = [250]
mass_range = [72]
sigma_range = [0.47] #let's keep sigma constant for now. 
epsilon_range = [1]

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

#Let's spawn a temporary folder and create a local copy of the database that we will work with
os.makedirs('.tmp_run', exist_ok=True)
os.chdir('.tmp_run')

shutil.copy2('../../simulations_CELL.csv', 'simulations_CELL.csv')

#Then ask  the user to confirm
proceed = DatabaseMaintenance.check_repeat_sim('simulations_CELL.csv', dir_names)

if proceed:
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
        nsteps = Gromacs_MDP.write_run_mdp_file('run.mdp', '50000', mdp_settings) # here we set the number of steps explicitly. 50000000 steps = 1 us
        time_per_sim = GromacsRun.calculate_time_per_sim(nsteps)
        print(f"This simulation will take roughly {time_per_sim} minutes")
        
        crash_status = GromacsRun.run_GMX_basic(filename='CELL.gro', timeout=15, nr_of_threads=12) #timeout in seconds. Default is 2400 s for a 1 us simulation

        #First, we move the simulation files to temporary folders in '.tmp_run'
        DatabaseMaintenance.move_files_by_ext('.', f'log/{dir_name}/.', ('.log'))
        DatabaseMaintenance.move_files_by_ext('.', f'params/{dir_name}/.', ('.top', '.itp'))
        DatabaseMaintenance.move_files_by_ext('.', f'traj/{dir_name}/.', ('.trr', '.edr', '.xtc'))
        DatabaseMaintenance.move_files_by_ext('.', f'settings/{dir_name}/.', ('.mdp', '.tpr'))
        DatabaseMaintenance.move_files_by_ext('.', f'coord/{dir_name}/.', ('.gro', '.xyz'))
    
        #Then we update the database in '.tmp_run'
        bead_names = [f"{key}" for key in ff_dict['bead_types'] if key.startswith('M')]
        entries_db = [{
        'job_name': jobname,
        'dir_name': dir_name,
        'beads': bead_names,
        'mass': row['mass'],
        'LJ_params': {
            'sigma': row['sigma'],
            'epsilon': row['epsilon']
        },
        'integrator': mdp_settings.get('integrator', ''),
        'tcoupl': mdp_settings.get('tcoupl', ''),
        'pcoupl': mdp_settings.get('pcoupl', ''),
        'crash': crash_status
        }]
        
        #Update the local database in '.tmp_run'
        DatabaseMaintenance.save_new_entries_db('simulations_CELL.csv', entries_db) 
else:
    print("\nThe simulation setup was stopped upon your request.")
    os.chdir('../')
    shutil.rmtree('.tmp_run')
    
#When the simulations are done, ask whether the user is really sure that they want to merge this version of simulations with the existing database
DatabaseMaintenance.confirm_and_cleanup('simulations_CELL.csv')
#if yes, update all the simulations with the database. Otherwise delete all the simulation files and restore the old database. 