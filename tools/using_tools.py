#You need to run this script in the directory where you want to run the test simulations
# make sure "cell_layer_tools.py" is in the same directory
# A 'CELL.gro' file needs to be in that directory too. ensure the atom names are
# 'M1', 'M2' and 'C' as these are read from the force field. Make sure it has valid 
# box dimensions

from cell_layer_tools import Gromacs_MDP
from cell_layer_tools import Gromacs_TOP
from cell_layer_tools import GromacsRun

import os
import shutil
import numpy as np
import pandas as pd
from itertools import product #for Cartesian product

#some possible MDP options we could use
langevin = {
    'integrator': 'sd',
    'tcoupl': 'no',
    'tau_t': '1.0',
    'ref_t': '298',
    'ld-seed' : '-1',            
}

berendsen_NVT = {
    'tcoupl': 'Berendsen',
    'tc-grps': 'System',
    'tau_t': '1.0',
    'ref_t': '298',
}

berendsen_NPT = {
    'tcoupl': 'Berendsen',
    'tc-grps': 'System',
    'tau_t': '1.0',
    'ref_t': '298',
    'Pcoupl': 'Berendsen',
    'Pcoupltype': 'isotropic',
    'tau_p': '12.0',
    'compressibility': '3e-4',
    'ref_p': '1',    
}

gen_vel_disabled = {
    'gen_vel': 'no',
    'gen_seed': '-1'
}

"""
forcefield_two_bead_types = {
    'M1_mass': '72',
    'M2_mass': '72',
    'C_mass': '72',
    'M1_M1_sigma': '0.47',
    'M1_M1_epsilon': '2',
    'M1_M2_sigma': '0.47',
    'M1_M2_epsilon': '2',
    'M2_M2_sigma': '0.47',
    'M2_M2_epsilon': '2',
}
"""

#a basic force field description
forcefield_one_bead_type = {
    'M1_mass': '72',
    'C_mass': '72',
    'M1_M1_sigma': '0.47',
    'M1_M1_epsilon': '2',
}

#a basic bond description

bond_dict = {
    ('C', 'M1'): {'func': 1, 'r0': 1.6, 'fk': 5000}
}

#We can set combinations of bonded parameters or LJ settings in the following manner. 
# or set it simply to a constant value if we want it constant. 
fk_values_m1 = np.linspace(1000, 10000, 10)
r0_values_m1 = np.linspace(1.0, 2.0, 10)
sigma_m1_m1 = np.linspace(0.1, 1.0, 10)
epsilon_m1_m1 = np.linspace(1.0, 5.0, 10)

combinations = []

for fk in fk_values_m1:
    for r0 in r0_values_m1:
        for sigma in sigma_m1_m1:
            for epsilon in epsilon_m1_m1:
                combination = {
                    'fk': fk,
                    'r0': r0,
                    'M1_M1_sigma': sigma,
                    'M1_M1_epsilon': epsilon
                }
                combinations.append(combination)

df = pd.DataFrame(combinations)
nr_comb = len(df)
print("Parameter combinations:", nr_comb)

for index, row in df.iterrows():
    folder_name = f"sim-{index + 1}"
    os.makedirs(folder_name)
    os.chdir(folder_name)
    print(os.getcwd())  # Print the current simulation directory
    os.system("cp ../CELL.gro .")
    new_forcefield_dict = {
        'M1_mass': str(forcefield_one_bead_type['M1_mass']),
        'C_mass': str(forcefield_one_bead_type['C_mass']),
        'M1_M1_sigma': str(row['M1_M1_sigma']),
        'M1_M1_epsilon': str(row['M1_M1_epsilon'])
    }

    new_bond_dict = {
        ('C', 'M1'): {
            'func': '1',
            'r0': str(row['r0']),
            'fk': str(row['fk'])
        }
    }

    #in each folder, write the required files for GROMACS
    Gromacs_MDP.write_min_mdp_file('min.mdp', 'cg', '5000', '5', '0.01') 
    Gromacs_MDP.write_run_mdp_file('run.mdp', '5000', langevin) #50000000 for 1 us

    Gromacs_TOP.parse_forcefield_one_bead_type(new_forcefield_dict)
    Gromacs_TOP.build_GMX_itp_from_gro("CELL.gro", new_bond_dict)
    Gromacs_TOP.parse_topol(nr_of_cells='1')
    
    # And execute the simulation based on the saved .mdp settings. 
    GromacsRun.run_GMX_basic(filename='CELL.gro', timeout=15, nr_of_threads=12, top_name='topol.top') #2400 is default for 1 us

    os.chdir('..')

