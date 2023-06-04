'''Examples of dictionaries that can be used to configure simulations using 'CELL_tools.py'
'''

####################################################################################################################
#                                           FORCE FIELD SETTINGS
#       Note: works with any number of bead types, as long as all the bonded types and LJ params are in the dict
#       The force field dict can also easily be modified by a for loop.
####################################################################################################################

dict_two_bead_types = {
    'bead_types': {
        'M1': {'name': 'M1', 'mass': 72},
        'M2': {'name': 'M2', 'mass': 72},
        'C': {'name': 'C', 'mass': 72},
    },
    'bond_types': {
        'C_M1': {'r0': 1.6, 'fk': 1000},
        'C_M2': {'r0': 1.7, 'fk': 2000}
    },
    'nonbond_params': {
        'M1_M1': {'sigma': 0.47, 'epsilon': 2},
        'M2_M2': {'sigma': 0.47, 'epsilon': 2},
        'M1_M2': {'sigma': 0.47, 'epsilon': 2}
    }
}

dict_four_bead_types = {
    'bead_types': {
        'M1': {'name': 'M1', 'mass': 72},
        'M2': {'name': 'M2', 'mass': 72},
        'M3': {'name': 'M3', 'mass': 72},
        'M4': {'name': 'M4', 'mass': 72},
        'C': {'name': 'C', 'mass': 72}
    },
    'bond_types': {
        'C_M1': {'r0': 1.6, 'fk': 1000},
        'C_M2': {'r0': 1.7, 'fk': 2000},
        'C_M3': {'r0': 1.8, 'fk': 3000}
    },
    'nonbond_params': {
        'M1_M1': {'sigma': 0.47, 'epsilon': 2},
        'M2_M2': {'sigma': 0.47, 'epsilon': 2},
        'M1_M2': {'sigma': 0.47, 'epsilon': 2},
        'M1_M3': {'sigma': 0.47, 'epsilon': 2},
        'M2_M3': {'sigma': 0.47, 'epsilon': 2},
        'M3_M3': {'sigma': 0.47, 'epsilon': 2},
        'M1_M4': {'sigma': 0.47, 'epsilon': 2},
        'M2_M4': {'sigma': 0.47, 'epsilon': 2},
        'M3_M4': {'sigma': 0.47, 'epsilon': 2},
        'M4_M4': {'sigma': 0.47, 'epsilon': 2}
    }
}

dict_single_bead_type = {
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

####################################################################################################################
#                                           PACKMOL SETTINGS
#              Supports any number of bead types. Note that the tolerance (4) and radius (18)
#              have been finetuned for 180 beads in total in the membrane. If you want to use more
#              or less, you might have to play around with these settings.
#              Currently only works with a single center bead!
####################################################################################################################

packmol_dict = {
        'settings': {
            'tolerance' : '4', #tolerance in angstrom
            'radius' : '18', #radius in angstrom
            'shape' : 'sphere' #cube and ellipsoid are also supported, but untested 
    },
        'beads': {
        'M1': {'name': 'M1', 'number': 45},
        'M2': {'name': 'M2', 'number': 45},
        'M3': {'name': 'M3', 'number': 45},
        'M4': {'name': 'M4', 'number': 45},
        'C': {'name': 'C', 'number': 1},
    }
}

####################################################################################################################
#                                               MDP OPTIONS
#                   Some examples of MDP options you might want to use. See also:
#                   https://manual.gromacs.org/current/user-guide/mdp-options.html
####################################################################################################################

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

