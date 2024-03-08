"""
Demonstration of wrapper to setup and run multiple simulations of different parameters with CellSimGMX. 
Needs to be executed inside a virtual environment with cellsimgmx installed ('which cellsimgmx')

v1.0
"""

import json
import numpy as np

####################################################################################################################
#                                           A. FORCE FIELD PARSER
#            Parses forcefield.itp based on a dictionary containing the bead names, mass and LJ identities of the particles
#                Note: currently the only hardcoded item is the name of the single center bead (hardcoded as 'C').
#                       Otherwise any arbitrary bead names and parameters are accepted by the class. 
#             Any number of membrane beads is allowed, assuming the LJ interaction pairs are defined in the dictionary. 
####################################################################################################################

class ForceFieldParser:
    """ 
        This class takes a dictionary including information about the force field (but is agnostic
        to the number of beads in the force field, as long as all the LJ and bonded terms are defined)
        For examples of force field dictionaries, check below. 
    """
    def __init__(self, dictionary):
        self.dictionary = dictionary
        self.itp_content = ""
        self.generate_forcefield_itp()

    def generate_atomtypes_section(self):
        """
        For code readability, sections are separated in helper functions which each take care of 
        one component of the 'forcefield.itp'
        """
        self.itp_content += "[ atomtypes ]\n"
        for bead_type, bead_info in self.dictionary['bead_types'].items():
            name = bead_info['name']
            mass = bead_info['mass']
            self.itp_content += f"   {bead_type.ljust(7)}    {mass}     0.0      A      0.0       0.0\n"

    def generate_nonbond_params_section(self):
        self.itp_content += "\n[ nonbond_params ]\n"
        nonbond_params = self.dictionary['nonbond_params']
        for key, values in nonbond_params.items():
            bead_type1, bead_type2 = key.split('_')  # Split over the underscore to obtain the bead types
            sigma = values['sigma']
            epsilon = values['epsilon']
            self.itp_content += f"  {bead_type1}  {bead_type2}    1      {sigma}  {epsilon}\n"

    def generate_bondtypes_section(self):
        self.itp_content += "\n[ bondtypes ]\n"
        for bond_type, bond_info in self.dictionary['bond_types'].items():
            r0 = bond_info['r0']
            fk = bond_info['fk']
            split_string = bond_type.split("_")  # Split over underscore to obtain bead types
            self.itp_content += f"   {split_string[0]}     {split_string[1]}  1     {r0}    {fk}\n"

    def generate_forcefield_itp(self):
        """
        Assembles the force field based on the helper functions and the provided dictionary
        """
        self.itp_content = "[ defaults ]\n"
        self.itp_content += "    1           2          yes          1.0     1.0\n\n"
        self.generate_atomtypes_section() #call helper functions to assemble force field
        self.generate_nonbond_params_section()
        self.generate_bondtypes_section()

        with open("forcefield.itp", "w") as itp:
            itp.write(self.itp_content)
            
# specifying the basic forcefield. This is a template force field for the simplest possible model of a CELL
# (a single center bead and one type of membrane particle). It is not necessary to change this
ff_dict_single_bead_type = {
    'bead_types': {
        'M1': {'name': 'M1', 'mass': 72},
        'C': {'name': 'C', 'mass': 72}
    },
    'bond_types': {
        'C_M1': {'r0': 1.85, 'fk': 500}
    },
    'nonbond_params': {
        'M1_M1': {'sigma': 0.47, 'epsilon': 2.5}
    }
}

# template forcefield for two hypothetical membrane particles. 
ff_dict_two_bead_types = {
    'bead_types': {
        'M1': {'name': 'M1', 'mass': 72},
        'M2': {'name': 'M2', 'mass': 72},
        'C': {'name': 'C', 'mass': 72},
    },
    'bond_types': {
        'C_M1': {'r0': 1.85, 'fk': 500},
        'C_M2': {'r0': 1.90, 'fk': 500}
    },
    'nonbond_params': {
        'M1_M1': {'sigma': 0.47, 'epsilon': 2.0},
        'M2_M2': {'sigma': 0.47, 'epsilon': 3.0},
        'M1_M2': {'sigma': 0.47, 'epsilon': 2.0}
    }
}

####################################################################################################################
#                                           B. JSON PARSER
#         Takes as input a JSON file, reads it and allows you to modify it based on a dict of JSON settings.                
####################################################################################################################

class JSONParser:
    """Overwrites input.JSON file with new parameters
    """
    def __init__(self, json_path, modified_values=None, output_path = "input.json"):
        self.json_path = json_path
        self.output_path = output_path
        self.load_json()

        if modified_values:
            self.update_values(modified_values)
            
        self.save_changes()

    def load_json(self):
        with open(self.json_path, 'r') as json_file:
            self.data = json.load(json_file)

    def update_values(self, modified_values):
        # self.data is formatted as [category][param]["value"] while modified_values is given as [param]["value"]
        for category in self.data:
            for param, value in modified_values.items():
                if param in self.data[category]:
                    self.data[category][param]["value"] = value

    def save_changes(self):
        with open(self.output_path, 'w') as json_file:
            json.dump(self.data, json_file, indent=2)
            

####################################################################################################################
#                                           C. Configuring the simulations you want to run
#                                           Configuring input for the wrapper               
####################################################################################################################

#link to a template JSON file. You don't need to change any parameters, if you don't specify changes, simply the 
# file you link here will be used for all simulations
JSON_PATH = "/wrk/matspunt/coding/CellSimGMX/cellsimgmx/ff/input.json"

# 1. Modifying JSON options: this can be done through a single dictionary, overwriting your template for EACH simulation
modified_json_options = {'nr_of_particles': 200, 'nearest_neighbour_springs': 4}

# or, you can do this in a loop for multiple values. For instance specifying 10 different options for 'nr_of_particles'
nr_of_particles = np.arange(100, 1100, 100)

for i in nr_of_particles:
    modified_json_options['nr_of_particles'] = i
    print(modified_json_options)

#2. Modifying force field options. For instance, specifying 10 different masses from 10 - 100 
mass_increments = np.arange(10, 101, 10)

for mass in mass_increments:
    ff_dict_single_bead_type['bead_types']['M1']['mass'] = mass
    #do something with the dict
    print(ff_dict_single_bead_type)

parser = ForceFieldParser(ff_dict_single_bead_type)
json_parser = JSONParser(JSON_PATH, modified_json_options, output_path = "input.json")