# -*- coding: utf-8 -*-
# Copyright (C) 2023  Mats Punt mats.punt(at)helsinki.fi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

import argparse
import logging
import json
import glob
import sys
import os

class CLIParser:
    """
    This class accepts command line arguments of the programme and stores them. Can be passed on to other classes to access CLI strings given by user. 
    """
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="CellSimGMX: A 3D discrete element framework simulator for epithelial tissues using GROMACS")
        
        #specifying directories for programme
        self.parser.add_argument("--ff-dir", "-ff", help="Absolute path to the directory where the forcefield is stored (accepts arbitrary .itp file names)")
        self.parser.add_argument("--input-dir", "-in", help="Absolute path to the directory where the simulation input is stored (accepts arbitrary .json file name)")
        self.parser.add_argument("--output-dir", "-out", help="Absolute path where files should be generated and simulation should be run")
        
        #Miscellaneous
        self.parser.add_argument("--verbose", "-v", action="store_true", help="Optional argument. When enabled, prints detailed logging. Useful for debugging problems. ")

    def parse_args(self):
        return self.parser.parse_args()

class JSONParser:
    """
    This class handles parsing of the input.JSON file. An instance of this class can be passed to other classes to access the 
    simulation settings. Currently supports four classes of settings "Cell", "Tissue", "Matrix" and "Simulation", as seen in the example input.JSON. 
    
    Attributes:
        cli_parser (CLIParser): An instance of the CLIParser class
        json_path (str): The path to the JSON input file.
        data (dict): The loaded JSON data.
    """

    _instance = None

    # Class holds data so use a Singleton pattern to ensure only a single instance is created
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.json_values = None 
        return cls._instance
    
    def __init__(self):
        self.cli_parser = CLIParser()
        #these do not need to be in the Singleton pattern as we only use them to obtain 'json_values'
        self.json_path = None
        self.data = None
        
    def load_json(self):
        """
        Loads the input.json as data, verifies input complies with JSON standard. 
        """
        args = self.cli_parser.parse_args()
        if args.input_dir:
            if not os.path.exists(args.input_dir):
                logging.error(f"Error: Input dir '{args.input_dir}' does not exist. Please check for typos and try again")
                sys.exit(1)
                
            if not os.access(args.input_dir, os.R_OK):
                logging.error(f"Error: Input directory '{args.input_dir}' is not readable. Do you have permissions and is the directory correct?")
                sys.exit(1)

        else:
            logging.error("Error: Input dir for settings not specified. Please specify using '--input-dir' and try again")
            sys.exit(1)

        json_files = glob.glob(f"{args.input_dir}/*.json")
        
        if len(json_files) > 1:
            logging.error("Multiple .JSON files were found in your input dir which is causing conflicts. Remove the incorrect ones and try again: \n") 
            for i in json_files:
                logging.error(str(i))
            sys.exit(1)
            
        if len(json_files) > 0: 
            self.json_path = json_files[0]
            try:
                with open(self.json_path, 'r') as file:
                    self.data = json.load(file)
                    logging.info(f"Using simulation input from: {self.json_path}")
            except json.JSONDecodeError:
                logging.error("The supplied JSON file contains formatting errors, please check it complies with the template supplied with the CellSimGMX package and correct any typos")
                sys.exit(1)
            except FileNotFoundError:
                logging.warning("No input.JSON could be found in your '--input-dir' entry. Are you in the right directory and is your file extension correct?")
                sys.exit(1)
        else:
            logging.warning("No input.JSON could be found in your '--input-dir' entry. Are you sure you supplied a correct path?")
            sys.exit(1)

    def extract_json_values(self):
        """
        Extracts values from input.JSON and stores them in variables. Variable names are by construction IDENTICAL to JSON variable names. This should make working with them straightforward. 
        Does a few basic sanity checks on input (can be expanded later if needed). 
        """
        
        args = self.cli_parser.parse_args()
        
        if self.data is None:
            logging.error("Some unknown problem occurred reading the input.JSON. Terminating programme...")
            sys.exit(1)

        json_values = {}
        for category, category_items in self.data.items():
            for key, value_data in category_items.items():
                value = value_data.get("value")
                if value is not None:
                    var_name = f"{key}".lower()
                    json_values[var_name] = value

                    # Do some basic checks on the input to guide the user
                    if var_name == "number_of_cells":
                        if not (isinstance(value, int) and 1 <= value <= 2000):
                            logging.warning(f"You have specified {value} cells in the system, this is possible but will significantly slow down simulations!")

                    if var_name == "initial_packing_shape":
                        if not (isinstance(value, str) and value in ['spherical', 'cuboid', 'ellipsoid']):
                            logging.error(f"{value} is not a valid cell shape!")
                            sys.exit(1)

                    if var_name == "tissue_packing":
                        if not (isinstance(value, str) and value in ['hexagonal', 'monolayer']):
                            logging.error(f"{value} is not a valid packing for tissues!")
                            sys.exit(1)
                            
                    ###### May be extended if needed!!!!
                    ###### Correct later for final JSON file!
                            
        if args.verbose:
            print(f"\nVERBOSE MODE ENABLED. Parameters read from: {self.json_path}.") 
            print("These are the parameters you have configured\n")
            for var_name, value in json_values.items():
                print(f"{var_name} : {value}")

        # Always save parameter information to the logfile
        logging.info("CellSimGMX has found the following parameters to be configured:\n")
        for var_name, value in json_values.items():
            logging.info(f"{var_name} : {value}")

        self.json_values = json_values #store the JSON values in the Class

class ForcefieldParserGMX:
    """
    This class handles parsing of the predefined force field in GROMACS format and stores it in dicts

    Attributes:
        atomtypes (dict)
        nonbond_params (dict)
        bondtypes (dict)
    """
    
    #Use a Singleton pattern too here such that only one instance is created that is recycled
    _instance = None 

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(ForcefieldParserGMX, cls).__new__(cls)
            cls._instance.atomtypes = {}
            cls._instance.nonbond_params = {}
            cls._instance.bondtypes = {}
        return cls._instance
    
    def __init__(self):
        self.cli_parser = CLIParser()
        self.itp_path = None

    def parse_GMX_ff(self):
        """
            Looks for a GMX compatible .itp file in the user-specific directory. Then uses helper functions to parse the 
            different .itp entries into dictionaries. 
        """
        
        args = self.cli_parser.parse_args()

        if args.ff_dir:
            if not os.path.exists(args.ff_dir):
                logging.error(f"Force field directory '{args.ff_dir}' does not exist. Please check for typos and try again")
                sys.exit(1)
            
            if not os.access(args.ff_dir, os.R_OK):
                logging.error(f"Force field directory '{args.ff_dir}' is not readable. Do you have permissions and is the directory correct?")
                sys.exit(1)
        else:
            logging.error("Force field directory not specified. Please specify using '--ff-dir' and try again ")
            sys.exit(1)

        itp_files = glob.glob(f"{args.ff_dir}/*.itp")

        if len(itp_files) > 1:
            logging.error("Multiple .itp files were found in the force field directory which is causing conflicts. Aborting run.")
            print("Remove the incorrect files and please try again:")
            for i in itp_files:
                print(str(i))
            sys.exit(1)

        if len(itp_files) > 0:
            self.itp_path = itp_files[0]
            try:
                with open(self.itp_path, 'r') as file:
                    self.itp_data = file.read()
                logging.info(f"Using force field input from: {self.itp_path}")
            except FileNotFoundError:
                logging.error("No .itp files could be found in your '--ff-dir' entry. Are you in the right directory and is your file extension correct?")
        else:
            logging.error("No .itp files could be found in your '--ff-dir' entry. Are you in the right directory and is your file extension correct?")

        # open the force field .itp and extract the different '[ ]' categories
        # it makes little sense to check the formatting here ( I would also not know how to do it systematically?) as the user 
        # does not need to edit the force field themselves
        with open(self.itp_path, 'r') as itp:
            section = None
            for line in itp:
                line = line.strip()
                if line.startswith('[') and line.endswith(']'):
                    section = line[1:-1].strip()
                elif section == 'atomtypes':
                    self.parse_atomtype_GMX_helper(line)
                elif section == 'nonbond_params':
                    self.parse_nonbond_params_GMX_helper(line)
                elif section == 'bondtypes':
                    self.parse_bondtype_GMX_helper(line)
        
        if args.verbose:
            atom_names = list(self.atomtypes.keys())
            print(f"\nVERBOSE MODE ENABLED. Succesfully read the force field from {self.itp_path}. ")
            print("Your system contains the following particle types:")
            print(atom_names)
            
            print("\nAnd the following LJ interaction pairs:")
            interaction_pairs = list(self.nonbond_params.keys())
            print(interaction_pairs)
            
            print("\nAnd the following bonded types:")
            bonded_types = list(self.bondtypes.keys())
            print(bonded_types)

        logging.info("CellSimGMX has found the following atom types in the specified force field:\n")
        atom_names = list(self.atomtypes.keys())
        logging.info(atom_names)
        
        logging.info("And the following LJ interaction pairs:")
        interaction_pairs = list(self.nonbond_params.keys())
        logging.info(interaction_pairs)
        
        logging.info("And the following bonded types:")
        bonded_types = list(self.bondtypes.keys())
        logging.info(bonded_types)

    def parse_atomtype_GMX_helper(self, line):
        """
            Parses atom types from the .itp so .gro files can be built from them. 
        """
        if not line.startswith(';'):
            columns = line.split()
            if len(columns) >= 6:
                atomtype = columns[0]
                mass = float(columns[1])
                charge = float(columns[2])
                ptype = columns[3]
                sigma = float(columns[4])
                epsilon = float(columns[5])
                self.atomtypes[atomtype] = {'mass': mass, 'charge': charge, 'ptype': ptype, 'sigma': sigma, 'epsilon': epsilon}

    def parse_nonbond_params_GMX_helper(self, line):
        #skip lines with comments
        if not line.startswith(';'):
            columns = line.split()
            if len(columns) >= 5:
                i = columns[0]
                j = columns[1]
                func = int(columns[2])
                sigma = float(columns[3])
                epsilon = float(columns[4])
                self.nonbond_params[(i, j)] = {'func': func, 'sigma': sigma, 'epsilon': epsilon}

    def parse_bondtype_GMX_helper(self, line):
        #skip lines with comments
        if not line.startswith(';'):
            columns = line.split()
            if len(columns) >= 5:
                i = columns[0]
                j = columns[1]
                func = int(columns[2])
                r0 = float(columns[3])
                fk = float(columns[4])
                self.bondtypes[(i, j)] = {'func': func, 'r0': r0, 'fk': fk}