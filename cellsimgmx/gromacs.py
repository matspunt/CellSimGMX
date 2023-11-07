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

import logging 
import sys
import subprocess as sp
from distutils.spawn import find_executable

from cellsimgmx import CLIParser
from cellsimgmx import JSONParser
from cellsimgmx import SystemConstructor

class SimulationPreparation:
    """
    Ensures and checks all necessary files are there to run GMX simulations. Creates .mdp / runfiles based
    on input.JSON for the minimization, equilibration and production part of the simulation. 
    
    write_mdp_files():
        Creates mdp files. 
                
    construct_tpr():
        ???? --> does this belong here or in next class?    
    
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        self.system = SystemConstructor()
        self.gmx_path = None
        self.args = self.cli_parser.parse_args()
        #only run contents of this class if a simulation was requested
        if self.args.no_sim == True:
            logging.warning("You disabled simulations using the '--no-sim' flag, exiting now...")

        if self.args.no_sim == False:
            print("INFO: Starting simulation routine. ")
            logging.info(f"Succesfully made it to the simulation routine")
            self.gmx_path = find_executable("gmx") #check whether gmx is available in the path and print to the user.
            if self.gmx_path is not None:
                print(f"INFO: gromacs found at: '{self.gmx_path}'")
                logging.info(f"Found a gromacs executable at: '{self.gmx_path}'")
            else:
                logging.error("GMX executable not found. Please make sure it's installed and sourced in PATH as 'gmx'")
                sys.exit(1)
            #Execute the functions that setup simulation input!
            #self.create_mdp_files()
            #self.construct(tpr)
            
    def write_mdp_files(self, mdptype):
        """
            Contains all necessary logic to parse mdps. 
        """
        
        #standard mdp file
        standard_minimization_mdp = {
            'integrator': 'steep',
            'nsteps': '10000',
            'emtol': '10',
            'emstep': '0.01',
            'nstxout-compressed': '1000',
            'cutoff-scheme': 'Verlet',
            'vdw_type': 'cutoff',
            'vdw-modifier': 'Potential-shift-verlet',
            'rvdw': '1.1'
            #don't have electrostatics in system so can ignore all electrostatics .mdp settings
            }
        
        conjugate_descent_mdp = {
            'integrator': 'cg',
            'nsteps': '10000',
            'emtol': '10',
            'emstep': '0.001',
            'nstcgsteep': '100'
        }
        
        fep_minimization_mdp = {
            'free_energy': 'yes',
            'couple-lambda0': 'vdw',
            'vdw_lambdas': '1.0'
        }
        
        
class RunSimulation:
    """
    Executes the .tpr, tracking the simulation progress. Important progress information is printed to the terminal. 
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()

### How to handle minimization and equilibration??!