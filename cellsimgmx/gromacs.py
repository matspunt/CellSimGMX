# -*- coding: utf-8 -*-
# Copyright (C) 2024 Mats Punt mats.punt(at)helsinki.fi
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
    
    write_mdp_minim():
        Creates mdp file of minimization routine based on input.JSON
        
    write_mdp_eq_prod():
        Creates mdp file of equilibration and production routine based on input.JSON 
    
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        self.json_values = self.json_parser.json_values
        self.system = SystemConstructor()
        self.gmx_path = None
        self.args = self.cli_parser.parse_args()
        # only run contents of this class if a simulation was requested!
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
                        
        #execute functions
        self.write_mdp_minim()
        self.write_mdp_eq_prod()
        
    def write_mdp_minim(self):
        """
            Simple function to construct mdp of minimization routine, based on JSON input
        """
        
        minimization = self.json_values["minimization"]
        
        steep_mdp = {
            #basic minimization routine
            'integrator': 'steep',
            'nsteps': '10000',
            'emtol': '10',
            'emstep': '0.01',
            'nstxout-compressed': '1000',
            'cutoff-scheme': 'Verlet',
            'vdw_type': 'cutoff',
            'vdw-modifier': 'Potential-shift-verlet',
            'rvdw': '1.1',
            'rcoulomb': '1.1',
            #don't have electrostatics in system so can ignore all electrostatics .mdp settings
            }
        
        cg_mdp = {
            # more stringent minimization to deal with poor initial configurations
            'integrator': 'cg',
            'nsteps': '10000',
            'emtol': '5',
            'emstep': '0.01',
            'nstxout-compressed': '1000',
            'cutoff-scheme': 'Verlet',
            'vdw_type': 'cutoff',
            'vdw-modifier': 'Potential-shift-verlet',
            'rvdw': '1.1',
            'rcoulomb': '1.1',
            #don't have electrostatics in system so can ignore all electrostatics .mdp settings     
        }
        
        em_dict = steep_mdp if minimization == "steep" else cg_mdp

        with open(f"{self.args.output_dir}/mdps/em.mdp", mode='w') as mdp:
            mdp.write(f"; {minimization.capitalize()} descent minimization\n")
            for key, value in em_dict.items():
                mdp.write(f"{key} = {value}\n")

    def write_dict_mdp_helper(self,dict, mdp_name):
        """Helper function to write dictionary contents to MDP file

        Args:
            dictionary (dict): Dict containing MDP options
            filename (str): Name of the mdp_file you want to write to
        """
        with open(mdp_name, "w") as file:
            for key, value in dict.items():
                file.write(f"{key} = {value}\n")
        
    def write_mdp_eq_prod(self):
        """
        Writes equilibration and production MDP files based on requested simulation routine
        """
        
        ensemble = self.json_values["ensemble"]
        timestep = self.json_values["timestep"] #sets production only!
        number_of_steps = self.json_values["number_of_steps"]
        number_of_steps = tuple(map(int, number_of_steps.split(', '))) # [0] is equilibration nsteps, [1] is production nsteps
                
        def extend_mdp(base_mdp, extra_mdp_options):
            return {**base_mdp, **extra_mdp_options}

        NVE_mdp = {
            #basic NVE mdp, only options which are different from default GMX settings are included
            'integrator': 'md',
            'dt': '0.02', # timestep for equilibration
            'nsteps': number_of_steps[0], # sets number of steps for equilibration
            'nstxout-compressed': '1000',
            'compressed-x-precision': '100',
            'nstlist': '20',
            'ns_type': 'grid',
            'pbc': 'xyz',
            'vdw_type': 'cutoff',
            'vdw-modifier': 'Potential-shift-verlet',
            'rvdw': '1.1',
            'rcoulomb': '1.1',
            'constraints': 'none'
            # tcoupl and pcoupl are disabled by default (don't need to specify)
        }
        
        self.write_dict_mdp_helper(NVE_mdp, "mdps/NVE_eq.mdp")

        #EQUILIBRATION with Berendsen coupling
        NVT_eq_mdp = extend_mdp(NVE_mdp, {'tcoupl': 'berendsen', 'tc-grps': 'System', 'tau-t': '1.0', 'ref-t': '310'})
        NpT_eq_mdp = extend_mdp(NVT_eq_mdp, {'pcoupl': 'berendsen', 'pcoupltype': 'isotropic', 'tau_p': '2.0', 'ref_p': '1.0', 'compressibility': '4.5e-05', ';refcoord-scaling': 'com'}) #refcoord-scaling is required for simulations with pressure coupling and position restraints

        #PRODUCTION with v-rescale thermostat and c-rescale barostat
        NVE_mdp['dt'] = timestep
        NVE_mdp['nsteps'] = number_of_steps[1]
        NVT_prod_mdp = extend_mdp(NVE_mdp, {'tcoupl': 'v-rescale', 'tc-grps': 'System', 'tau-t': '1.0', 'ref-t': '310'})
        NpT_prod_mdp = extend_mdp(NVT_prod_mdp, {'pcoupl': 'C-rescale', 'pcoupltype': 'isotropic', 'tau_p': '5.0', 'ref_p': '1.0', 'compressibility': '4.5e-05', ';refcoord-scaling': 'com'})
        
        #Save required mdps dependent on requested ensemble
        if ensemble == "NVE":
            NVE_mdp['dt'] = timestep
            NVE_mdp['dt'] = number_of_steps[1]
            self.write_dict_mdp_helper(NVE_mdp, "mdps/production.mdp")

        elif ensemble == "NVT":
            self.write_dict_mdp_helper(NVT_eq_mdp, "mdps/NVT_eq.mdp")
            self.write_dict_mdp_helper(NVT_prod_mdp, "mdps/production.mdp")

        elif ensemble == "NpT":
            self.write_dict_mdp_helper(NVT_eq_mdp, "mdps/NVT_eq.mdp")
            self.write_dict_mdp_helper(NpT_eq_mdp, "mdps/NpT_eq.mdp")
            self.write_dict_mdp_helper(NpT_prod_mdp, "mdps/production.mdp")

class RunSimulation:
    """
    Generates .tpr and runs the simulation based on mdps that are found (i.e. the required simulations are inferred from the mdps that are saved)
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        
        # based on MDP files that are found, generate .tpr