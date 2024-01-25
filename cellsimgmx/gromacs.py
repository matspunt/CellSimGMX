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
        
    Todo: add JSON parsing options of Simulation header
    
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
            
    def build_mdp(self):
        """
            Function to construct mdp objects based on JSON information. 
        """
        
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
            'rvdw': '1.1'
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
            'rvdw': '1.1'
            #don't have electrostatics in system so can ignore all electrostatics .mdp settings     
        }

        #assert Type in ['min', 'NVE_eq', 'NVT_eq', 'NPT_eq', 'prod'], "type should be be EM, NVT, NPT, MD"
        
        def addMDPOption(name: str, value, comment=''):
            """Formatting function for adding a parameter.

            Args:
                name (str): parameter name.
                value (any): parameter value.
                comment (str, optional): parameter inline comment. Defaults to ''.
            """

            if comment == '':
                file.write("{:20s} = {:13s}\n".format(name, str(value)))
            else:
                file.write("{:20s} = {:13s} ; {:13s}\n".format(name, str(value), comment))

        if Type in ['EM']:
            dt = 0.01
            addParam('integrator', 'steep', 'Use steep for EM.')
            addParam('emtol', 1000, 'Stop when max force < 1000 kJ/mol/nm.')
            addParam('emstep', dt, 'Time step (ps).')

        if Type in ['NVT', 'NPT', 'MD']:
            dt = 0.002
            addParam('integrator', 'md')
            addParam('dt', dt, 'Time step (ps).')

        addParam('nsteps', nsteps, "{:d} ps.".format(int(dt * nsteps)))

        # OUTPUT

        addTitle("Output control")
        addParam('nstxout-compressed', nstxout, 'Write .xtc frame every {:d} ps.'.format(int(dt * nstxout)))

        # NEIGHBOUR SEARCHING

        addTitle("Neighbour searching")
        addParam('cutoff-scheme', 'Verlet', 'Related params are inferred by GROMACS.')

        # BONDED

        if Type in ['NVT', 'NPT', 'MD']:
            addTitle("Bond parameters")
            addParam('constraints', 'h-bonds', 'Constrain H-bond vibrations.')
            addParam('constraint_algorithm', 'lincs', 'Holonomic constraints.')
            addParam('lincs_iter', 1, 'Related to accuracy of LINCS.')
            addParam('lincs_order', 4, 'Related to accuracy of LINCS.')

        # ELECTROSTATICS

        addTitle("Electrostatics")
        addParam('coulombtype', 'PME', 'Particle Mesh Ewald electrostatics.')
        addParam('rcoulomb', 1.2, 'CHARMM is calibrated for 1.2 nm.')
        addParam('fourierspacing', 0.14)

        # VAN DER WAALS

        addTitle("Van der Waals")
        addParam('vdwtype', 'cut-off', 'Twin range cut-off with nblist cut-off.')
        addParam('rvdw', 1.2, 'CHARMM is calibrated for 1.2 nm.')
        addParam('vdw-modifier', 'force-switch', 'Specific for CHARMM.')
        addParam('rvdw-switch', 1.0, 'Specific for CHARMM.')

        # TEMPERATURE COUPLING

        if Type in ['NVT', 'NPT', 'MD']:
            addTitle("Temperature coupling")
            addParam('tcoupl', 'v-rescale')
            addParam('tc-grps', 'SYSTEM')
            addParam('tau-t', 0.5, 'Coupling time (ps).')
            addParam('ref-t', 300, 'Reference temperature (K).')

        # PRESSURE COUPLING

        if Type in ['NPT', 'MD']:
            addTitle('Pressure coupling')
            addParam('pcoupl', 'C-rescale', 'Use C-rescale barostat.')
            addParam('pcoupltype', 'isotropic', 'Uniform scaling of box.')
            addParam('tau_p', 5.0, 'Coupling time (ps).')
            addParam('ref_p', 1.0, 'Reference pressure (bar).')
            addParam('compressibility', 4.5e-05, 'Isothermal compressbility of water.')

            if Type == 'NPT' or posRes:
                addParam('refcoord_scaling', 'com', 'Required with position restraints.')

        # PERIODIC BOUNDARY CONDITIONS

        addTitle("Periodic boundary condition")
        addParam('pbc', 'xyz', 'Apply periodic boundary conditions.')

        # WRAP UP

        file.close()

        
        with open("param.mdp", 'w') as mdp:
            mdp.write(standard_minimization_mdp)
            #format MDP correctly, not simple dict > str conversion!!
        mdp.close()
                
class RunSimulation:
    """
    Executes the .tpr, tracking the simulation progress. Important progress information is printed to the terminal. 
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()

### How to handle minimization and equilibration??!