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

from cellsimgmx import CLIParser
from cellsimgmx import JSONParser

class SimulationPreparation:
    """
    Ensures and checks all necessary files are there to run GMX simulations. Creates .mdp / runfiles based
    on input.JSON for the minimization, equilibration and production part of the simulation. 
    
    create_mdp_files:
        Creates mdp files. 
        
    gather_simulation_files:
        Checks whether all necessary files are there, whether they are readable etc. 
        
    construct_tpr:
        ???? --> does this belong here or in next class?    
    
    Returns dict:
        Containing file names/paths of all files required to start the simulation. 
    
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        
class RunSimulation:
    """
    Executes the .tpr, tracking the simulation progress. Important progress information is printed to the terminal. 
    
    
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()