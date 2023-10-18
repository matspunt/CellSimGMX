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

import numpy as np

from cellsimgmx import CLIParser
from cellsimgmx import JSONParser
from cellsimgmx import CellTopology

class TissueConstructor:
    """
    This class uses the individual CELL-{timestamp}.gro object directly to build a tissue object in various
    configurations, which can later be used to add a matrix to if needed. 

    Methods:
    --------
    
    check_for_tissue():
        Checks whether 'tissues' are enabled in the input, otherwise passes class. 
    
    read_gro_cell():
        Reads .gro file of individual cell. Stores contents in 'self.' Reason for doing it like this instead of 
        accessing the data directly from CellConstructor is that this will make it easier to create non-uniform
        tissues later (by running CellConstructor multiple times) which then only requires minor edits to this class. 
    
    replicate_cell_on_grid():
        Creates basic standard packing by replicating on grid, stores information in 'self.tissue'. 

    replicate_cell_monolayer():
        Fits the requested number of cells on a monolayer in x- and y-direction (z = 1), stores information in 'self.tissue'. 
        
    replicate_cell_hexagonal():
        Creates a hexagonal geometry in the y-direction. The amount of shearing can be configured if needed. 
        Stores information in 'self.tissue'
        
    replicate_cell_disordered():
        Puts the requested number of cells in disordered fashion (like 'gmx insert-molecules') together based 
        on a contact distance logic. Stores information in 'self.tissue'. 
            Note: might lead to unnecessarily large box volumes for NVT.        
    """

    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        self.cell = CellTopology()
        self.cell_gro_content = {}
        self.tissue_gro_content = {}
      
          
            
    def read_gro_cell(self):
        """
        Reads .gro file of an individual cell, stores .gro file information in a dict 'self.cell_gro_content'
        """
        args = self.cli_parser.parse_args()
        
        with open(f'{args.output_dir}/{self.cell.groname}', 'r') as f:
            lines = f.readlines()

            nr_atoms = int(lines[1].strip())

            for line_index in range(2, 2 + nr_atoms):
                line = lines[line_index]

                atom_name = line[10:15].strip()
                x, y, z = map(float, [line[20:28], line[28:36], line[36:44]])
                
                self.cell_gro_content[line_index - 2] = {
                    "name": atom_name,
                    "coords": [x, y, z]
                }
                #'self.cell_gro_content' holds cell information as
                # key = atom index, 'name' = atomname, 'coords' = x,y,z
                
    def replicate_cell_on_grid(self):
        """
        Replicates cells uniformly on the smallest possible 3D (cubic) grid. 
        
        Todo 19/10: test by building gro from this directly --> and then later migrating to System class!!!!!
        """
        json_values = self.json_parser.json_values
        nr_of_cells = json_values["number_of_cells"]
        
        #offset between cells in nm
        offset = 4 # value can be played with, or can be exposed to JSON if needed. 
        
        total_offset = offset * (nr_of_cells - 1)

        # this expression tries to fit the requested number of cells on the smallest grid
        grid_size_x = int(np.ceil(nr_of_cells ** (1/3)))
        grid_size_y = int(np.ceil(nr_of_cells ** (1/3)))
        grid_size_z = int(np.ceil(nr_of_cells ** (1/3)))

        atom_index = 1  # Restart atom index at 1. 

        for i in range(grid_size_x):
            for j in range(grid_size_y):
                for k in range(grid_size_z):
                    for atom_index_cell, atom_info in self.cell_gro_content.items():
                        x, y, z = atom_info['coords']
                        new_coords = (
                            i * offset + x,
                            j * offset + y,
                            k * offset + z
                        )
                        self.tissue_gro_content[atom_index] = {
                            'name': atom_info['name'],
                            'coords': new_coords
                        }
                        atom_index += 1
    
class MatrixConstructor:
    """
    Depending on the input.JSON options, this class either uses the individual CELL-{timestamp}.gro object
    or tissue-{packing} object to create a matrix object. 

    Methods:
    --------
    
    process_matrix_settings():
        Checks the input for whether matrix is enabled, and checks the other matrix options.
        
    build_matrix_cell():
        Builds a matrix object for a single cell (if so requested). 
        
    build_matrix_tissue():
        Builds a matrix object for a tissue (if so requested). 
                
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        
class SystemConstructor:
    """
    Uses the tissue or matrix object from TissueConstructor / MatrixConstructor classes to make the final 'system-{timestamp}.gro'
    which will be input for the 'simulation' module. Then, based on the system specifications, creates a 'system-{timestamp}.top'
    file which links the necessary .itp files and prepares the input required for a GMX simulation. 
    
    fit_box_around_system():
        Fits a simulation box around the final system
        
    build_gro_file_system():
        Constructs the .gro file of the system
                
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        