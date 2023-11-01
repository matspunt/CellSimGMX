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
import logging
import datetime

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
        Checks whether 'tissues' are enabled in the input, otherwise passes class. Also figures out packing. 
    
    read_gro_cell():
        Reads .gro file of individual cell. Stores contents in 'self.' Reason for doing it like this instead of 
        accessing the data directly from CellConstructor is that this will make it easier to create non-uniform
        tissues later (by running CellConstructor multiple times) which then only requires minor edits to this class. 
    
    replicate_cell_on_grid():
        Creates basic standard packing by replicating on grid, stores information in 'self.tissue_coords'. 

    replicate_cell_monolayer():
        Fits the requested number of cells on a monolayer in x- and y-direction (z = 1), stores information in 'self.tissue_coords'. 
        
    replicate_cell_hexagonal():
        Creates a hexagonal geometry in the y-direction. The amount of shearing can be configured if needed. 
        Stores information in 'self.tissue_coords'
        
    replicate_cell_disordered():
        Puts the requested number of cells in disordered fashion (like 'gmx insert-molecules') together based 
        on a contact distance logic. Stores information in 'self.tissue_coords'. 
            Note: might lead to unnecessarily large box volumes for NVT.        
    """

    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        self.cell = CellTopology()
        self.cell_gro_content = {}
        self.tissue_coords = {}
        self.centered_tissue_coords = {}

        #Figure out whether we're dealing with a 'cell' or 'tissue' simulation, then continue based on that
        json_values = self.json_parser.json_values
        self.simulation_type = json_values["simulation_type"]
        
        #we need to check if this class is even needed for the requested simulation settings
        self.check_for_tissue()
        
    def check_for_tissue(self):
        """
        In this function we define what functions of this class need to be executed. If no tissue simulation was requested, do nothing.
        Otherwise pack the tissue in the requested manner. 
        """
        
        json_values = self.json_parser.json_values
        tissue_packing = json_values["tissue_packing"]
        
        if self.simulation_type == "cell":
            logging.warning(f"You requested '{self.simulation_type}' for simulation, ignoring 'tissue' logic")
            pass #do nothing in the tissue class, we already have the cell object and .gro!
        else:
            self.read_gro_cell()
            if tissue_packing == "grid":
                logging.info(f"You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                print(f"NOTE: You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                self.replicate_cell_on_grid()
                
            if tissue_packing == "hexagonal":
                logging.info(f"You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                print(f"NOTE: You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                self.replicate_cell_hexagonal()
                
            if tissue_packing == "monolayer":
                logging.info(f"You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                print(f"NOTE: You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                self.replicate_cell_monolayer()
                
            if tissue_packing == "disordered":
                logging.info(f"You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                print(f"NOTE: You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                self.replicate_cell_disordered()
        
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
                
    def fit_box_around_coord(self):
        """
        Fits closest fitting box around the coordinates. Note: will be moved to System class later!
        """
        
        min_x = min_y = min_z = float("inf")
        max_x = max_y = max_z = float("-inf")

        for atom_data in self.tissue_coords.values():
            x, y, z = atom_data["coords"]
            # Find minimum and maximum coordinates to figure out where box needs to go
            min_x = min(min_x, x)
            min_y = min(min_y, y)
            min_z = min(min_z, z)
            max_x = max(max_x, x)
            max_y = max(max_y, y)
            max_z = max(max_z, z)

        #Determine lengths of box vectors
        box_size_x = abs(max_x - min_x)
        box_size_y = abs(max_y - min_y)
        box_size_z = abs(max_z - min_z)
        
        return box_size_x, box_size_y, box_size_z
    
    def write_gro_tissue(self):
        """
        Note: this function will be moved to System class later, just for testing tissue packing functionality
        """
        
        args = self.cli_parser.parse_args()
               
        total_coords = len(self.tissue_coords)
        
        #Sum up all coordinates, then divide over the total number of atoms/coordinates
        sum_x = sum(float(atom_data["coords"][0]) for atom_data in self.tissue_coords.values())
        sum_y = sum(float(atom_data["coords"][1]) for atom_data in self.tissue_coords.values())
        sum_z = sum(float(atom_data["coords"][2]) for atom_data in self.tissue_coords.values())
        
        #Find geometric center in each Cartesian direction
        center_x = sum_x / total_coords
        center_y = sum_y / total_coords
        center_z = sum_z / total_coords
        
        #find the required box size
        box_size = self.fit_box_around_coord()
        
        # make the box vectors bigger to have some vacuum on edges to mitigate PBC effects in NVT
        edge_offset = 3 # in nm in .gro file, means 1.5 on each side
        
        extra_box_size = (
            box_size[0] + edge_offset,
            box_size[1] + edge_offset,
            box_size[2] + edge_offset
        )

        # Note: (0,0,0) is recognized as the box origin, thus need to translate the coordinates based
        # on the box dimensions. Calculate the translation factor of the positions to the geometric center
        trans_x = -(center_x - (extra_box_size[0] / 2))
        trans_y = -(center_y - (extra_box_size[1] / 2))
        trans_z = -(center_z - (extra_box_size[2] / 2))
        
        # Move the atom coordinates based on the calculated translation factor
        for atom_index, atom_data in self.tissue_coords.items():
            name = atom_data["name"]
            x, y, z = atom_data["coords"]

            # Apply the translation to the coordinates
            x_translated = x + trans_x
            y_translated = y + trans_y
            z_translated = z + trans_z

            self.centered_tissue_coords[atom_index] = {
                "name": name,
                "coords": [x_translated, y_translated, z_translated]
    }

        now = datetime.datetime.now()
        gro_header = "GRO file of TISSUE with {} cells at {}\n".format(str(len(self.tissue_coords)), now.strftime("%H:%M:%S"))
        
        self.groname = "TISSUE-{}.gro".format(now.strftime("%H-%M-%S"))
        
        with open(f"{args.output_dir}/{self.groname}", mode='w') as gro:
            gro.write(gro_header)
            gro.write(str(len(self.centered_tissue_coords)) + "\n")
            
            for atom_index, atom_data in self.centered_tissue_coords.items():
                atomname = atom_data["name"]
                coords = atom_data["coords"]
                x, y, z = coords
                resname = "CELL"
                atom_number = atom_index + 1
                
                # Format the coordinates in GRO format
                line = "{0:>5}{1:<5}{2:>5}{3:>5}{4:>8.3f}{5:>8.3f}{6:>8.3f}\n".format(
                    atom_number, resname, atomname, atom_number, x, y, z
                )

                gro.write(line)
            gro.write("{:>10.5}{:>10.5}{:>10.5f}\n".format(extra_box_size[0], extra_box_size[1], extra_box_size[2]))
        gro.close()
        
        logging.warning(f"Built a .gro file '{args.output_dir}/{self.groname}' of a tissue. ")
                 
    def replicate_cell_on_grid(self):
        """
        Replicates cells uniformly on the smallest possible 3D (cubic) grid. 
        """
        
        json_values = self.json_parser.json_values
        nr_of_cells = json_values["number_of_cells"]
        
        #offset between cells in nm
        offset = 4 # exact value can be played with, or can be exposed to JSON if needed. 
        
        # this expression tries to fit the requested number of cells on the smallest cubic grid
        grid_size_x = int(np.ceil(nr_of_cells ** (1/3)))
        grid_size_y = int(np.ceil(nr_of_cells ** (1/3)))
        grid_size_z = int(np.ceil(nr_of_cells ** (1/3)))

        atom_index = 0  # Restart atom index at 1. 

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
                        self.tissue_coords[atom_index] = {
                            'name': atom_info['name'],
                            'coords': new_coords
                        }
                        atom_index += 1
                        #self.tissue_coords is a dict formatted as: 11582: {'name': 'M2', 'coords': (11.697, 10.191, 11.756)}
                        
        self.write_gro_tissue()
    
    def replicate_cell_monolayer(self):
        """
        Replicates cells as a monolayer (cubic packing with a single layer in Z)
        """
        
        json_values = self.json_parser.json_values
        nr_of_cells = json_values["number_of_cells"]
        
        #offset between cells in nm
        offset = 4 # exact value can be played with, or can be exposed to JSON if needed. 

        # Determine the number of cells along each axis
        grid_size_x = int(np.ceil(np.sqrt(nr_of_cells)))
        grid_size_y = int(np.ceil(nr_of_cells / grid_size_x))
        grid_size_z = 1  # Limiting the grid size to 1 layer in the z-direction

        atom_index = 0  # Restart atom index at 1.

        # Generate the new coordinates based on the offset in x and y
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
                        self.tissue_coords[atom_index] = {
                            'name': atom_info['name'],
                            'coords': new_coords
                        }
                        atom_index += 1
        
        self.write_gro_tissue()
        
    def replicate_cell_hexagonal(self):
        print("You have enabled hexagonal!")
        ###Todo: expand on this class!
        
    def replicate_cell_disordered(self):
        print("You have enabled disordered!")
        
        
        
    
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
        