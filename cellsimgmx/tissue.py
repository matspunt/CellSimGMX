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

import numpy as np
import logging
import datetime
import sys
import re

from cellsimgmx import CLIParser
from cellsimgmx import JSONParser
from cellsimgmx import ForcefieldParserGMX
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
    """

    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        self.cell = CellTopology()
        self.cell_gro_content = {}
        self.tissue_coords = {}

        #Figure out whether we're dealing with a 'cell' or 'tissue' simulation, then continue based on that
        json_values = self.json_parser.json_values
        self.simulation_type = json_values["simulation_type"]
        
        #we need to check if this class is even needed for the requested simulation settings
        self.check_for_tissue()
        
    def check_for_tissue(self):
        """
        In this function we define what functions of this class need to be executed. If no tissue simulation was requested, do nothing.
        Otherwise pack the tissue in the requested manner. 
        
        Returns:
            self.<tissue_packing_function> - Dependent on user input
            pass - If tissue simulation is disabled. 
        """
        
        json_values = self.json_parser.json_values
        tissue_packing = json_values["tissue_packing"]
        
        if self.simulation_type == "cell":
            logging.warning(f"You requested '{self.simulation_type}' for simulation, ignoring 'tissue' related input")
            self.read_gro_cell() # the cell coordinates are all we need!
        else:
            self.read_gro_cell()
            if tissue_packing == "grid":
                logging.info(f"You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                print(f"INFO: You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                self.replicate_cell_on_grid()
                
            if tissue_packing == "hexagonal":
                logging.info(f"You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                print(f"INFO: You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                self.replicate_cell_hexagonal()
                
            if tissue_packing == "monolayer":
                logging.info(f"You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                print(f"INFO: You requested '{self.simulation_type}' for simulation, packing as '{tissue_packing}'")
                self.replicate_cell_monolayer()
        
    def read_gro_cell(self):
        """
        Reads .gro file of an individual cell, stores .gro file information. 
        
        Outputs:
            self.cell_gro_content (dict) - Cell coordinates where key = index, name = atomname, coord = x, y, z
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
                    for _, atom_info in self.cell_gro_content.items():
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
                    for _, atom_info in self.cell_gro_content.items():
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
        
    def replicate_cell_hexagonal(self):
        """
        Hexagonal packing based on a shearing factor and offset (both configurable). 
        """
        json_values = self.json_parser.json_values
        nr_of_cells = json_values["number_of_cells"]
        
        #offset between cells in nm
        offset = 4 # exact value can be played with, or can be exposed to JSON if needed. 
        hexagonal_shearing = 0.5 # the higher this is set, the more displaced layers are. 0.5 is hexagonal packing
        
        total_offset = offset * (nr_of_cells - 1) #required for uneven layers

        # Determine the number of cells along each axis
        grid_size_x = int(np.ceil(nr_of_cells ** (1/3)))
        grid_size_y = int(np.ceil(nr_of_cells ** (1/3)))
        grid_size_z = int(np.ceil(nr_of_cells ** (1/3)))
        
        atom_index = 0  # Restart atom index at 1.

        for i in range(grid_size_x):
            for j in range(grid_size_y):
                for k in range(grid_size_z):
                    for _, atom_info in self.cell_gro_content.items():
                        atom = atom_info["name"]
                        x, y, z = atom_info["coords"]
                        new_x_coord = x + i * offset - total_offset / 2
                        # Apply shearing to uneven layers only!!
                        new_y_coord = y + j * offset + (k % 2) * hexagonal_shearing * offset - total_offset / 2
                        new_z_coord = z + k * offset - total_offset / 2

                        self.tissue_coords[atom_index] = {
                            "name": atom, "coords": [new_x_coord, new_y_coord, new_z_coord]}
                        atom_index += 1
                             
class MatrixConstructor:
    """
    Depending on the input.JSON options, this class either uses the individual CELL-{timestamp}.gro object
    or tissue_coords object to create a matrix object. Combines individual object to build System coords.

    Methods:
    --------
    
    process_matrix_settings():
        Checks the input for whether matrix is enabled, and checks the other matrix options.
        
    build_matrix_cell():
        Builds a matrix object for a single cell (if so requested). 
        
    build_matrix_tissue():
        Builds a matrix object for a tissue (if so requested). 
        
    build_matrix_top():
        Generates a matrix topology (.itp) based on the requested matrix type. 
                
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        self.ff_parser = ForcefieldParserGMX()
        self.tissue = TissueConstructor() # tissue contains both the contents of CELL.gro and the tissue dicts + settings.
        json_values = self.json_parser.json_values
        self.matrix_on_off = json_values["matrix_on_off"]
        self.simulation_type = json_values["simulation_type"]
        self.tissue_packing = json_values["tissue_packing"]
        self.matrix_coords = {} # if enabled, contains JUST the matrix information
        self.system_coords = {} # will contain the final system information after class is executed!
        
        self.process_matrix_settings()
        
    def process_matrix_settings(self):
        """Checks whether the matrix is enabled. If not, do nothing and continue to System creation. If yes, 
        checks whether a cell or tissue simulation was requested, parses coordinates accordingly. If a 
        disordered packing was requested, lets user know this is not supported (yet). 
        
        Returns:
            self.build_matrix_cell (func) - In case cell simulation with matrix is requested
            self.build_matrix_tissue (func) - In case tissue simulation with matrix is requested. 
            self.system_coords (dict) - In case no matrix was requested. Contains cell or tissue coordinates
                depending on 'simulation' type JSON setting.
        """

        if self.matrix_on_off == "off" and self.simulation_type == "tissue":
            logging.warning(f"You have disabled the matrix ('{self.matrix_on_off}'), ignoring all matrix settings. ")
            #set system coordinates to that of the previously constructed tissue
            self.system_coords = self.tissue.tissue_coords.copy()
        if self.matrix_on_off == "off" and self.simulation_type == "cell":
            logging.warning(f"You have disabled the matrix ('{self.matrix_on_off}'), ignoring all matrix settings. ")
            #set system coordinates to that of an individual cell
            self.system_coords = self.tissue.cell_gro_content.copy()
        else:
            if self.matrix_on_off == "on" and self.simulation_type == "cell":
                logging.info(f"Currently fitting '{self.simulation_type}' object on a matrix...")
                print(f"INFO: Matrix is set to '{self.matrix_on_off}'. Building the matrix object. ")              
                self.build_matrix_cell()
            if self.matrix_on_off == "on" and self.simulation_type == "tissue":
                logging.info(f"Currently fitting '{self.simulation_type}' object on a matrix...")
                print(f"INFO: Matrix is set to '{self.matrix_on_off}'. Building the matrix object. ")             
                self.build_matrix_tissue()
                sys.exit(1)           
                
    def build_matrix_cell(self):
        """
        Generates coordinates and topology for a matrix object based on a single cell. 
        
        Parameters:
            self.tissue.cell_gro_content (dict): The individual cell coordinates
            
        Returns:
            self.matrix_coords (dict): Atom names and coords of matrix element (dict key = atom index). 
                This dict is required for constructing the topology. 
            self.system_coords (dict): Combined  dict of cell + matrix. 
                This dict is required for parsing the final system coordinates. 

        Todo:
             - Allow different matrix bead names than MX1
             - Multiple layers stack oddly
        """
        args = self.cli_parser.parse_args()

        json_values = self.json_parser.json_values
        matrix_layers = json_values["nr_of_layers"]
        #matrix_beads = json_values["matrix_beads"] - currently unused!
        
        matrix_z_offset = 0.4 # in nm, distance between cell and first matrix layer

        # we fit the matrix under the cell in z-direction, thus need to look for 'lowest'
        # z coordinate in the system (in this case, a single cell)
        lowest_z = min(atom['coords'][2] for atom in self.tissue.cell_gro_content.values())
        
        #want to center the matrix with respect to the cell on top
        center_x = sum(atom['coords'][0] for atom in self.tissue.cell_gro_content.values()) / len(self.tissue.cell_gro_content)
        center_y = sum(atom['coords'][1] for atom in self.tissue.cell_gro_content.values()) / len(self.tissue.cell_gro_content)

        #we can know the x,y length of the matrix plane too from the coordinates
        # and then multiply them by a matrix length factor afterwards, to have a larger matrix.
        max_x = max(atom['coords'][0] for atom in self.tissue.cell_gro_content.values())
        max_y = max(atom['coords'][1] for atom in self.tissue.cell_gro_content.values())
        matrix_x = max_x * 4
        matrix_y = max_y * 4
        
        atom_index = 0

        for layer in range(matrix_layers):
            z_offset = lowest_z - matrix_z_offset - (layer+1) * matrix_z_offset
            # Generate matrix particles centered at (0, 0, 0) with adjusted z coordinate
            for x in range(int(matrix_x) + 1):
                for y in range(int(matrix_y) + 1):
                    particle = {
                        'name': 'MX1',
                        'coords': [x - matrix_x / 2 + center_x, y - matrix_y / 2 + center_y, z_offset],
                    }
                    self.matrix_coords[atom_index] = particle
                    atom_index += 1
                
        #after preparing the matrix coordinates, we finalize the self.system_coords dict
        # continue from last index of CELL
        last_index = max(self.tissue.cell_gro_content.keys())
        self.system_coords = self.tissue.cell_gro_content.copy()

        next_index = last_index + 1
        # then add the matrix dict but index the atom indices again!
        for atom_index, atom_info in self.matrix_coords.items():
            new_atom_info = {
                'name': atom_info['name'],
                'coords': atom_info['coords']
            }
            self.system_coords[next_index] = new_atom_info
            next_index += 1
        
        logging.info(f"Succesfully constructed a matrix object of {len(self.matrix_coords)} particles")
        
        if args.verbose:
            print(f"\nSuccesfully constructed a matrix object of '{len(self.matrix_coords)}' particles\n")  
        
        self.build_matrix_top()
        
    def build_matrix_tissue(self):
        """
        Generates coordinates and topology for a matrix object based on a tissue.  
        
        Parameters:
            self.tissue.tissue_coords (dict): The individual cell coordinates
            
        Returns:
            self.matrix_coords (dict): Atom names and coords of matrix element (dict key = atom index)
                This dict is required for constructing the topology. 
            self.system_coords (dict): Combined  dict of tissue + matrix. 
                This dict is required for parsing the final system coordinates. 

        Todo:
             - Allow different matrix bead names than MX1
             - Multiple layers stack oddly
        """
        
        args = self.cli_parser.parse_args()
        
        json_values = self.json_parser.json_values
        matrix_layers = json_values["nr_of_layers"] #only 1 currently works well.
        #matrix_beads = json_values["matrix_beads"] - currently unused!
        
        matrix_z_offset = 0.4 # in nm, distance between cell and first matrix layer

        # logic is exactly the same as for single CELL, so removed comments
        lowest_z = min(atom['coords'][2] for atom in self.tissue.tissue_coords.values())
        center_x = sum(atom['coords'][0] for atom in self.tissue.tissue_coords.values()) / len(self.tissue.tissue_coords)
        center_y = sum(atom['coords'][1] for atom in self.tissue.tissue_coords.values()) / len(self.tissue.tissue_coords)

        max_x = max(atom['coords'][0] for atom in self.tissue.tissue_coords.values())
        max_y = max(atom['coords'][1] for atom in self.tissue.tissue_coords.values())
        matrix_x = max_x * 2
        matrix_y = max_y * 2

        atom_index = 0

        for layer in range(matrix_layers):
            z_offset = lowest_z - matrix_z_offset - (layer + 1) * matrix_z_offset
            # Generate matrix particles centered at (0, 0, 0) with adjusted z coordinate
            for x in range(int(matrix_x) + 1):
                for y in range(int(matrix_y) + 1):
                    particle = {
                        'name': 'MX1',
                        'coords': [x - matrix_x / 2 + center_x, y - matrix_y / 2 + center_y, z_offset],
                    }
                    self.matrix_coords[atom_index] = particle
                    atom_index += 1
        
        #after preparing the matrix coordinates, we finalize the self.system_coords dict
        # continue from last index of CELL
        last_index = max(self.tissue.tissue_coords.keys())
        self.system_coords = self.tissue.tissue_coords.copy()
        
        next_index = last_index + 1
        # then add the matrix dict but index the atom indices again!
        for atom_index, atom_info in self.matrix_coords.items():
            new_atom_info = {
                'name': atom_info['name'],
                'coords': atom_info['coords']
            }
            self.system_coords[next_index] = new_atom_info
            next_index += 1
            
        logging.info(f"\nSuccesfully constructed a matrix object of '{len(self.matrix_coords)}' particles\n")
        
        if args.verbose:
            print(f"\nSuccesfully constructed a matrix object of '{len(self.matrix_coords)}' particles\n")  
        
        self.build_matrix_top()
        
    def build_matrix_top(self):
        """
        Constructs the topology of a matrix based on the 'self.matrix_coords' object. Applies position restraints
        on each atom in the matrix. Agnostic to bead names (as long as they occur in the force field!). As much
        information as possible is taken from the forcefield.itp given by the user. 
        
        Works with both matrix for a single cell, and tissue. 
        
        Outputs:
            A topology MX-{timestamp}.itp 
        """
        
        args = self.cli_parser.parse_args()
        now = datetime.datetime.now()
        
        self.itpname = "MX-{}.itp".format(now.strftime("%H-%M-%S"))
        
        with open(f"{args.output_dir}/toppar/{self.itpname}", "w") as itp:
            #Write the topology header
            header = "; Topology file for matrix built at {}\n".format(now.strftime("%H:%M:%S"))
            itp.write(header)
            ff_itp =  f"; Copied the force field into 'toppar' from:\n;#include \"{self.ff_parser.itp_path}\"\n"
            itp.write(ff_itp)
            
            # Write the [moleculetype] directive and the [atoms] directive based on the matrix_coords dict
            itp.write("\n[ moleculetype ]\n; Name        nrexcl\n  MX          1\n\n[ atoms ]\n; nr type resnr residue atom cgnr  charge   mass\n")

            for index, atom_info in self.matrix_coords.items():  
                resname = "MX"
                atom_nr = index + 1
                atomname = atom_info['name']
                mass = self.ff_parser.atomtypes[atomname]['mass']
                itp.write("  {:<3s}  {:<3s}   1    {:<3s}    {:<3s}  {:<3s} 0.0000  {:<3s}\n".format(str(atom_nr), atomname, resname, atomname, str(atom_nr), str(mass)))

            # matrix particles are not bound to each other - they solely interact through LJ
            # but we do apply position restraints on each matrix atom 
            # the value of position restraint is tuned on the .mdp level
            itp.write("\n#ifdef POSRES\n[ position_restraints ]\n; i  funct          fcx              fcy                 fcz\n")
            for index, _ in self.matrix_coords.items():  
                atom_nr = index + 1
                posres = "POSRES_MATRIX_FC"
                itp.write("  {:<3s}   1    {:<6s}    {:<6s}      {:<6s}\n".format(str(atom_nr), posres, posres, posres))
            
        logging.info(f"A matrix topology file '{args.output_dir}/{self.itpname}' has been built")
        itp.close()
        
        if args.verbose:
            print(f"\nVERBOSE MODE ENABLED. A matrix topology file '{args.output_dir}/{self.itpname}' has been built")
        
class SystemConstructor:
    """
    Uses the tissue or matrix object from TissueConstructor / MatrixConstructor classes to make the final 'system-{timestamp}.gro'
    which will be input for the 'simulation' module. Then, based on the system specifications, creates a 'system-{timestamp}.top'
    file which links the necessary .itp files and prepares the input required for a GMX simulation. 
    
    fit_box_around_system():
        Fits a simulation box around the final system. 
        
    write_gro_system():
        Constructs the .gro file of the system
        
    construct_system_topology():
        Creates a file 'system.top' that contains the necessary information about the system and links the correct
        .itps
                
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        self.centered_system_coords = {}
        self.matrix = MatrixConstructor()
        self.write_gro_system()
        #self.construct_system_topology()
        
    def fit_box_around_coord(self):
        """
        Fits closest fitting (cubic) box around the given coordinates.
        
        Returns:
            box_size_x, box_size_y, box_size_y (floats) - Vectors of the box.
        """
        
        min_x = min_y = min_z = float("inf")
        max_x = max_y = max_z = float("-inf")

        for atom_data in self.matrix.system_coords.values():
            x, y, z = atom_data["coords"]
            # Find minimum and maximum coordinates to figure out where box needs to go
            min_x = min(min_x, x)
            min_y = min(min_y, y)
            min_z = min(min_z, z)
            max_x = max(max_x, x)
            max_y = max(max_y, y)
            max_z = max(max_z, z)

        box_size_x = abs(max_x - min_x) #length of box vectors
        box_size_y = abs(max_y - min_y)
        box_size_z = abs(max_z - min_z)
        
        return box_size_x, box_size_y, box_size_z
    
    def write_gro_system(self):
        """
        This function wrtites the coordinates (.GRO) of particles of a System. A 'System' is defined
        as per the input.JSON. It can contain a cell, a tissue, a tissue with matrix etc dependent
        on input logic. GRO writing is independent of input settings however. 
        
        By default, atoms are centered in the box. If desired, a customizable offset can be added
        to increase the vacuum in the box to counteract PBC effects, or make minimization work. 
        
        Outputs:
            SYSTEM-{timestamp}.gro - Containing the coordinates of the final system. 
        """
        
        json_values = self.json_parser.json_values
        box_coord_offset = json_values["box_coord_offset"] #used for box fitting
        nr_of_particles = json_values["nr_of_particles"] + 1 #including center-particle
        resid = 0
        # in case set to 3.0 for instance, this means 1.5 on each side
        # in case set to 0, the tightest possible box is simply fit. 
        
        centering = True # use to enable or disable centering in the box (useful for debugging)
        
        args = self.cli_parser.parse_args()
        
        total_coords = len(self.matrix.system_coords)
        
        #Sums up all coordinates, then divides over the total number of atoms/coordinates...
        sum_x = sum(float(atom_data["coords"][0]) for atom_data in self.matrix.system_coords.values())
        sum_y = sum(float(atom_data["coords"][1]) for atom_data in self.matrix.system_coords.values())
        sum_z = sum(float(atom_data["coords"][2]) for atom_data in self.matrix.system_coords.values())
        
        #... to find geometric center in each Cartesian direction
        center_x = sum_x / total_coords
        center_y = sum_y / total_coords
        center_z = sum_z / total_coords
        
        #calculate the required box size (add offset later, in case of 0, nothing happens)
        box_size = self.fit_box_around_coord()
        
        extra_box_size = (
            box_size[0] + box_coord_offset,
            box_size[1] + box_coord_offset,
            box_size[2] + box_coord_offset
        )

        # Note: (0,0,0) is recognized as the box origin in GMX, thus need to translate the coordinates based
        # on the box dimensions. Calculate the translation factor of the positions to the geometric center      
        trans_x = -(center_x - (extra_box_size[0] / 2))
        trans_y = -(center_y - (extra_box_size[1] / 2))
        trans_z = -(center_z - (extra_box_size[2] / 2))
        
        # Move the atom coordinates based on the calculated translation factor
        for atom_index, atom_data in self.matrix.system_coords.items():
            name = atom_data["name"]
            x, y, z = atom_data["coords"]

            x_translated = x + trans_x
            y_translated = y + trans_y
            z_translated = z + trans_z
                
            self.centered_system_coords[atom_index] = {
                "name": name,
                "coords": [x_translated, y_translated, z_translated]
            }

        now = datetime.datetime.now()
        gro_header = "GRO file of SYSTEM with {} particles written at {}\n".format(str(len(self.matrix.system_coords)), now.strftime("%H:%M:%S"))
        
        self.groname = "SYSTEM-{}.gro".format(now.strftime("%H-%M-%S"))
        
        with open(f"{args.output_dir}/{self.groname}", mode='w') as gro:
            gro.write(gro_header)
            gro.write(str(len(self.centered_system_coords)) + "\n")
            
            if centering == True:
                logging.info("Centering of System coordinates is enabled!")
                #take coordinates from translated (centered) dict
                for atom_index, atom_data in self.centered_system_coords.items():
                    atomname = atom_data["name"]
                    coords = atom_data["coords"]
                    x, y, z = coords
                    #look for atomnames beginning with MX, this is Matrix res
                    if re.match(r'^MX[1-5]$', atomname):
                        resname = "MX  "
                    else:
                        # rest will be CELL
                        resname = "CELL"
                    
                    # resid increases per every increase in the number of particles per cell
                    if atom_index % nr_of_particles == 0:
                        resid += 1

                    atom_number = atom_index + 1
                    
                    # Format the coordinates in GRO format
                    line = "{0:>5}{1:<5}{2:>5}{3:>5}{4:>8.3f}{5:>8.3f}{6:>8.3f}\n".format(
                        resid, resname, atomname, atom_number, x, y, z
                    )

                    gro.write(line)
                    
            if centering == False:
                logging.info("Centering of System coordinates is disabled!")
                #take coordinates from non-translated dict
                for atom_index, atom_data in self.matrix.system_coords.items():
                    atomname = atom_data["name"]
                    coords = atom_data["coords"]
                    x, y, z = coords
                    #look for atomnames beginning with MX, this is Matrix res
                    if re.match(r'^MX[1-5]$', atomname):
                        resname = "MX  "
                    else:
                        # rest will be CELL
                        resname = "CELL"
                    atom_number = atom_index + 1
                    
                    # Format the coordinates in GRO format
                    line = "{0:>5}{1:<5}{2:>5}{3:>5}{4:>8.3f}{5:>8.3f}{6:>8.3f}\n".format(
                        atom_number, resname, atomname, atom_number, x, y, z
                    )

                    gro.write(line)
            gro.write("{:>10.5}{:>10.5}{:>10.5f}\n".format(extra_box_size[0], extra_box_size[1], extra_box_size[2]))
        gro.close()
        
        logging.info(f"Built a .gro file '{args.output_dir}/{self.groname}' of the final requested system ")
        print(f"INFO: Coordinate generation finished, '{args.output_dir}/{self.groname}' is saved. ")
        
        self.construct_system_topology()
        logging.info(f"Topology and subtopologies '{args.output_dir}/system.top' built of the final system ")
        print(f"INFO: All input creation is completed. ")
        
    def construct_system_topology(self):
        """
        Constructs the top file required to run the simulation, links itps and forcefield
        """
        
        args = self.cli_parser.parse_args()
               
        json_values = self.json_parser.json_values
        matrix_on_off = json_values["matrix_on_off"]
        simulation_type = json_values["simulation_type"]
        
        with open(f"{args.output_dir}/system.top", "w") as top:
            forcefield =  f"#include \"{args.output_dir}/toppar/forcefield.itp\"\n"
            top.write(forcefield)
            cell_itp = f"#include \"{args.output_dir}/toppar/{self.matrix.tissue.cell.itpname}\"\n"
            top.write(cell_itp)
            if matrix_on_off == 'on':
                matrix_itp = f"#include \"{args.output_dir}/toppar/{self.matrix.itpname}\"\n"
                top.write(matrix_itp)
            if simulation_type == "cell":
                number_of_cells = 1
            if simulation_type == "tissue":
                number_of_cells = json_values["number_of_cells"]
            top_info =  f"\n[ system ]\nCellSimGMX system\n\n[ molecules ]\nCELL {number_of_cells}" 
            top.write(top_info)
            #the matrix individual particles are all written to .itp because of posres so there is only
            # a single matrix 'molecule' in the simulation
            matrix = "\nMX      1"
            if matrix_on_off == 'on':
                top.write(matrix)
        top.close()