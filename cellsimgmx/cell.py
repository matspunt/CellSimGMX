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
import matplotlib.pyplot as plt
from scipy.spatial import KDTree #for n-nearest neighbour search
import datetime
import random
import logging
import os
import sys

from cellsimgmx import CLIParser
from cellsimgmx import JSONParser
from cellsimgmx import ForcefieldParserGMX

class CellConstructor:
    """
    This class builds coordinates of a single Cell, the basic building block of a System based on the input.JSON. We use the number of particles
    and cell radius to build an initial .gro file into the desired shape, and then in CellTopology we assign the atom names, bonds etc. to make the topology. 
    
    Methods:
    pack_particles_on_shape():
        Pack particles onto a specified shape within a cell of given radius.
        
    Returns:
    self.particles():
        A list of the particle positions where (0,0,0) is always the center particle. 
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        self.particles = [(0, 0, 0)] # initialize with center particle, we store particle positions in this list
             
    def pack_particles_on_shape(self):
        """
            Takes .JSON input and uses the information there (nr_of_particles, cell_radius, and packing_shape) to create coordinates of a cell. 
            The packing for ellipsoids is currently bad. Note that there is no direct error handling here, since the user is not expected to change
            the number of particles or radius themselves --> these will only be changed in the development phase. Since the parameters can be
            printed in -verbose mode that should make finding problems clear enough. Returns the particle positions. 
        """
    
        args = self.cli_parser.parse_args()
        
        json_values = self.json_parser.json_values # don't need  to check for validity here, if problems in JSON, programme wil already have terminated
        
        nr_of_particles = json_values["nr_of_particles"]
        cell_radius = json_values["cell_radius"]
        shape = json_values["initial_packing_shape"]
        
        logging.info(f"You have specified the shape '{shape}' for a single cell")

        # Pack the particles on sphere using the Golden Spiral method 
        # see also: https://extremelearning.com.au/evenly-distributing-particles-on-a-sphere/
        # there are better methods, but good enough for GMX minimization
        if shape == 'spherical':
            theta = np.pi * (3.0 - np.sqrt(5.0))  # the Golden angle in radians
            for i in range(nr_of_particles):
                y = 1 - (i / float(nr_of_particles - 1)) * 2  # such that y-coordinates are between -1,1
                radius = np.sqrt(1 - y * y) * cell_radius # then scale coords by requested radius (in nm)
                phi = i * theta
                x = np.cos(phi) * radius
                z = np.sin(phi) * radius
                self.particles.append((x, y * cell_radius, z))
                                
        elif shape == 'cuboid':
            # Generate a face centered cubic (FCC) lattice
            # define particles per lattice to fit them on
            nr_of_particles = nr_of_particles + 1 # to account for center particle
            n = int(np.ceil(nr_of_particles**(1/3)))
            for i in range(n):
                for j in range(n):
                    for k in range(n):
                        if len(self.particles) < nr_of_particles:
                            x = (i * 2 - n + 1) * (cell_radius*0.2) #need to scale by a factor of ~1/5th to obtain similar volume as sphere
                            y = (j * 2 - n + 1) * (cell_radius*0.2)
                            z = (k * 2 - n + 1) * (cell_radius*0.2)
                            self.particles.append((x, y, z))
                           
        elif shape == 'ellipsoid':
            # Ellipsoids are much more complicated and formal methods require solving of the elliptic integral: https://mathworld.wolfram.com/EllipticIntegraloftheSecondKind.html
            # A formal Python implementation can be found here: https://github.com/maxkapur/param_tools
            # Acceptance/rejection criteria scale poorly but since the number of particles we need is comparatively small this is good enough 
            # (and I wouldn't even know how to implement a geometrical method ^^)
            
            while len(self.particles) < nr_of_particles + 1 : #account for center particle here!
                # First generate random spherical coordinates
                # and scale by the coordinates by the poles of the ellipsoid
                theta = np.random.uniform(0, 2 * np.pi)
                phi = np.arccos(2 * np.random.uniform(0, 1) - 1)
                
                ellips_scaling = [1.5,2.5, 4.0] #x, y, z
                
                x = cell_radius * ellips_scaling[0] * np.sin(phi) * np.cos(theta)
                y = cell_radius * ellips_scaling[1] * np.sin(phi) * np.sin(theta)
                z = cell_radius * ellips_scaling[2] * np.cos(phi)
                
                new_particle = (x, y, z)
                 
                # Check whether the particle is not too close to existing particles
                is_not_too_close = True
                for particle in self.particles:
                    #calculate Euclidean distance between each particle
                    distance = np.sqrt(sum((np.array(particle) - np.array(new_particle))**2))
                    if distance < 1.0:
                        is_not_too_close = False
                        break
                
                if is_not_too_close:
                    self.particles.append(new_particle)
                    
        logging.info(f"An object of shape '{shape}' has been constructed.\n")
    
        if args.output_dir:
            if not os.path.exists(args.output_dir):
                logging.error(f"Error: Output dir '{args.output_dir}' does not exist. Please check for typos and try again")
                sys.exit(1)
                
            if not os.access(args.input_dir, os.R_OK):
                logging.error(f"Error: Output directory '{args.output_dir}' is not readable. Do you have permissions and is the directory correct?")
                sys.exit(1)
                
        else:
            logging.error("Error: Output dir for settings not specified. Please specify using '--output-dir' and try again")
            sys.exit(1)
                        
        if args.verbose:
            print(f"\nVERBOSE MODE ENABLED. Currently creating a single cell object with {nr_of_particles} particles...\n") 
        
            ### if --verbose flag enabled, save an image plot of the packed cell as a quick reference
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            x = [particle[0] for particle in self.particles]
            y = [particle[1] for particle in self.particles]
            z = [particle[2] for particle in self.particles]

            colors = ['r' if (xi, yi, zi) == (0, 0, 0) else 'b' for xi, yi, zi in zip(x, y, z)]

            ax.scatter(x, y, z, c=colors, marker='o')
            ax.set_title('Packing of beads on cell surface')

            now = datetime.datetime.now()
            figname = "CELL_packing-{}.png".format(now.strftime("%H-%M-%S"))
            plt.savefig(f"{args.output_dir}/{figname}")
            print(f"{figname} has been created. Please ensure the packing is visually correct, should you encounter any problems at a later stage. ")
            
        return self.particles
    
class CellTopology:
    """
    Uses particle positions from CellConstructor and assigns them particle (atom) names in different logic. Then saves a .gro and .itp based on the force field.
    These form the basic building blocks of a tissue 

    Methods:
    assign_atom_names():
        Takes list of particle positions and assigns them atom names. 
        
    find_nearest_neighbours():
        Looks for nearest neighbours based on the input.json and stores that information as the self.atomname indices to later parse the bonds accordingly. 
        
    build_gro_file_cell():
        Saves a .gro file with the right atom names. Can be simulated individually if needed. Otherwise, can be used to build the tissue. 
        
    build_top_from_cell():
        Creates an .itp file for a single cell based on the atomname and force field information, recognizes center bead. If nearest neighbours are set, also parses
        these bonds in the .itp
    
    """
    
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        self.ff_parser = ForcefieldParserGMX()
        self.cell = CellConstructor()
        self.particles = None
        self.atomnames = {(0,0,0): "C"} #define a dictionary where the keys will be the atom names and data are the particle positions
        self.nnneighbours = {} # dict that stores the indices and atomnames of n-nearest neighbours in the membrane
        self.gro = None
        self.itp = None
        
    def assign_atom_names(self):
        args = self.cli_parser.parse_args()
        # access particle positions
        self.particles = self.cell.pack_particles_on_shape()
        
        #first need to check for selection of membrane and junction beads
        json_values = self.json_parser.json_values
        
        membrane_type = json_values["membrane_type"]
        membrane_beads = json_values["membrane_beads"]
        membrane_beads = membrane_beads.split(", ")
        
        #First we make all positions membrane type based on input preferences!
        if membrane_type == "even":
            particles = [particle for particle in self.particles if particle != (0, 0, 0)] #exclude center particle

            #count how often each membrane bead should be used, assuming even distribution
            count_each_memb_type = len(particles) // len(membrane_beads)
            atomnames_list = membrane_beads * count_each_memb_type
            random.shuffle(atomnames_list)
            
            #then assign a random atom name to all membrane particles as a dict, where the key is the coordinate
            self.atomnames.update({particles[i]: atomnames_list[i] for i in range(len(particles))})
            
        if membrane_type == "segmented":
            particles = [particle for particle in self.particles if particle != (0, 0, 0)] #exclude center particle

            count_each_memb_type = len(particles) // len(membrane_beads)
            atomnames_list = membrane_beads * count_each_memb_type
            atomnames_list = sorted(atomnames_list, key=lambda atomname: membrane_beads.index(atomname)) #sort the atomnames in the order of the input (such that 'M1' comes first, then 'M2' etc.)
            self.atomnames.update({particles[i]: atomnames_list[i] for i in range(len(particles))}) #update atomnames, this will form segments of the cell. 
            
        # Then, if junction beads are requested in the input, we can exchange some membrane beads for them. 
        junction_beads = json_values["junction_beads"]
        junction_bead_ratio = json_values["junction_bead_ratio"]
        
        if junction_beads != "off":
            junction_beads = junction_beads.split(", ")
            
            particles = len(self.atomnames) - 1 
            
            #calculate nr of junction beads required and split that by the different types of junction beads specified
            nr_of_junction_beads = junction_bead_ratio * particles
            nr_of_junction_beads = int(nr_of_junction_beads) #cast to integer

            count_per_junction = nr_of_junction_beads // len(junction_beads)
            
            junctionnames_list = junction_beads * count_per_junction
       
            # define a random selection of our 'self.atomnames' dict with membrane beads that we are going to replace
            # note: we are not editing the keys (particle coordinates), only the atom name!
            replace_atomnames = random.sample(list(self.atomnames.keys()), len(junctionnames_list))
            
            for coord, new_atomname in zip (replace_atomnames, junctionnames_list):
                self.atomnames[coord] = new_atomname 
            
        else:
            #if no junction beads were found, continue as normal and make no edits to 'self.atomnames'
            pass
        
        logging.info(f"The atom names have been assigned, the resulting atom names in a single cell are:")
        just_atomnames = list(map(lambda x: x[1], self.atomnames.items()))
        logging.info(f"{just_atomnames}")
        
        if args.verbose:
            print(f"\nVERBOSE MODE ENABLED. Final atom names have been assigned, the resulting atom names in a single cell are:\n")
            print(*(map(lambda x: x[1], self.atomnames.items())))
            print(f"\nMembrane beads have been packed in a '{membrane_type}' fashion, if you enabled junction beads, they are randomly distributed over the surface.")
            
            
    def find_nearest_neighbours(self):
        """
            Uses output from assign_atom_names() to determine the n-nearest neighbours on the surface membrane, if enabled by the user in 'input.JSON'. 
            
            Output:
                A dict (self.nnneighbours) with as keys the indices of the connected neighbours and as value their atomnames. 
        """
        args = self.cli_parser.parse_args()
        
        json_values = self.json_parser.json_values
        nearest_neighbour_springs = json_values["nearest_neighbour_springs"]
        
        if nearest_neighbour_springs != "off":
            coords, atomnames = zip(*self.atomnames.items()) #store information separately for indexing later in the loop

            tree = KDTree(coords)

            for index, (coord, atomname) in enumerate(self.atomnames.items()):
                _, neighbour_indices = tree.query(coord, k=3)  # Adjust k as needed
                key = f"{index + 1} {atomname}"
                self.nnneighbours[key] = [(i+1, atomnames[i]) for i in neighbour_indices if i != index]
                
            print(self.nnneighbours)
            
            """"
            SKIP CENTER BEAD IN NEIGHBOUR LIST!!!!!!
            """
            
        else:
            #if nearest neighbour springs are disabled then do nothing in this function
            pass 
 
    def build_gro_file_cell(self):
        """
        Uses output from assign_atom_names() to save the atomnames and coordinates in a .gro file
        
        Output:
            A GROMACS compatible coordinate (.gro) file on disk formatted as CELL-{timeprint}.gro

        """
        args = self.cli_parser.parse_args()
        
        now = datetime.datetime.now()
        gro_header = "GRO file of CELL with {} particles at {}\n".format(str(len(self.particles)), now.strftime("%H:%M:%S"))
        
        groname = "CELL-{}.gro".format(now.strftime("%H-%M-%S"))

        with open(f"{args.output_dir}/{groname}", mode='w') as gro:
            gro.write(gro_header)
            gro.write(str(len(self.particles)) + "\n")

            for i, (coord, atomname) in enumerate(self.atomnames.items()):  
                atomname = atomname
                x, y, z = coord
                resname = "CELL"
                atom_number = i + 1
                
                # Format the coordinates in GRO format
                line = "{0:>5}{1:<5}{2:>5}{3:>5}{4:>8.3f}{5:>8.3f}{6:>8.3f}\n".format(
                    atom_number, resname, atomname, atom_number, x, y, z
                )

                gro.write(line)
            gro.write("{:>10.5f}{:>10.5f}{:>10.5f}\n".format(5.000, 5.000, 5.000))  # Add box (arbitrary size)
        
        logging.info(f"Built a .gro file {groname} of a single cell. ")
        
        if args.verbose:
            print(f"\nVERBOSE MODE ENABLED. A coordinate file {groname} has been saved. Please inspect it in case problems arise later. \n")

    def build_top_from_cell(self):
        """
            Uses output from assign_atom_names() and find_nearest_neigbors() to build a topology based on the specified force field. 
            Basically links self.atomnames dict information (the beads in the cell) to ForceFieldParserGMX() class information. 
            
            Outputs:
                CELL-{timeprint}.itp file on disk
        """
        args = self.cli_parser.parse_args()
        
        now = datetime.datetime.now()
        itpname = "CELL-{}.itp".format(now.strftime("%H-%M-%S"))
        
        with open(f"{args.output_dir}/{itpname}", "w") as itp:
            #Write the topology header
            header = "; Topology file for a single CELL generated at {}\n".format(now.strftime("%H:%M:%S"))
            itp.write(header)
            ff_itp =  f"; Using forcefield from:\n#include \"{self.ff_parser.itp_path}\"\n"
            itp.write(ff_itp)
            
            # Write the [moleculetype] directive and the [atoms] directive based on the atomnames dict
            itp.write("\n[ moleculetype ]\n; Name        nrexcl\n  CELL        1\n\n[ atoms ]\n; nr type resnr residue atom cgnr  charge   mass\n")
            for i, (coord, atomname) in enumerate(self.atomnames.items()):  
                resname = "CELL"
                atom_nr = i + 1 
                mass = self.ff_parser.atomtypes[atomname]['mass']
                itp.write("  {:<3s}  {:<3s}   1    {:<3s}    {:<3s}  {:<3s} 0.0000  {:<3s}\n".format(str(atom_nr), atomname, resname, atomname, str(atom_nr), str(mass)))

            # writing [bonds] directive
            # First, connect all membrane beads to the center based on their atomnames. That is, each atom is connected to C bead (first atom)
            # do this by comparing the self.atomnames dict to the bondtype dict read from the force field and look for bonds with a C bead in them
            itp.write("\n[ bonds ]\n; i j func   r0   fk\n; center - membrane bonds\n")
            
            for i, (coord, atomname) in enumerate(self.atomnames.items()):
                if atomname != 'C': #skip the center bead
                    parsed = False #because of the two if statements required in both dicts, the for loop prints multiple times. Workaround (I am a genius coder ;D)
                    for bond, entry in self.ff_parser.bondtypes.items():
                        if "C" in bond and atomname in bond:
                            if not parsed:
                                atom_index = i + 1
                                #by construction, the Center bead is indexed as 1
                                itp.write(" 1  {:<3s} {:<3s} {:<3s} {:<3s} \n".format(str(atom_index), str(entry['func']), str(entry['r0']), str(entry['fk'])))
                                parsed = True
        
            #Finally, we take the self.nnneighbour dict
                        
        itp.close()
