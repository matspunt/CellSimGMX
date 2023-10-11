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

"""
    Todo: 
    
    1) make cuboid and ellipsoid packing functional!
    2) Add topology parser from old code
    3) 
    
"""

import numpy as np
import matplotlib.pyplot as plt
import datetime
import logging
import tempfile

from cellsimgmx import CLIParser
from cellsimgmx import JSONParser
from cellsimgmx import ForcefieldParserGMX

class CellConstructor:
    """
    This class builds a .gro file of a single Cell, the basic building block of a System based on the input.JSON. We we use the number of particles
    and cell radius to build an initial .gro file into the desired shape, and then in CellTopology we assign the atom names, bonds etc. to make the topology. 
    
    Methods:
    pack_particles_on_shape():
        Pack particles onto a specified shape within a cell of given radius.

    build_gro_from_cell_object():
        Generate a Gromacs compatible GRO file from a cell object defined in previous function. 
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        self.particles = [(0, 0, 0)] # initialize with center particle, we store particle positions in this list
        self.tmpgro = None #store a temporary gro file of the cell that can be used to apply 
             
    def pack_particles_on_shape(self):
        """
            Takes .JSON input and uses the information there (nr_of_particles, cell_radius, and packing_shape) to create coordinates of a cell. 
            The packing for ellipsoids is currently bad. Note that there is no direct error handling here, since the user is not expected to change
            the number of particles or radius themselves --> these will only be changed in the development phase. Since the parameters can be
            printed in -verbose mode that should make finding problems easy enough imo. 
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
            """
            The following works but the packing is really crappy! Come up with a better way than this in future. 
            """
            while len(self.particles) < nr_of_particles + 1 : #account for center particle here!
                # First generate random spherical coordinates
                theta = np.random.uniform(0, 2 * np.pi)
                phi = np.arccos(2 * np.random.uniform(0, 1) - 1)
                
                # Then scale the coordinates unevenly to obtain some ellipsoid shape
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
                    if distance < 0.3: #this parameter needs to be tuned!
                        is_not_too_close = False
                        break
                
                #accept if outside of threshold
                if is_not_too_close:
                    self.particles.append(new_particle)
    
        logging.info(f"A single cell of shape '{shape}' has been constructed.")

        if args.verbose:
            print(f"\nVERBOSE MODE ENABLED. Currently creating a single cell object of {nr_of_particles}...\n") 
        
            ### if --verbose flag enabled, save an image plot of the packed cell as a quick reference
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            x = [point[0] for point in self.particles]
            y = [point[1] for point in self.particles]
            z = [point[2] for point in self.particles]

            colors = ['r' if (xi, yi, zi) == (0, 0, 0) else 'b' for xi, yi, zi in zip(x, y, z)]

            ax.scatter(x, y, z, c=colors, marker='o')
            ax.set_title('Packing of beads on cell surface')

            now = datetime.datetime.now()
            figname = "CELL_packing-{}.png".format(now.strftime("%H-%M-%S"))
            plt.savefig(figname)
            print(f"{figname} has been created. Please ensure the packing is visually correct, should you encounter any problems at a later stage. ")
            
    def build_tmp_gro_from_object(self):
        total_particles = len(self.particles)
        gro_header = "GRO file of CELL with {} particles at {}\n".format(str(total_particles), datetime.datetime.now().strftime("%H:%M:%S"))

        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.gro') as gro:
            # this file is accessible through gro.name
            print(gro.name)
            gro.write(gro_header)
            gro.write(str(total_particles) + "\n")

            # Loop through the list of particle coordinates
            for i, coord in enumerate(self.particles):
                atom_name = 'A'  # Doesn't matter, is going to be replaced later anyway
                resname = 'CELL'
                atom_number = i + 1
                x, y, z = coord

                # Format the coordinates in GRO format
                line = "{0:>5}{1:<5}{2:>5}{3:>5}{4:>8.3f}{5:>8.3f}{6:>8.3f}\n".format(
                    atom_number, resname, atom_name, atom_number, x, y, z
                )

                gro.write(line)
            gro.write("{:>10.5f}{:>10.5f}{:>10.5f}\n".format(5.000, 5.000, 5.000))  # Add arbitrary box

    #class CellTopology:

        #In 'CellConstructor' we have built a cell object and .gro file with non specific atom types. Now we need to process the .gro file
