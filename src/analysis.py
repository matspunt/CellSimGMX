import os
import math
import numpy as np

class SingleCellAnalysis:
    """
    Analysis module that reads a .gro file containing a single cell and calculates various properties of the surface particles. 
    """
    
    def helper_gro_parser(self, gro_path):
        """
        Parses a .gro file and returns a list of tuples containing the parsed information.
        
        gro_path (str): The path to the .gro file.
        """
        abspath = os.path.abspath(gro_path)
        gro_list = []

        with open(gro_path, "r") as gro:
            gro.readline()
            gro_lines = gro.readlines()[1:-1]
            for line in gro_lines:
                line_split = line.split()
                atom_name = line_split[1]
                atom_index = int(line_split[2])
                x = float(line_split[3])
                y = float(line_split[4])
                z = float(line_split[5])

                gro_list.append((atom_index, atom_name, x, y, z))
                
        return gro_list
    
    def calculate_particle_density(self, gro_path, subset_area, bin_size, sigma, bond_radius):
        """
        Calculates the averaged density of particles within a specified subset of the surface of a single CELL.

        Parameters:
        - gro_path (str): The path to the .gro file.
        - subset_area (float): The area of the subset region.
        - sigma (float): The sigma value used to calculate the particle radius.

        Returns:
        - averaged_density (float): The averaged density of particles within the subset area.
        """
        
        ### Should we dynamically calculate radius by calculating distance from Center bead?
        ### Fuck this to be honest. 
        
        atom_radius = (1.12246 * sigma) / 2 # atom radius is based on minimum VdW distance

        gro_list = self.helper_gro_parser(gro_path)
        gro_list.pop() # removes N-bead
        
        total_atoms = len(gro_list)
        
        total_area = 4 * math.pi * (bond_radius**2) #approximate the CELL as a perfect sphere. This is absolutely wrong, how to approach surface?

        for atom in gro_list:
            x, y, z = atom[2], atom[3], atom[4]
             
        num_bins = math.ceil(total_area / bin_size) #round up number of bins
    
        bin_counts = np.zeros(num_bins, dtype=int)
    
        # In the next logic, we assign particles to bins based on their position
        for atom in gro_list:
        
            theta = math.acos(atom[2] / bond_radius)
            phi = math.atan2(atom[1], atom[0])
        
            bin_index = int((theta / math.pi) * num_bins)
        
            # If an atom is found, add it to the bin
            bin_counts[bin_index] += 1
    
            # Calculate the average density in each bin
            bin_densities = bin_counts / bin_size
    
        return bin_densities


##UNIT TESTING - CURRENTLY DOES NOT WORK. 
analysis = SingleCellAnalysis()
density = analysis.calculate_particle_density("/wrk/matspunt/projects/cell_model/param_C180/analysis_testing/two_bead_type/sim/CELL.gro", 20, 4, 0.47, 1.85)
