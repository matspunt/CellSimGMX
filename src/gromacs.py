import os
import datetime
import subprocess as sp
from distutils.spawn import find_executable
import shutil
import sys

from settings_parser import JSONParser
from settings_parser import ForcefieldParserGMX

class Gromacs_IO(JSONParser, ForcefieldParserGMX):
    """
    Methods:
    --------
    
    estimate_box_size_xyz(self, xyz_file):
        Estimates the box size based on coordinates in a XYZ file.
    
    convert_xyz_to_gro(self, xyz_file, gro_file, edge_offset):
        Converts an XYZ file to a GRO file with a user specified edge offset, 
        and centers its coordinates in the middle of a (cubic) box
        
    def build_GMX_top_single_CELL(self, gro_path):
        Builds .itp from a ssingle .gro file based on the JSON input settings
    """

    def __init__(self, json_directory, forcefield_directory):
        JSONParser.__init__(self, json_directory)
        ForcefieldParserGMX.__init__(self, forcefield_directory)

    def estimate_box_size_xyz(self, xyz_file):
        """
        Estimates the required box size based on the coordinates in the XYZ file. This finds
        the tightest box that can be fitted based on the min, max coordinates in each
        direction. Only supports cubic boxes. 

        Returns:
            Tuple containing the box size in nm for each dimension (x, y, z).
        """
        try:
            with open(xyz_file, 'r') as xyz:
                # Skip the first two lines
                next(xyz)
                next(xyz)

                min_x = min_y = min_z = float("inf")
                max_x = max_y = max_z = float("-inf")
                for line in xyz:
                    columns = line.split()
                    name = columns[0]  # atom_names
                    
                    x, y, z = map(float, columns[1:]) #coordinates

                    # Find minimum and maximum coordinates to figure out where box needs to go
                    min_x = min(min_x, x)
                    min_y = min(min_y, y)
                    min_z = min(min_z, z)
                    max_x = max(max_x, x)
                    max_y = max(max_y, y)
                    max_z = max(max_z, z)

            # Calculate the optimal box size in nm (divide by 10!) based on min and max coordinates
            box_size_x = abs(max_x - min_x) / 10.0
            box_size_y = abs(max_y - min_y) / 10.0
            box_size_z = abs(max_z - min_z) / 10.0

            return box_size_x, box_size_y, box_size_z

        except FileNotFoundError:
            raise Exception(f"ERROR: Input file '{xyz_file}' not found. Please check the file path and try again.")
        
    def convert_xyz_to_gro(self, xyz_file, gro_file, edge_offset=3.0):
        """
        Converts an XYZ file to a GRO file and positions the atoms in the new box, then adjusts
        the box size based on the user-defined offset to deal with PBC problems. 

        Args:
            edge_offset(float): Set to 0, the coordinates would perfectly clip the box edges. This offset mimics the
            behavior of the smallest distance to box edges in CHARMM-GUI. The default is set to 2.5 nm.
        """
        try:
            with open(xyz_file, 'r') as xyz:
                # Skip the first two lines
                next(xyz)
                next(xyz)

                atoms = [line.split() for line in xyz]

        except FileNotFoundError:
            raise Exception(f"ERROR: Input file '{xyz_file}' not found. Please check the file path and try again.")

        # Calculate geometric center coordinates in x,y and z directions, this is a bit ugly but it works :D
        total_coords = len(atoms)
        sum_x = sum(float(coord[1]) for coord in atoms)
        sum_y = sum(float(coord[2]) for coord in atoms)
        sum_z = sum(float(coord[3]) for coord in atoms)
        center_x = sum_x / total_coords / 10.0
        center_y = sum_y / total_coords / 10.0
        center_z = sum_z / total_coords / 10.0

        box_size = self.estimate_box_size_xyz(xyz_file) #estimate tightest fitting cubic box
        
        # Calculate the box size with the specified edge offset
        new_box_size = (
            box_size[0] + edge_offset,
            box_size[1] + edge_offset,
            box_size[2] + edge_offset
        )

        # Note: (0,0,0) is recognized as the box origin, we thus need to translate the coordinates based
        # on the box dimensions! Calculate the ratio of translation of the positions to the geometric 
        # center of the XYZ coordinates
        trans_x = -(center_x - (new_box_size[0] / 2))
        trans_y = -(center_y - (new_box_size[1] / 2))
        trans_z = -(center_z - (new_box_size[1] / 2))
        
        # Move the atoms based on the calculated translation factor
        translated_atoms = []
        for atom in atoms:
            name, x, y, z = atom
            x = (float(x) / 10.0) + trans_x
            y = (float(y) / 10.0) + trans_y
            z = (float(z) / 10.0) + trans_z
            translated_atoms.append((name, x, y, z))

        with open(gro_file, 'w') as gro:
            now = datetime.datetime.now()
            header = "GRO file built from {} at {}\n".format(xyz_file, now.strftime("%H:%M:%S"))
            gro.write(header)
            gro.write(f'{len(translated_atoms)}\n')

            for i, atom in enumerate(translated_atoms):
                name, x, y, z = atom
                gro.write("{:>5d}CELL{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n".format(i + 1, name, i + 1, x, y, z))

            gro.write("{:>10.5f}{:>10.5f}{:>10.5f}\n".format(*new_box_size))

        print(f"SUCCESS: The GRO file has been created as '{gro_file}'.")

    def build_GMX_top_single_CELL(self, gro_path):
        """
            Builds itp based on a given coordinate file (in .gro format). Currently assumes the Center bead is at the end of the .gro
            (in how the bonds are drawn) but this can be easily made general if needed. 
        """
        #access required information from the input files to build topology
        self.parse_GMX_ff()
        self.load_json_CELL()
        self.extract_inputs_JSON()
        
        with open(gro_path, "r") as gro:
            abspath = os.path.abspath(gro_path)
            gro_list = [line.split() for line in gro.readlines()[2:-1]]  # Read and split lines, excluding the first and last lines

            with open("CELL.top", "w") as top:
                #Write the topology header
                now = datetime.datetime.now()
                header = "; Topology file generated from {} at {}\n".format(abspath, now.strftime("%H:%M:%S"))
                top.write(header)

                #write the forcefield that is used
                ff_itp =  f"\n; Using forcefield from:\n#include \"{self.itp_path}\"\n"
                top.write(ff_itp)
                
                # Write the [moleculetype] directive and the [atoms] directive based on the .GRO
                top.write("\n[ moleculetype ]\n; Name        nrexcl\n  CELL        1\n\n[ atoms ]\n; nr type resnr residue atom cgnr  charge\n")
                for lines in gro_list:
                    atom_type, atom_nr = str(lines[1]), str(lines[2]) #save the atom name and index from the .GRO
                    atom_name = atom_type[0]
                    #we don't need to parse the mass explicitly, since GMX will take it from our force field, but this can be easily added, if needed
                    top.write("  {:<3s}  {:<3s}  1    CELL    {:<3s}  {:<3s}  0.0000 \n".format(atom_nr, atom_type, atom_name, atom_nr))
                    
                # Write the [ bonds ] directive based on the force field and .GRO file 
                top.write("\n[ bonds ]\n; i j func  r0 fk\n")
                #we extract the atom names (except Center bead) from the bondedtypes dictionary
                atom_names = set()
                for key in self.bondtypes.keys():
                    atom_names.add(key[1])
                    
                #now we look for the atom in the .gro and parse the right bonded type accordingly
                for lines in gro_list:
                    if lines[1] in atom_names:
                        atom_type, atom_nr = str(lines[1]), str(lines[2]) #save the atom name and index from the .GRO
                        n_atom = str(len(gro_list))  # if the last atom is always the C-bead, we can easily find it like this
                        values = self.bondtypes[('C', atom_type)] 
                        top.write("  {:<2s} {:<3s} {:<3s} {:<3s} {:<3s} \n".format(atom_nr, n_atom, str(values['func']), str(values['r0']), str(values['fk'])))
                        
                # Write system information for simulation at the end (I see no reason why to separate the topology in a .itp and .top file)
                system =  f"\n[ system ]\nCELL MODEL\n\n[ molecules ]\nCELL    {self.Simulation_nr_cells}" 
                top.write(system)
            top.close()
            print("A topology file has been constructed")

### WORKING WITH THE FUNCTIONS
json_directory = "/wrk/matspunt/coding/CELL_MODEL/src"
forcefield_directory = "/wrk/matspunt/coding/CELL_MODEL/src"
gromacs_io = Gromacs_IO(json_directory, forcefield_directory)

#gromacs_io.convert_xyz_to_gro("CELL.xyz", "CELL.gro", edge_offset=3.0)
#gromacs_io.convert_xyz_to_gro("CELL_27_standard.xyz", "CELL_27_standard.gro", edge_offset=3.0)
gromacs_io.convert_xyz_to_gro("CELL_8_standard.xyz", "CELL_8_standard.gro", edge_offset=3.0)
#gromacs_io.convert_xyz_to_gro("CELL_2_monolayer.xyz", "CELL_2_monolayer.gro", edge_offset=3.0)


