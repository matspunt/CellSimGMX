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
    
    convert_xyz_to_gro(xyz_file, gro_file):
        Converts an XYZ file to a GRO file         

    build_GMX_top(gro_path):
        Builds a GROMACS topology file (.itp) based on a given coordinate file (in .gro format)
        and the force field that is supplied. 
    """

    def __init__(self, json_directory, forcefield_directory):
        JSONParser.__init__(self, json_directory)
        ForcefieldParserGMX.__init__(self, forcefield_directory)
        self.GMX_path = find_executable("gmx") #check whether GROMACS is available in the environment
        
        if self.GMX_path is not None:
            print(f"\n\nUsing GROMACS version from: {self.GMX_path} \n")
        else:
            print(f"\n\nError: cannot find GROMACS (gmx) binary. Is it installed and sourced correctly? \n")
            sys.exit(1)
    
    @staticmethod #make it static so we can call it without instancing the class. 
    def convert_xyz_to_gro(xyz_file, gro_file):
        """Converts an XYZ file to a GRO file (hardcodes resname as 'CELL' but is otherwise non-specific)

        Args:
            xyz_file (str): The path to the input XYZ file.
            gro_file (str): The path to the output GRO file.
        """
        try:
            with open(xyz_file, 'r') as file:
                # skip first two lines
                file.readline()
                file.readline()

                atoms = []
                for line in file:
                    name = line[0:4].strip()
                    x, y, z = map(float, line[13:].split())
                    atoms.append((name, x, y, z))
        except FileNotFoundError:
            raise Exception(f"ERROR: input file '{xyz_file}' not found. Please check that PACKMOL built a suitable coordinate file and try again.  ")

        with open(gro_file, 'w') as gro:
            now = datetime.datetime.now()
            header = "GRO file built from {} at {}\n".format(xyz_file, now.strftime("%H:%M:%S"))
            gro.write(header)
            gro.write(f'  {len(atoms)}\n')

            for i, atom in enumerate(atoms):
                name, x, y, z = atom
                x /= 10.0  # .xyz file format works with Angstroms instead of nanometers
                y /= 10.0
                z /= 10.0
                gro.write(f'{i + 1:>5d}CELL{name:>5s}{i + 1:>5d}{x:>8.3f}{y:>8.3f}{z:>8.3f}\n') #.GRO only saves 3 decimals (single precision)

            gro.write(f'{0:>10.5f}{0:>10.5f}{0:>10.5f}\n') #add an empty box at the end of the .GRO

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
            


