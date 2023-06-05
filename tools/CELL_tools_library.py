import datetime
from distutils.spawn import find_executable
import glob
import json
import os
import subprocess as sp
import sys
import tempfile
import shutil
import numpy as np
import pandas as pd
import csv

####################################################################################################################
#                                           A. DATABASE
# We maintain a central 'database' in the form of a .csv file. We include some checks of the existing entries before
# allowing the user to proceed to prevent mistakes in the simulation definitions, and unnecessarily repeating simulations
####################################################################################################################

import pandas as pd

class DatabaseMaintenance:
    @staticmethod
    def check_repeat_sim(csv_file, dir_name_list):
        
        #with this counter we track the number of overlapping simulations and warn the user
        counter = 0

        with open(csv_file, 'r') as file:
            reader = csv.reader(file)
            for row in reader:
                if len(row) >= 2:
                    dir_name = row[1]
                    if dir_name in dir_name_list:
                        counter += 1
                        # let's introduce an index to deal with repeats. If the user decides to continue, the directory name is 
                        # changed by appending 'repeat{index}' at the end
                        index = 1
                        new_dir_name = f"{dir_name}_repeat{index}"
                        while new_dir_name in dir_name_list:
                            index += 1
                            new_dir_name = f"{dir_name}_repeat{index}"
                        dir_name_list[dir_name_list.index(dir_name)] = new_dir_name

        print(f"{counter} repeating simulations out of {len(dir_name_list)} requested have been found.")
        user_input = input("Repeating simulations will be saved in the .csv with a 'repeatX' suffix. Do you want to proceed? (y/n): ")

        if user_input.lower() == 'y':
            print("Proceeding with the simulation setup!")
            return True
        else:
            print("Stopping the run based on your request...")
            return False

    @staticmethod
    def save_new_entries_db(csv_file, new_entries):
        df = pd.DataFrame(new_entries, columns=['job_name', 'dir_name', 'beads', 'mass', 'LJ_params', 'integrator', 'tcoupl', 'pcoupl'])

        df.to_csv(csv_file, mode='a', header=False, index=False)

    @staticmethod
    def move_files_to_parentdir(source_dir, destination_dir, extensions):
        os.makedirs(destination_dir, exist_ok=True)
        for root, dirs, files in os.walk(source_dir):
            for file in files:
                _, file_ext = os.path.splitext(file)
                if file_ext.lower() in extensions:
                    src_path = os.path.join(root, file)
                    dest_path = os.path.join(destination_dir, file)
                    shutil.move(src_path, dest_path)
    
    @staticmethod          
    def del_files_with_extension(directory, extension):
        for root, dirs, files in os.walk(directory):
            for file in files:
                _, file_ext = os.path.splitext(file)
                if file_ext.lower() == extension:
                    file_path = os.path.join(root, file)
                    os.remove(file_path)

####################################################################################################################
#                                           B. FORCE FIELD PARSER
#            Parses forcefield.itp based on a dictionary containing the bead names, mass and LJ identities of the particles
#                Note: currently the only hardcoded item is the name of the single center bead (hardcoded as 'C').
#                       Otherwise any arbitrary bead names and parameters are accepted by the class. 
#             Any number of membrane beads is allowed, assuming the LJ interaction pairs are defined in the dictionary. 
####################################################################################################################


class ForceFieldParser:
    """
        This class takes a dictionary including information about the force field (but is agnostic
        to the number of beads in the force field, as long as all the LJ and bonded terms are defined)
        For examples of force field dictionaries, check below. 
    """
    def __init__(self, dictionary):
        self.dictionary = dictionary
        self.itp_content = ""

    def generate_atomtypes_section(self):
        """
        For code readability, sections are separated in helper functions which each take care of 
        one component of the 'forcefield.itp'
        """
        self.itp_content += "[ atomtypes ]\n"
        for bead_type, bead_info in self.dictionary['bead_types'].items():
            name = bead_info['name']
            mass = bead_info['mass']
            self.itp_content += f"   {bead_type.ljust(7)}    {mass}     0.0      A      0.0       0.0\n"

    def generate_nonbond_params_section(self):
        self.itp_content += "\n[ nonbond_params ]\n"
        nonbond_params = self.dictionary['nonbond_params']
        for key, values in nonbond_params.items():
            bead_type1, bead_type2 = key.split('_')  # Split over the underscore to obtain the bead types
            sigma = values['sigma']
            epsilon = values['epsilon']
            self.itp_content += f"  {bead_type1}  {bead_type2}    1      {sigma}  {epsilon}\n"

    def generate_bondtypes_section(self):
        self.itp_content += "\n[ bondtypes ]\n"
        for bond_type, bond_info in self.dictionary['bond_types'].items():
            r0 = bond_info['r0']
            fk = bond_info['fk']
            split_string = bond_type.split("_")  # Split over underscore to obtain bead types
            self.itp_content += f"   {split_string[0]}     {split_string[1]}  1     {r0}    {fk}\n"

    def generate_forcefield_itp(self):
        """
        Assembles the force field based on the helper functions and the provided dictionary
        """
        self.itp_content = "[ defaults ]\n"
        self.itp_content += "    1           2          yes          1.0     1.0\n\n"
        self.generate_atomtypes_section() #call helper functions to assemble force field
        self.generate_nonbond_params_section()
        self.generate_bondtypes_section()

        with open("forcefield.itp", "w") as file:
            file.write(self.itp_content)


####################################################################################################################
#                                      C. CONSTRUCTING COORDINATE FILES
#              Note that any number of membrane particles are allowed but that the NAMES
#              of the particles match those in the force field definition (these are necessary to build 
#              the final CELL.top). Also currently hardcoded for a single bead named 'C'
####################################################################################################################


class InputSingleCell():
    """
    This class preprocesses the necessary inputs for creating a single cell using PACKMOL. The inputs from 'packmol_dict'
    are used to create an .xyz file with the necessary number of atoms and packing radius 
    """

    def __init__(self, packmol_dict):
        self.packmol_dict = packmol_dict
        
    def preprocessing_PACKMOL_single_cell(self):
        """
        Saves a .XYZ file for each individual bead in the CELL for PACKMOL input. The name of the beads is read from the 'beads' parameter. 
        PACKMOL needs the actual coordinate files before it can run. This approach is dumb, but effective. 
        
        Then builds the .inp string which contains the PACKMOL input settings based on the settings in the 'input.json'

        Returns:
            inp (str): Packmol input file string for building a single CELL
        """
        #extract PACKMOL settings from the dict
        tolerance = self.packmol_dict['settings']['tolerance']
        radius = self.packmol_dict['settings']['radius']
        shape = self.packmol_dict['settings']['shape']
        
        #create placeholder .xyz files for the different beads. 
        bead_names = [bead_data['name'] for bead_data in self.packmol_dict['beads'].values()]
        
        xyz = "1\n\nX         10.00000       10.00000       10.00000"
        
        #Extract the bead names, and build their .xyz files (these will be removed later)
        for bead in bead_names:
            file_name = f"{bead}.xyz"
            with open(file_name, "w") as f:
                f.write(xyz.replace("X", bead))

        inp = f"tolerance {tolerance}\nfiletype xyz\noutput CELL.xyz\nmovebadrandom\n"
        
        for bead in bead_names: #create the PACKMOL settings for each bead type 
            if bead != "C":
                dr = float(radius) - 0.1 #after testing: for optimal packing, the inner radius should be 0.1 Angstrom smaller than the cell radius
                #the center bead name is assumed to be called as "C.xyz"!!!! We assume a fixed topology
                inp += f"\nstructure {bead}.xyz\n  number {self.packmol_dict['beads'][bead]['number']}\n  atoms 1\n    inside {shape} 0. 0. 0. {radius}\n  end atoms\n  atoms 1\n    outside {shape} 0. 0. 0. {dr}\n  end atoms\nend structure\n\n"
        inp += f"structure C.xyz\n  number 1\n  center\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n"
        
        return inp
    
class PackmolExecuterSingleCell(InputSingleCell):
    """
    This class executes PACKMOL for a single cell using the preprocessed inputs inherited from the "InputSingleCell" class.
    """

    packmol_path = find_executable("packmol") # define packmol_path as a class attribute

    def __init__(self, packmol_dict):
        """
        Initialize the parent class (InputSingleCell) and its attributes.
        """
        super().__init__(packmol_dict) # inherit functionality from "InputSingleCell" preprocessor class

    @staticmethod
    def check_packmol_path():
        """
        Checks if the packmol binary is present in the system and raises an error if not.
        """
        if PackmolExecuterSingleCell.packmol_path is not None:  # access packmol_path using the class name
            print(f"Using PACKMOL from: {PackmolExecuterSingleCell.packmol_path}")
        else:
            print(f"Error: cannot locate PACKMOL binary. Is it installed as 'packmol' and are you in the right environment?\n")
            sys.exit(1)
 
    def run_packmol_single_CELL(self):
        """
        Runs Packmol with the specified input parameters, and builds the coordinate .xyz file. Kills the program if problems in PACKMOL were found. 
        """
        #self.check_packmol_path()
        input_file = self.preprocessing_PACKMOL_single_cell() #inherit the inp from the 'InputSingleCell' class to keep code modular
        # see: https://github.com/mosdef-hub/mbuild/issues/419, using tempfile as workaround
        packmol_inp = tempfile.NamedTemporaryFile(mode="w", delete=False, prefix="packmol-", suffix=".inp")
        packmol_inp.write(input_file)
        packmol_inp.close()
        now = datetime.datetime.now()
        log_file = "PACKMOL_build-{}.log".format(now.strftime("%d-%m-%H-%M-%S"))
        with open(log_file, "w") as f:
            try:
                proc = sp.run("{} < {}".format(self.packmol_path, packmol_inp.name), stdout=sp.PIPE, universal_newlines=True, shell=True, timeout=60)
                stdout = proc.stdout
                f.write(stdout)
                #if 'Success!' in proc.stdout:
                    #print("SUCCES: The coordinate file has been built by PACKMOL. A logfile is saved as:" + log_file)
                # remove all .xyz files except for CELL.xyz
                for file in glob.glob('*.xyz'):
                    if file != 'CELL.xyz':
                        os.remove(file)
                if 'ENDED WITHOUT PERFECT PACKING' in proc.stdout:
                    print("WARNING: PACKMOL has finished but complains about imperfect packing. Please validate the resulting structure, and check the logfile at: " + log_file)
                    print("If the packing is obviously wrong, please validate your input parameters and try again.")
                    # remove all .xyz files except for CELL.xyz
                for file in glob.glob('*.xyz'):
                    if file != 'CELL.xyz':
                        os.remove(file)
            except sp.TimeoutExpired:
                print("ERROR: Packmol took longer than 60 seconds to run, this is highly unusual and suggests the packing could not be resolved.")
                print("Most likely the cell radius and/or the number of particles are incorrectly defined in your input.JSON")
                for file in glob.glob('*.xyz'):
                        os.remove(file)
 
####################################################################################################################
#                                           D. GRO-CONVERTER & ITP PARSER
#                       ITPs are parsed based on the constructed coordinates and force field!
# 
####################################################################################################################

class Gromacs_IO:
    """
    Methods:
    --------
    convert_xyz_to_gro(xyz_file, gro_file):
        Converts an XYZ file to a GRO file         

    build_GMX_top(gro_path):
        Builds a GROMACS topology file (.itp) based on a given coordinate file (in .gro format)
        and the force field that is supplied. 
    """

    @staticmethod
    def check_GMX_path():
        """
        Checks if the GROMACS binary (gmx) is present in the system and returns its path.
        """
        gmx_path = find_executable("gmx")  # check whether GROMACS is available in the environment
        if gmx_path is not None:
            print(f"Using GROMACS version from: {gmx_path}\n")
            return gmx_path
        else:
            print(f"\n\nError: cannot find GROMACS (gmx) binary. Is it installed and sourced correctly? \n")
            sys.exit(1)

    @staticmethod #make it static so we can call it without instancing the class. 
    def convert_xyz_to_gro(xyz_file, gro_file, box_size):
        """Converts an XYZ file to a GRO file (hardcodes resname as 'CELL' but is otherwise non-specific)

        Args:
            xyz_file (str): The path to the input XYZ file.
            gro_file (str): The path to the output GRO file.
            box_size (float): The size of the cubic box.
    """
        try:
            with open(xyz_file, 'r') as xyz:
                # skip first two lines
                xyz.readline()
                xyz.readline()

                atoms = []
                for line in xyz:
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
                gro.write(f'{i + 1:>5d}CELL{name:>5s}{i + 1:>5d}{x:>8.3f}{y:>8.3f}{z:>8.3f}\n')  # .GRO only saves 3 decimals (single precision)

            gro.write(f'{box_size:>10.5f}{box_size:>10.5f}{box_size:>10.5f}\n')  # add the box at the end of the .gro file
            # note that due to the XYZ --> GRO conversion the CELL is centered at the origin (0,0,0) and not the box center
            
    @staticmethod
    def build_GMX_top_single_CELL(gro_path, forcefield_path):
        """
        Builds itp based on a given coordinate file (in .gro format) and force field file.
        Currently assumes the Center bead is at the end of the .gro (in how the bonds are drawn),
        but this can be easily made general if needed.
        """

        with open(gro_path, "r") as gro:
            abspath = os.path.abspath(gro_path)
            gro_list = [line.split() for line in gro.readlines()[2:-1]]  # Read and split lines, excluding the first and last lines

            with open("CELL.top", "w") as top:
                # Write the topology header
                now = datetime.datetime.now()
                header = "; Topology file generated from {} at {}\n".format(abspath, now.strftime("%H:%M:%S"))
                top.write(header)

                # Read the force field file and extract bond types
                bondtypes = {}
                with open(forcefield_path, "r") as forcefield:
                    lines = forcefield.readlines()

                    # Extract bond types from [bondtypes] section
                    bondtypes_section = False
                    for line in lines:
                        line = line.strip()
                        if line == "[ bondtypes ]":
                            bondtypes_section = True
                        elif bondtypes_section and line:
                            atom_type1, atom_type2, func, r0, fk = line.split()
                            bondtypes[(atom_type1, atom_type2)] = {
                                "func": int(func),
                                "r0": float(r0),
                                "fk": float(fk),
                            }

                # Write the forcefield that is used
                ff_itp = f"\n#include \"{forcefield_path}\"\n"
                top.write(ff_itp)

                # Write the [moleculetype] directive and the [atoms] directive based on the .GRO
                top.write("\n[ moleculetype ]\n; Name        nrexcl\n  CELL        1\n\n[ atoms ]\n; nr type resnr residue atom cgnr  charge\n")
                for lines in gro_list:
                    atom_type, atom_nr = str(lines[1]), str(lines[2])  # save the atom name and index from the .GRO
                    atom_name = atom_type[0]
                    # we don't need to parse the mass explicitly, since GMX will take it from our force field, but this can be easily added, if needed
                    top.write("  {:<3s}  {:<3s}  1    CELL    {:<3s}  {:<3s}  0.0000 \n".format(atom_nr, atom_type, atom_name, atom_nr))

                # Write the [ bonds ] directive based on the force field and .GRO file
                top.write("\n[ bonds ]\n; i j func  r0 fk\n")
                # we extract the atom names (except Center bead) from the bondtypes dictionary
                atom_names = set()
                for key in bondtypes.keys():
                    atom_names.add(key[1])

                # now we look for the atom in the .gro and parse the right bonded type accordingly
                for lines in gro_list:
                    if lines[1] in atom_names:
                        atom_type, atom_nr = str(lines[1]), str(lines[2])  # save the atom name and index from the .GRO
                        n_atom = str(len(gro_list))  # if the last atom is always the C-bead, we can easily find it like this
                        values = bondtypes[('C', atom_type)]
                        top.write("  {:<2s} {:<3s} {:<3s} {:<3s} {:<3s} \n".format(atom_nr, n_atom, str(values['func']), str(values['r0']), str(values['fk'])))

                # Write system information for simulation at the end
                # Count the number of cells from the provided .gro file by counting the occurrences of the C-bead
                # This is a bit hacky but minimizes the amount of parameters the user neeeds to input
                nr_cells = sum(1 for line in gro_list if line[1] == 'C')
                system = f"\n[ system ]\nCELL MODEL\n\n[ molecules ]\nCELL    {nr_cells}"
                top.write(system)

####################################################################################################################
#                                           D. MDP parser
#                                           Parsing of MDPS
# 
####################################################################################################################

class Gromacs_MDP:
    @staticmethod
    def write_min_mdp_file(filename, em_algorithm, nsteps, emtol, emstep):
            #we don't need to change any other minimization settings than the ones exposed here. 
        min_mdp = {
            'integrator': em_algorithm,
            'nsteps': nsteps,
            'emtol': emtol,
            'emstep': emstep,
            'nstcomm': '100',
            'nstxout': '0',
            'nstvout': '0',
            'nstfout': '0',
            'nstlog': '1000',
            'nstenergy': '100',
            'nstxout-compressed': '1000',
            'compressed-x-precision': '100',
            'compressed-x-grps': 'System',
            'energygrps': 'System',
            'cutoff-scheme': 'Verlet',
            'nstlist': '20',
            'ns_type': 'grid',
            'pbc': 'xyz',
            'verlet-buffer-tolerance': '0.005',
            'coulombtype': 'cutoff',
            'coulomb-modifier': 'Potential-shift-verlet',
            'rcoulomb': '1.1',
            'epsilon_r': '15',
            'vdw_type': 'cutoff',
            'vdw-modifier': 'Potential-shift-verlet',
            'rvdw': '1.1',
            'tcoupl': 'no',
            'pcoupl': 'no',
            'constraints': 'none',
            }
        
        with open(filename, "w") as f:
            for key, value in min_mdp.items():
                f.write(f"{key} = {value}\n")
            
        #if not os.path.exists('toppar'):
            #os.makedirs('toppar')
            #shutil.move(filename, 'toppar/'+filename)
        #else:
            #shutil.move(filename, 'toppar/'+filename)       
              
    @staticmethod
    def write_run_mdp_file(filename, nsteps, params_dict):

        run_mdp = {
            'integrator': 'md',
            'dt': '0.02',
            'nsteps': nsteps,
            'nstcomm': '100',
            'nstxout': '0',
            'nstvout': '0',
            'nstfout': '0',
            'nstlog': '1000',
            'nstenergy': '100',
            'nstxout-compressed': '1000',
            'compressed-x-precision': '100',
            'compressed-x-grps': 'System',
            'energygrps': 'System',
            'cutoff-scheme': 'Verlet',
            'nstlist': '20',
            'ns_type': 'grid',
            'pbc': 'xyz',
            'verlet-buffer-tolerance': '0.005',
            'coulombtype': 'cutoff',
            'coulomb-modifier': 'Potential-shift-verlet',
            'rcoulomb': '1.1',
            'epsilon_r': '15',
            'vdw_type': 'cutoff',
            'vdw-modifier': 'Potential-shift-verlet',
            'rvdw': '1.1',
            'tc-grps': 'System',
            'tau_t': '1.0',
            'ref_t': '310',
            'Pcoupltype': 'isotropic',
            'tau_p': '12.0',
            'compressibility': '3e-4',
            'ref_p': '1',
            'constraints': 'none',
            'constraint_algorithm': 'Lincs',
            'gen_vel': 'yes',
            'gen_seed': '-1'
        }
        
        #where params_dict is a dict with the desired change in parameters which overrules the existing .mdp
        run_mdp.update(params_dict)
        
        with open(filename, "w") as f:
            for key, value in run_mdp.items():
                f.write(f"{key} = {value}\n")
        
        #move everything to a folder "toppar" to keep .itp and .mdp nicely contained. 
        #if not os.path.exists('toppar'):
            #os.makedirs('toppar')
            #shutil.move(filename, 'toppar/'+filename)
        #else:
            #shutil.move(filename, 'toppar/'+filename)


####################################################################################################################
#                                           E. RUNNING GROMACS SIMULATIONS
#                                           The following section contains the 'gmx grompp' and
#                                           'gmx mdrun' wrappers required to run simulations with the params
####################################################################################################################

class GromacsRun:
    @staticmethod
    def run_GMX_basic(filename='CELL.gro', timeout=2400, nr_of_threads=12, top_name='CELL.top'):
        """
            Basic logic to start and run GROMACS. Timeout = 2400 s, is appropriate for 
            about a microsecond simulation to deal with crashes/divisions by zero.
            Assumes single minimization and single run step. 
            Assumes 'min.mdp' and 'run.mdp' name for mdps
        """
        GMX = find_executable("gmx")
        if not GMX:
            raise RuntimeError("Cannot find GROMACS executable 'gmx' in PATH")
        try:
            # we need to extract the box size to obtain the input for 'gmx editconf'
            # so the CELL is centered in the box
            # This is shitty but it works...
            with open(filename, 'r') as gro:
                lines = gro.readlines()
                box_size = lines[-1].strip()
            
            editconf= f'gmx editconf -f {filename} -o {filename} -box {box_size} >/dev/null 2>&1'
            os.system(editconf)
            
            grompp_min_cmd = ['gmx', 'grompp', '-p', top_name, '-f', 'min.mdp', '-c', filename, '-o', '1-min', '-maxwarn', '1']
            sp.run(grompp_min_cmd, check=True, stdout=open('gmx_run.log', 'a'), stderr=sp.STDOUT)
            
            mdrun_min_cmd = ['gmx', 'mdrun', '-nt', str(nr_of_threads), '-pin', 'on', '-deffnm', '1-min', '-v']
            sp.run(mdrun_min_cmd, check=True, timeout=timeout, stdout=open('gmx_run.log', 'a'), stderr=sp.STDOUT)

            grompp_run_cmd = ['gmx', 'grompp', '-p', top_name, '-f', 'run.mdp', '-c', '1-min.gro', '-o', '2-run', '-maxwarn', '1']
            sp.run(grompp_run_cmd, check=True, stdout=open('gmx_run.log', 'a'), stderr=sp.STDOUT)

            mdrun_run_cmd = ['gmx', 'mdrun', '-nt', str(nr_of_threads), '-pin', 'on', '-deffnm', '2-run', '-v']
            sp.run(mdrun_run_cmd, check=True, timeout=timeout, stdout=open('gmx_run.log', 'a'), stderr=sp.STDOUT)

        except sp.TimeoutExpired:
            print('The simulation timed out and was terminated. Log saved as "gmx_run.log"')
            return
        except sp.CalledProcessError as e:
            print(f'The simulation failed with error {e.returncode}, most likely there are problems in the input parameters. Check the log at "gmx_run.log"')
            return
        else:
            print(f'This simulation completed without any errors, the log information is saved in "gmx_run.log"\n')
        