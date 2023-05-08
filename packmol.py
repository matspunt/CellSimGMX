import datetime
from distutils.spawn import find_executable
import glob
import json
import os
import subprocess as sp
import sys
import tempfile

class PackmolInput:
    """
    A class for reading input options for Packmol simulation from a standard JSON file.

    Methods:
    --------
    __init__(self, filename):
        Initializes the class, actually reads the file. 
    construct_bead_xyz(bead_identities)
        Based on the bead 'force field', constructs the XYZ files required for PACKMOL
    """
    
    def __init__(self, filename):
        with open(filename) as f:
            data = json.load(f)
        self.tolerance = data["packmol_tolerance"]
        self.number_of_cells = (data["number_of_cells"])
        self.beads = data["beads"]
        self.init_shape = data["initial_packing_shape"]
        if self.init_shape not in ["sphere", "cube", "ellipsoid"]:
            sys.exit("Error: initial_packing_shape must be 'sphere', 'cube', or 'ellipsoid'")
        self.cell_radius = float(data["cell_radius"])
    
    def construct_bead_xyz(self):
        """
        Saves a .XYZ file for each individual bead in the CELL for PACKMOL input.
        The name of the beads is read from the 'beads' parameter.
        """
        xyz = "1\n\nX         10.00000       10.00000       10.00000"
        bead_identities = []
        for bead, count in self.beads.items():
            for i in range(count):
                bead_identities.append(bead)
        for bead in bead_identities:
            file_name = f"{bead}.xyz"
            with open(file_name, "w") as f:
                f.write(xyz.replace("X", bead))
        
    def prepare_inp(self):
        """
        Prepares a Packmol input file formatted based on the input.json file.

        Returns:
            inp (str): Packmol input file string. Is used in PackmolExecuter.run_packmol()
        """
        inp = f"tolerance {self.tolerance}\nfiletype xyz\noutput CELL.xyz\nmovebadrandom\n"
        for bead, number in self.beads.items():
            if bead != "N":
                dr = self.cell_radius - 0.1
                inp += f"\nstructure {bead}.xyz\n  number {number}\n  atoms 1\n    inside {self.init_shape} 0. 0. 0. {self.cell_radius}\n  end atoms\n  atoms 1\n    outside {self.init_shape} 0. 0. 0. {dr}\n  end atoms\nend structure\n\n"
        inp += f"structure N.xyz\n  number {self.beads['N']}\n  center\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n"
        return inp

class PackmolExecuter:
    """
    This class inherits PACKMOL input settings from the "PackmolInput" class and uses it to run PACKMOL. 

    Methods:
    ----------
    check_packmol_path():
        Checks if the packmol binary is present in the system and raises an error if not.
    run_packmol():
        Runs Packmol with the specified input parameters and builds the single CELL model based on JSON input. 
    """
    
    def __init__(self):
        self.packmol_path = find_executable("packmol")
        json_file = str(glob.glob("*json"))[2:-2] #look for a JSON file in the given directory
        self.packmol_input = PackmolInput(json_file) #inherit from the PackmolInput class to read the JSON file
    
    def check_packmol_path(self):
        """
        Checks if the packmol binary is present in the system and raises an error if not.
        """
        if self.packmol_path is not None:
            print(f"PACKMOL located at: {self.packmol_path} \n")
        else:
            print(f"Error: cannot locate PACKMOL binary. Is it installed as 'packmol' and are you in the right environment?\n")
            sys.exit(1)
 
    def run_packmol(self):
        """
        Runs Packmol with the specified input parameters, and builds the coordinate .xyz file. Kills the program if problems in PACKMOL were found. 
        """
        self.check_packmol_path()
        self.packmol_input.construct_bead_xyz()
        input_file = self.packmol_input.prepare_inp() #inherit the inp from the 'PackmolInput' class
        # see: https://github.com/mosdef-hub/mbuild/issues/419, use tempfile as workaround
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
                if 'Success!' in proc.stdout:
                    print("\n\nSUCCES: The coordinate file has been built as 'CELL.xyz'. A logfile is saved at: " + log_file)
                # remove all .xyz files except for CELL.xyz
                for file in glob.glob('*.xyz'):
                    if file != 'CELL.xyz':
                        os.remove(file)
                if 'ENDED WITHOUT PERFECT PACKING' in proc.stdout:
                    print("\n\nWARNING: PACKMOL has finished but complains about imperfect packing. Please validate the resulting structure 'CELL.xyz_FORCED', and check the logfile at: " + log_file)
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

#We can read the JSON with the PACKMOL input class
execute_packmol = PackmolExecuter()
execute_packmol.run_packmol()

