import datetime
from distutils.spawn import find_executable
import glob
import json
import os
import subprocess as sp
import sys
import tempfile

from settings_parser import JSONParser
from settings_parser import ForcefieldParserGMX

class InputSingleCell(JSONParser, ForcefieldParserGMX):
    """
    This class preprocesses the necessary inputs for creating a single cell using PACKMOL. The inputs from both the .JSON file and the forcefield
    are used to create an .xyz file with the necessary number of atoms etc. Currently hardcoded for 180 atoms. 
    """

    def __init__(self, json_directory, forcefield_directory):
        """
        Initialize JSONParser and ForcefieldParserGMX so its attributes become a member of 'self'
        """
        JSONParser.__init__(self, json_directory)
        ForcefieldParserGMX.__init__(self, forcefield_directory)
       
    def preprocessing_PACKMOL_single_cell(self):
        """
        Saves a .XYZ file for each individual bead in the CELL for PACKMOL input. The name of the beads is read from the 'beads' parameter. 
        PACKMOL needs the actual coordinate files before it can run. This approach is dumb, but effective. 
        
        Then builds the .inp string which contains the PACKMOL input settings based on the settings in the 'input.json'

        Returns:
            inp (str): Packmol input file string for building a single CELL
        """
        xyz = "1\n\nX         10.00000       10.00000       10.00000"
        self.parse_GMX_ff()
        #Now extract the bead names from the keys in the 'atomtypes' dict
        bead_names = list(self.atomtypes.keys())
        
        # and save placeholder .XYZ files for PACKMOL
        for bead in bead_names:
            file_name = f"{bead}.xyz"
            with open(file_name, "w") as f:
                f.write(xyz.replace("X", bead))
                
        #parse the JSON information and then construct the PACKMOL.inp
        self.load_json_CELL()
        self.extract_inputs_JSON()

        inp = f"tolerance {self.PACKMOL_tol}\nfiletype xyz\noutput CELL.xyz\nmovebadrandom\n"
        
        for bead in bead_names: #create the PACKMOL settings for each bead type 
            if bead != "C":
                dr = self.PACKMOL_cell_radius - 0.1 #for optimal packing, the inner radius should be 0.1 Angstrom smaller than the cell radius
                #the center bead name is assumed to be called as "C.xyz"!
                inp += f"\nstructure {bead}.xyz\n  number 90\n  atoms 1\n    inside {self.PACKMOL_init_shape} 0. 0. 0. {self.PACKMOL_cell_radius}\n  end atoms\n  atoms 1\n    outside {self.PACKMOL_init_shape} 0. 0. 0. {dr}\n  end atoms\nend structure\n\n"
        inp += f"structure C.xyz\n  number 1\n  center\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n"
        
        return inp

class PackmolExecuterSingleCell(InputSingleCell):
    """
    This class executes PACKMOL for a single cell using the preprocessed inputs inherited from the "InputSingleCell" class.
    """

    def __init__(self, json_directory, forcefield_directory):
        """
        Initialize the parent class (InputSingleCell) and its attributes.
        """
        super().__init__(json_directory, forcefield_directory) #inherit functionality from "InputSingleCell" preprocessor class
        self.packmol_path = find_executable("packmol") #check whether PACKMOL is available in the environment the user is calling the programme. 

    def check_packmol_path(self):
        """
        Checks if the packmol binary is present in the system and raises an error if not.
        """
        if self.packmol_path is not None:
            print(f"Using PACKMOL from: {self.packmol_path} \n")
        else:
            print(f"Error: cannot locate PACKMOL binary. Is it installed as 'packmol' and are you in the right environment?\n")
            sys.exit(1)
 
    def run_packmol_single_CELL(self):
        """
        Runs Packmol with the specified input parameters, and builds the coordinate .xyz file. Kills the program if problems in PACKMOL were found. 
        """
        self.check_packmol_path()
        input_file = self.preprocessing_PACKMOL_single_cell() #inherit the inp from the 'InputSingleCell' class
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


