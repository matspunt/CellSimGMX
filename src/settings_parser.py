import json
import glob
import sys
import os

class JSONParser:
    """
    This class handles parsing of the input.JSON file. An instance of this class can be passed to other classes to access the 
    simulation settings. 
    """
    def __init__(self, directory):
        self.directory = directory
        self.file_path = None
        self.data = None

    def load_json_CELL(self):
        """
        Any filename for the .JSON is accepted, as long as it has the right extension. Currently an explicit directory needs to
        be defined for this function to work. Ideally this can be configured through a flag at a later stage. 
        """
        json_files = glob.glob(f"{self.directory}/*.json")
        if len(json_files) > 1:
            print("Multiple .JSON files were found at " + self.file_path + " which is causing conflicts. Remove the wrong ones and try again.") 
        if len(json_files) > 0: 
            self.file_path = json_files[0]
            try:
                with open(self.file_path, 'r') as file:
                    self.data = json.load(file)
                    #print("Simulation settings were succesfully loaded from:" + self.file_path)
            except json.JSONDecodeError:
                print("The supplied JSON file contains errors, please check it complies with the template supplied with the CELL package. ")
            except json.FileNotFoundError:
                print("No input.JSON could be found in the defined directory.")
        else:
            print("No input.JSON could be found in the defined directory.")

    def extract_inputs_JSON(self):
        """
        Extracts the simulation input values from the loaded JSON data. Note that this might require some recasting into alternative datatypes!!! 

        The input values are extracted and stored as attributes of the class. This section will obviously change a lot as we refine the project/setup
        so there are only a few options added now. Fortunately, reading JSON files is super verbose so this is easy to expand on. 
        """
        if self.data is not None:
            #the JSON is divided in main categories that are used in different parts of the code (i.e. PACKMOL, GROMACS)
            #so it makes sense to parse them under separate denominators. This is overly verbose but that's probably a good thing. 
            PACKMOL_data = self.data.get("PACKMOL")
            if PACKMOL_data:
                self.PACKMOL_tol = PACKMOL_data.get("packmol_tolerance", {}).get("value")
                self.PACKMOL_cell_radius = PACKMOL_data.get("cell_radius", {}).get("value")
                self.PACKMOL_init_shape = PACKMOL_data.get("initial_packing_shape", {}).get("value")
                #let's explicitly check for the right shape because typos can easily occur here. Other settings are more ambiguous/up to user. 
                if self.PACKMOL_init_shape not in ["sphere", "cube", "ellipsoid"]:
                    sys.exit("Error: initial_packing_shape must be 'sphere', 'cube', or 'ellipsoid'")
                #print("PACKMOL_cell_radisu:", type(self.PACKMOL_cell_radius))

            Simulation_data = self.data.get("Simulation")
            if Simulation_data:
                self.Simulation_nr_cells = Simulation_data.get("number_of_cells", {}).get("value")

            GROMACS_data = self.data.get("GROMACS")
            if GROMACS_data:
                self.GROMACS_dt = GROMACS_data.get("timestep", {}).get("value")
        else:
            print("Please try again, Something went wrong with loading the .JSON file")

class ForcefieldParserGMX:
    """
    This class handles parsing of the predefined forcefield in GROMACS format ("forcefield.itp").
    Currently using dicts to store the .itp information, but perhaps this is not ideal. An instance of
    this class can be parsed to "Gromacs_IO" class to access force field information. 

    Args:
        directory (str): Directory where the force field is located. 

    Attributes:
        atomtypes (dict)
        nonbond_params (dict)
        bondtypes (dict)
    """
    
    def __init__(self, directory):
        self.directory = directory
        self.itp_path = None
        self.atomtypes = {}
        self.nonbond_params = {}
        self.bondtypes = {}

    def parse_GMX_ff(self):
        """
            First looks for a GMX compatible .itp file in the user-specific directory. Then uses helper functions to parse the 
            different .itp entries into dictionaries. 
        """
        try:
            itp_files = glob.glob(os.path.join(self.directory, "*.itp"))
            if len(itp_files) > 0:
                self.itp_path = itp_files[0]
                #print("Using the force field from:", self.itp_path)
            else:
                raise FileNotFoundError("ERROR: No forcefield.itp file was found in the indicated directory, exiting the programme. ")
        except Exception as e:
            print(e)
            sys.exit(1)
            
        #open the force field .itp and extract the different '[ ]' categories
        with open(self.itp_path, 'r') as itp:
            section = None
            for line in itp:
                line = line.strip()
                if line.startswith('[') and line.endswith(']'):
                    section = line[1:-1].strip()
                elif section == 'atomtypes':
                    self.parse_atomtype_GMX_helper(line)
                elif section == 'nonbond_params':
                    self.parse_nonbond_params_GMX_helper(line)
                elif section == 'bondtypes':
                    self.parse_bondtype_GMX_helper(line)

    def parse_atomtype_GMX_helper(self, line):
        """
            This function is required for PACKMOL. 
        """
        if not line.startswith(';'):
            columns = line.split()
            if len(columns) >= 6:
                atomtype = columns[0]
                mass = float(columns[1])
                charge = float(columns[2])
                ptype = columns[3]
                sigma = float(columns[4])
                epsilon = float(columns[5])
                self.atomtypes[atomtype] = {'mass': mass, 'charge': charge, 'ptype': ptype, 'sigma': sigma, 'epsilon': epsilon}

    def parse_nonbond_params_GMX_helper(self, line):
        """
            This function is probably not necessary as we don't need the LJ parameters to create topologies
            But I am leaving it in if we want to modify force field definitions on the fly. 
        """
        #skip lines with comments
        if not line.startswith(';'):
            columns = line.split()
            if len(columns) >= 5:
                i = columns[0]
                j = columns[1]
                func = int(columns[2])
                sigma = float(columns[3])
                epsilon = float(columns[4])
                self.nonbond_params[(i, j)] = {'func': func, 'sigma': sigma, 'epsilon': epsilon}

    def parse_bondtype_GMX_helper(self, line):
        #skip lines with comments
        if not line.startswith(';'):
            columns = line.split()
            if len(columns) >= 5:
                i = columns[0]
                j = columns[1]
                func = int(columns[2])
                r0 = float(columns[3])
                fk = float(columns[4])
                self.bondtypes[(i, j)] = {'func': func, 'r0': r0, 'fk': fk}
                
                
