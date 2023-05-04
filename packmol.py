class PackmolExecuter:
    """
    Provides an interface for calling PACKMOL through Python.
    
    Methods
    -------
    check_packmol_path()
        Checks if PACKMOL is installed, and prints the path to the binary if it is found.
        
    construct_bead_xyz(bead_identities)
        Based on the bead 'force field', constructs the XYZ files required for PACKMOL
        
    run_packmol(input_file)
        Runs PACKMOL with the constructed input file, or raises an exception if the run fails.
    """

    def __init__(self):
        self.packmol_path = find_executable("packmol")
        
    def check_packmol_path(self):
        if self.packmol_path is not None:
            print(f"Using PACKMOL from: {self.packmol_path}")
        else:
            print(f"Error: cannot locate PACKMOL binary. Are you in the right environment?")
      
    def construct_bead_xyz(self, bead_identities):
        """
        Saves a .XYZ file for each individual bead in the CELL for PACKMOL input

        Args:
            bead_identities (list): List of bead identities.
        """
    
        xyz = "1\n\nX         10.00000       10.00000       10.00000"    
    
        for bead in bead_identities:
            file_name = f"{bead}.xyz"
            with open(file_name, "w") as f:
                f.write(xyz.replace("X", bead))
 
    def run_packmol(self, input_file):
        """
        Runs PACKMOL with the given input file, or raises an exception if the run fails.
        
        Parameters
        ----------
        input_file : str
            The input file for PACKMOL, as a string.
        """
        try:
            self.check_packmol_path()
            # Directly parsing the .inp leads to an error with subprocess
            # see: https://github.com/mosdef-hub/mbuild/issues/419, use tempfile as workaround
            packmol_inp = tempfile.NamedTemporaryFile(mode="w", delete=False, prefix="packmol-", suffix=".inp")
            packmol_inp.write(input_file)
            packmol_inp.close()
            #create the bead name .xyz. THIS IS CURRENTLY HARDCODED FOR THREE BEADS!!!!!
            bead_list = ["A", "B", "N"] # HARDCODED
            self.construct_bead_xyz(bead_list)
            now = datetime.datetime.now()
            log_file = "PACKMOL_build-{}.log".format(now.strftime("%d-%m-%H-%M-%S"))
            with open(log_file, "w") as f:
                proc = sp.run("{} < {}".format(self.packmol_path, packmol_inp.name), stdout=sp.PIPE, universal_newlines=True, shell=True)
                stdout = proc.stdout
                f.write(stdout)

            if 'Solution written to file: CELL.xyz' in proc.stdout:
                print("\n\nThe coordinate .xyz has been successfully built. A logfile is saved at: " + log_file)
                # remove all .xyz files except for CELL.xyz
                for file in glob.glob('*.xyz'):
                    if file != 'CELL.xyz':
                        os.remove(file)
            else:
                error_message = f"\n\nPACKMOL run failed. Last 20 lines of logfile " + log_file + ":\n"
                error_message += "\n".join(proc.stdout.splitlines()[-20:])
                raise Exception(error_message)
        except:
            print("Stopping the run. Please ensure PACKMOL is installed properly as the 'packmol' binary.")
            return
