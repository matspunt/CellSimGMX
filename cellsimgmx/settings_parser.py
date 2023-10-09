import json
import glob
import sys
import os
import argparse

class CLIParser:
    """
    This class accepts command line options of the programme and holds them as child. Can be passed on to other classes to access CLI strings given by user. 
    """
    def __init__(self):
        self.parser = argparse.ArgumentParser(description="CellSimGMX: A 3D discrete element framework simulator for epithelial tissues using GROMACS")
        
        #specifying directories for programme
        self.parser.add_argument("--ff-dir", "-ff", help="Path to the directory where the forcefield is stored (accepts arbitrary .itp file names)")
        self.parser.add_argument("--input-dir", "-in", help="Path to the directory where the simulation input is stored (accepts arbitrary .json file name)")
        self.parser.add_argument("--output-dir", "-out", help="Path where files should be generated and simulation should be run")
        
        #Miscellaneous
        self.parser.add_argument("--verbose", "-v", action="store_true", help="Optional argument. When enabled, prints detailed logging. Useful for debugging problems. ")

    def parse_args(self):
        return self.parser.parse_args()

class JSONParser:
    def __init__(self):
        self.cli_parser = CLIParser()

    def process_json(self):
        args = self.cli_parser.parse_args()
        if args.ff_dir:
            print(f"Forcefield directory: {args.ff_dir}")
        if args.input_dir:
            print(f"Input directory: {args.input_dir}")



                
