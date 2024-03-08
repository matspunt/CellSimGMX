# -*- coding: utf-8 -*-
# Copyright (C) 2024 Mats Punt mats.punt(at)helsinki.fi
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

import logging
import datetime
from cellsimgmx import JSONParser, ForcefieldParserGMX
from cellsimgmx import CellTopology
from cellsimgmx import SimulationPreparation, ExecuteSimulations

def main():
    ### LOGGING DETAILS
    ### note, logfile is saved where programme is executed, change this to output directory at some point if possible
    now = datetime.datetime.now()
    logging.basicConfig(
    filename = "cellsimgmx-{}.log".format(now.strftime("%H-%M")),
    level=logging.INFO, #print >INFO msgs to logfile
    format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    terminal_handler = logging.StreamHandler()
    terminal_handler.setLevel(logging.WARNING)  # Only print WARNINGS or ERRORS to terminal
    terminal_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    logging.getLogger('').addHandler(terminal_handler)
    
    log_handler = None
    for handler in logging.getLogger('').handlers:
        if isinstance(handler, logging.FileHandler):
            #to obtain log file name
            log_handler = handler
            break
    #### END LOGGING DETAILS
    
    logging.info("Started programme execution. ")
    
    #extract input information from user settings and process the force field
    # Note: can clean up __main__ logic when code is production ready
    json_parser = JSONParser()
    ff_parser = ForcefieldParserGMX()
    json_parser.load_json()
    json_parser.extract_json_values()
    ff_parser.parse_GMX_ff()
     
    celltopology = CellTopology()
    celltopology.build_gro_file_cell()
    celltopology.build_top_from_cell()
    
    mdp_generation = SimulationPreparation()
    execute_sims = ExecuteSimulations()
        
    print(f"CellSimGMX finished. '{log_handler.baseFilename}' has been saved.")
    
if __name__ == "__main__":
    main() 