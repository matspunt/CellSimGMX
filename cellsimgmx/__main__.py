# -*- coding: utf-8 -*-
# Copyright (C) 2023  Mats Punt mats.punt(at)helsinki.fi
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
from cellsimgmx import CellConstructor

def main():
    now = datetime.datetime.now()
    
    ### LOGGING DETAILS
    logging.basicConfig(
    filename = "cellsimgmx-{}.log".format(now.strftime("%H-%M-%S")),
    level=logging.INFO, #print >INFO msgs to logfile
    format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    terminal_handler = logging.StreamHandler()
    terminal_handler.setLevel(logging.WARNING)  # Only print WARNINGS or ERRORS to terminal!
    terminal_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))

    logging.getLogger('').addHandler(terminal_handler)
    #### END LOGGING DETAILS
    
    logging.info("Started programme execution. ")
    
    #extract input information from user settings and process the force field
    json_parser = JSONParser()
    ff_parser = ForcefieldParserGMX()
    json_parser.load_json()
    json_parser.extract_json_values()
    ff_parser.parse_GMX_ff()
    
    # build an individual cell
    cellconstruct = CellConstructor()
    cellconstruct.pack_particles_on_shape()
    cellconstruct.build_tmp_gro_from_object()
    
    logging.info("Succesfully ended programme execution. No errors were found. \n")
    
if __name__ == "__main__":
    main()