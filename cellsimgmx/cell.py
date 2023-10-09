from cellsimgmx import CLIParser
from cellsimgmx import JSONParser
from cellsimgmx import ForcefieldParserGMX

class CellConstructor:
    """
    This class builds a .gro file of a single Cell, the basic building block of a System
    """
    def __init__(self):
        self.cli_parser = CLIParser()
        self.json_parser = JSONParser()
        self.forcefieldparsergmx = ForcefieldParserGMX()
        
#class CellTopology:
    