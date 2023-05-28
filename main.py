from packmol import PackmolExecuterSingleCell
from gromacs import Gromacs_IO

"""
    Currently (28/05), there are some redundant calls to the parser classes. Let's monitor how much that slows down the execution, if this is insignificant
    I prefer leaving it like this as it makes everything modular (i.e. GMX and PACKMOL classes are totally separated)
    If we want only a single call to the parser functions, I need to rethink how I inherit parsing functionality. Let's see...
"""

input_dir = "/wrk/matspunt/coding/CELL_MODEL" 

packmol = PackmolExecuterSingleCell(input_dir, input_dir)
packmol.run_packmol_single_CELL()

GMX = Gromacs_IO(input_dir, input_dir)
GMX.convert_xyz_to_gro("CELL.xyz", "CELL.gro")
GMX.build_GMX_top_single_CELL("CELL.gro")

