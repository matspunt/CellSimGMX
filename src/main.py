from packmol import PackmolExecuterSingleCell
from gromacs import Gromacs_IO

"""
    Currently (06/06), there are some redundant calls to the parser classes. Let's monitor how much that slows down the execution, if this is insignificant
    I prefer leaving it like this as it makes everything modular (i.e. GMX and PACKMOL classes are totally separated)
    If we want only a single call to the parser functions, I need to rethink how I inherit parsing functionality. Let's see...
    Using 'numba' in code to speedup? I don't really use numpy that much, since 'genfromtxt()' or 'loadtxt()' are quite slow. 
"""

input_dir = "/wrk/matspunt/coding/CELL_MODEL/src" 

packmol = PackmolExecuterSingleCell(input_dir, input_dir)
packmol.run_packmol_single_CELL()

GMX = Gromacs_IO(input_dir, input_dir)
GMX.convert_xyz_to_gro("CELL.xyz", "CELL.gro")
GMX.build_GMX_top_single_CELL("CELL.gro")


#Note: to package this as a Python package, see: https://packaging.python.org/en/latest/tutorials/packaging-projects/

