from packmol import PackmolExecuterSingleCell
from gromacs import Gromacs_IO

input_dir = "/wrk/matspunt/coding/CELL_MODEL/src" 

packmol = PackmolExecuterSingleCell(input_dir, input_dir)
packmol.run_packmol_single_CELL()

GMX = Gromacs_IO(input_dir, input_dir)
GMX.convert_xyz_to_gro("CELL.xyz", "CELL.gro")
GMX.build_GMX_top_single_CELL("CELL.gro")

#Note: to package this as a Python package, see: https://packaging.python.org/en/latest/tutorials/packaging-projects/

