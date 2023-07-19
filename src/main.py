from packmol import PackmolExecuterSingleCell
from gromacs import Gromacs_IO
from tissue import TissueConstruction

#we specify an input directory for our forcefield.itp and input.JSON
input_dir = "/wrk/matspunt/coding/CELL_MODEL/src" 

packmol = PackmolExecuterSingleCell(input_dir, input_dir)
packmol.run_packmol_single_CELL()

GMX = Gromacs_IO(input_dir, input_dir)
GMX.convert_xyz_to_gro("CELL.xyz", "CELL.gro", edge_offset=2.5) #edge offset defines the distance to the box edges (the box size is guessed from the coordinates)
GMX.build_GMX_top_single_CELL("CELL.gro") #build 'CELL.top'

print(f"\nNow constructing the tissues from the PACKMOL model:\n")
TissueConstruction.replicate_cell_on_grid('CELL.xyz', 27, 40) #filename, nr of cells, offset in Angstrom
TissueConstruction.replicate_cell_monolayer('CELL.xyz', 4, 40)
TissueConstruction.replicate_sheared_cell_layer('CELL.xyz', 54, 40, shearing=0.5) #shearing factor additional
TissueConstruction.random_packing_in_box('CELL.xyz', tolerance = 8, nr_of_cells = 40, box_x = 150, box_y = 150, box_z = 150) #for random packing, we have to give box dimensions

#Using the Gromacs module, we can convert the .xyz models to .gro
print(f"\nNow converting the tissues to .GRO format:\n")
GMX.convert_xyz_to_gro("CELL_27_standard.xyz", "CELL_27_standard.gro", edge_offset=3.0)
GMX.convert_xyz_to_gro("CELL_4_monolayer.xyz", "CELL_4_monolayer.gro", edge_offset=3.0)
GMX.convert_xyz_to_gro("CELL_54_sheared.xyz", "CELL_54_sheared.gro", edge_offset=3.0)
GMX.convert_xyz_to_gro("CELL_40_random_packing.xyz", "CELL_40_random_packing.gro", edge_offset=3.0)

#Note: to package this as a Python package, see: https://packaging.python.org/en/latest/tutorials/packaging-projects/

