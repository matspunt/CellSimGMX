import datetime

class Lammps_IO:
    @staticmethod
    def convert_xyz_to_lmp(input_file, output_file):
        with open(input_file, 'r') as xyz_file:
            # Read the number of atoms from the first line
            num_atoms = int(xyz_file.readline())
            xyz_file.readline()

            with open(output_file, 'w') as lmp:
                now = datetime.datetime.now()
                header = "GRO file built from {} at {}\n".format(xyz_file, now.strftime("%H:%M:%S"))
                lmp.write(header)
                lmp.write('LAMMPS coordinate file\n\n')
                lmp.write(f'{num_atoms} atoms\n')
                lmp.write('1 atom types\n\n')
                lmp.write('0.0 10.0 xlo xhi\n')
                lmp.write('0.0 10.0 ylo yhi\n')
                lmp.write('0.0 10.0 zlo zhi\n\n')

                lmp.write('Atoms\n\n')
                for i in range(1, num_atoms + 1):
                    line = xyz_file.readline().split()
                    atom_type = 1  #This works with just a single atom type
                    x, y, z = map(float, line[1:4])
                    lmp.write(f'{i} {atom_type} {x} {y} {z}\n')
                
#####unit testing
Lammps_IO().convert_xyz_to_lmp('CELL.xyz', 'CELL.lmp')