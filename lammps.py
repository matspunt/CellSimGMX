class Lammps_IO:
    @staticmethod
    def convert_xyz_to_lmp(input_file, output_file):
        with open(input_file, 'r') as xyz_file:
            # Read the number of atoms from the first line
            num_atoms = int(xyz_file.readline())

            # Skip the comment line
            xyz_file.readline()

            # Open the output LAMMPS file for writing
            with open(output_file, 'w') as lmp_file:
                # Write the header section
                lmp_file.write('LAMMPS coordinate file\n\n')
                lmp_file.write(f'{num_atoms} atoms\n')
                lmp_file.write('1 atom types\n\n')
                lmp_file.write('0.0 10.0 xlo xhi\n')
                lmp_file.write('0.0 10.0 ylo yhi\n')
                lmp_file.write('0.0 10.0 zlo zhi\n\n')

                # Write the atom properties section
                lmp_file.write('Atoms\n\n')
                for i in range(1, num_atoms + 1):
                    line = xyz_file.readline().split()
                    atom_type = 1  # Assuming a single atom type for simplicity
                    x, y, z = map(float, line[1:4])
                    lmp_file.write(f'{i} {atom_type} {x} {y} {z}\n')
                    
        print(f'Conversion complete. Output file: {output_file}')

#####unit testing
Lammps_IO().convert_xyz_to_lmp('CELL.xyz', 'CELL.lmp')